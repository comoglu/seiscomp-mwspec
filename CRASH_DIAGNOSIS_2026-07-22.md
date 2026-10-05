# scamp SIGSEGV in Mw(spec) combiner — diagnosis (2026-07-22)

## Update 2026-10-06 — most likely mechanism, and what v0.6.0 changes

Not reproduced (the playback data and the apport dump are gone), but the code
explains the backtrace:

- `Client::StreamApplication::readRecords()` (streamapplication.cpp ~352) calls
  `storeRecord(rec)` inside a `try`, and its `catch (...)` does `delete rec`
  then logs "Skipping invalid record …" (reading `rec` after the delete).
  If **anything inside `storeRecord()` throws after a processor stored the
  record as `_stream.lastRecord`**, the record is freed while still
  referenced. The next `store()` on that processor releases it → `0x0` frame
  in `decrementReferenceCount`, exactly the backtrace above.
- So the "Skipping invalid record (fsamp 10, nsamp 163…)" lines just before
  the crash are not bad records rejected up front: they are this catch firing
  for ordinary records. Mw(spec) is simply the processor that keeps
  `lastRecord` longest (25 s window), so it trips first. Which code threw is
  still unknown — it can be any processor fed that record (Mwpd, stock, us).
- SeisComP core fix (upstream): don't `delete` a record that `storeRecord()`
  has taken over (validate `endTime()` in a separate try), and log before
  deleting.

v0.6.0 plugin changes:
- `AmplitudeProcessor_MwSpecCombiner::feed()` never lets an exception escape
  (logs it, sets status Error), so Mw(spec) can no longer trigger the delete.
- The worker no longer dereferences `_environment.hypocenter` in
  `computeAmplitude()`: scamp does not keep a messaging Origin alive, so in
  real time that pointer could be dangling ~25 s later. Depth, origin time and
  hypocentral distance are copied in `setEnvironment()`.

Remaining exposure: another processor throwing for a record that Mw(spec)
already holds would still crash in Mw(spec)'s next `store()` until the core
is fixed. To find the thrower, re-run the playback with all processors and
check for "Skipping invalid record" lines on streams with valid data.

**Author:** captured while debugging why Mwpd wasn't appearing in a SeisComP 8.0.0
playback. Root cause turned out to be **this plugin**, not Mwpd.

## One-line summary

`scamp` segfaults (SIGSEGV) inside `AmplitudeProcessor_MwSpecCombiner::feed()`
when it forwards a record to a worker (`_c0.feed(record)`), because the worker's
`_stream.lastRecord` (a `boost::intrusive_ptr<const Record>`) is **dangling** —
the base `WaveformProcessor::store()` crashes while releasing the previously held
record. Only the fast magnitudes (mb / ML / MLv) survive because they finalize
before the crash; Mwpd, Mwp and Mw(spec) itself need a longer window and are lost
when the process dies.

## Environment

- SeisComP **8.0.0 Development**, API 18.0.0, GIT HEAD `9a2928e3`, `/opt/seiscomp`, Ubuntu 24.04.
- Installed plugin built from **this repo** (`/home/ubuntu/Projects/seiscomp-mwspec`), branch `develop` (`a57b8ac`, v0.5.0).
- Crash dump: `/var/crash/_opt_seiscomp_bin_scamp.1000.crash` (apport, Signal 11), 2026-07-22 20:06:44 local.
- Reproduced on a real-time playback (`msrtsimul2 -j 7`) of GA event `ga2026fvuvzv`. Also seen on earlier plugin versions per operator — this is **not** new.

## Backtrace (top frames)

```
#0  0x0000000000000000 in ?? ()
#1  Seiscomp::Core::BaseObject::decrementReferenceCount ()      baseobject.inl:51
#2  Seiscomp::Core::intrusive_ptr_release ()                    baseobject.inl:32
#3  boost::intrusive_ptr<Seiscomp::Record const>::~intrusive_ptr ()
#4  boost::intrusive_ptr<Seiscomp::Record const>::operator= ()  intrusive_ptr.hpp:162
#5  Seiscomp::Processing::WaveformProcessor::store ()           waveformprocessor.cpp:375   <-- "_stream.lastRecord = record;"
#6  MwSpec::AmplitudeProcessor_MwSpecCombiner::feed ()          combiner.cpp:255            <-- "_c0.feed(record);"
#7  AmpTool::handleRecord ()                                    amptool.cpp:1534   streamID = "S1.AUNHS..BHZ"
```

- Frame #0 is `0x0` (jumped to a null/garbage vtable/dtor) — classic **use-after-free / heap corruption** signature, not a plain null deref of our own code.
- The faulting operation is the assignment `_stream.lastRecord = record;` in the base `store()`. `operator=` first releases the *old* `lastRecord`; that release walks into freed memory → the process dies decrementing a refcount on an already-destroyed `Record`.
- It is reached only via the Mw(spec) combiner forwarding to `_c0` (`combiner.cpp:255`).

## Why Mwpd is not involved

Mwpd (the duration-amplitude plugin) computed 35 station amplitudes cleanly in
the same run, right up to the instant of the crash. The fault frames are entirely
`MwSpecCombiner` + base `WaveformProcessor`. The "no Mwpd" symptom is a
*side effect*: scamp dies ~26 s into processing, before Mwpd's T0 resolves.

## Context worth noting

Immediately before the crash the log is full of malformed records on the same
network being skipped:

```
[error/StreamApplication] Skipping invalid record for S1.AUKAT..BHZ (fsamp: 10.00, nsamp: 163)
[error/StreamApplication] Skipping invalid record for S1.AULHS..BHZ (fsamp: 10.00, nsamp: 200)
[error/StreamApplication] Skipping invalid record for S1.AUMTS..BHZ (fsamp: 10.00, nsamp: 207)
```

These `S1.AU*` BHZ records claim 10 Hz with odd/short sample counts. The crash
stream `S1.AUNHS..BHZ` is the same family. Suspect that a short/odd/gappy record
(or the sequence of them) drives the combiner/worker into a bad record-lifetime
state. Heap corruption from mishandling such a record is consistent with the
`0x0` top frame.

## Architecture (for the fix)

`AmplitudeProcessor_MwSpecCombiner` holds two **value-member** workers `_c0`,
`_c1` (not pointers), each a full `AmplitudeProcessor`. In the ctor
(`combiner.cpp:51+`) each gets `setPublishFunction(bind(&…::newAmplitude, this, …))`
and then `reset()`. `feed()` (line ~231) routes by component:
- `_nActive == 1` (vertical, the Mw(spec) case): forwards Z records to `_c0.feed(record)` — **this is the crashing path (line 255)**.
- otherwise routes the two horizontals to `_c0` / `_c1`.

`_c0.feed()` is the base `WaveformProcessor::feed → store`, and `store()` ends with
`_stream.lastRecord = record;`.

## Hypotheses (ranked)

1. **Worker `_stream.lastRecord` holds a raw/freed Record across a reset or reuse.**
   If `_c0` is reused for a new origin/stream (combiner reset, or the same
   combiner instance re-driven) without its `_stream.lastRecord` being cleared,
   and the previously referenced `Record` was already destroyed, the next
   `store()` assignment releases a dangling pointer. Check `reset()`
   (`combiner.cpp:117`) — it calls `_c0.reset()` / `_c1.reset()`; confirm the base
   `reset()` actually clears `_stream.lastRecord`, and that the combiner isn't
   feeding after a partial/failed reset.

2. **Copy/move of the combiner (or a worker) breaks the intrusive_ptr.** Value-member
   workers with bound `this` publish-functions are fragile: if the combiner object
   is ever copied/moved (e.g. stored by value in a container, or the factory
   returns by value), `_c0`'s bound `this` and its `_stream` refcounts can be left
   inconsistent. Verify the combiner is only ever heap-allocated via the factory
   and never copied.

3. **Double-feed / feeding after Finished.** `feed()` guards on `status() >
   Finished` and `_c0.isFinished()`, but if a worker transitions to a terminal
   state mid-record and is fed again (or fed a record whose buffer is already
   released by AmpTool), the stale `lastRecord` release can fire. Note AmpTool
   holds the record via `tmp` (`intrusive_ptr`) only for the duration of
   `handleRecord`; nothing keeps it alive afterwards.

4. **Malformed-record path.** The short `S1.AU*` records may produce a `store()`
   where `arr`/`_stream` is in an unexpected state. Even if StreamApplication
   "skips" some, others of the family reach the processor.

## How to continue

- **Reproduce deterministically:** replay `ga2026fvuvzv` (the sorted mseed used) with only `Mw(spec)` enabled in `scamp.cfg amplitudes` and the `S1.AU*` streams bound. It crashes within ~30 s.
- **Run under a sanitizer** for the exact offending access:
  ```
  # build the plugin with -fsanitize=address -g and run scamp under it
  seiscomp exec env LD_PRELOAD=libasan.so.8 ASAN_OPTIONS=halt_on_error=1 scamp --plugins mwspec --debug -I <mseed>
  ```
  ASan will name the freed `Record` and both the free and the use stacks.
- **Cheap instrumentation:** log `record` and `_c0._stream.lastRecord.get()` pointers at `combiner.cpp:255` before `_c0.feed()`, and in the worker's `store()` around line 375, to see the dangling value.
- **Likely fix location:** ensuring `_c0`/`_c1` clear `_stream.lastRecord` on `reset()` and are never fed a record whose lifetime isn't guaranteed for the store; and confirming the combiner instance is never copied.

## Operational workaround (already applied)

Remove `Mw(spec)` from `amplitudes = …` in `scamp.cfg`, `update-config`, restart
scamp. scamp then runs stably and Mwpd (and everything else) computes normally.
Re-enable Mw(spec) once the record-lifetime bug is fixed.

## Cross-refs

- Base crash line: `common/libs/seiscomp/processing/waveformprocessor.cpp:375` (`_stream.lastRecord = record;`).
- Trigger: `combiner.cpp:255` (`_c0.feed(record);`).
- Combiner reset: `combiner.cpp:117`. Combiner ctor / worker publish binding: `combiner.cpp:51-61`.
