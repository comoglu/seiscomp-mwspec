# mwspec — spectral moment magnitude `Mw(spec)` for SeisComP

A native SeisComP magnitude plugin that computes a **per-station moment
magnitude** from a Brune ω²-fit of the displacement spectrum. It is a faithful
port of the Seisan **SPEC / AUTOMAG** spectral-magnitude algorithm
(L. Ottemöller) and fills the one magnitude type Seisan provides that SeisComP
did not: spectral `Mw` computed independently at each station (rather than
scaled from mB/Mwp or derived from a full moment-tensor inversion).

Registers an amplitude processor and a magnitude processor of type
**`Mw(spec)`** with `scamp` / `scmag` / `scolv`.

> **Status:** v0.6.0. P (vertical) and S (N+E vector-sum) phases; configurable
> velocity/Q/density model, distance gate, per-station corrections and
> calibration. The displacement spectrum is **bit-validated against Seisan's own
> `spectrum()` routine** (×1.000 over a population of real channels), and the
> Brune fit sits inside the SCEC/USGS Community Stress-Drop ensemble.

## Method (summary)
For each pick the phase window is deconvolved to ground **displacement**, the
displacement amplitude spectrum is corrected for path attenuation
(`Q(f)=Q0·f^α` and near-surface `κ`), a usable S/N band is selected, and the
Brune model `log10A = log10Ω0 − log10(1+(f/fc)²)` is fitted by grid search:

```
M0 = 4π·ρ·c³·R·Ω0 / (Fs·Rθφ)          Mw = ⅔·log10(M0) − 6.06
```

with the layered velocity/Q/density model interpolated to the source depth and
the geometric-spreading distance R from Seisan `spec_dist` (hypocentral for P;
Herrmann–Kijko for regional S). See `brune.{h,cpp}` for the (framework-
independent, unit-tested) science core.

## Build & install (in-tree drop-in)
This builds as a normal in-tree SeisComP plugin: drop it into a SeisComP
**source** tree, register it in the parent `CMakeLists.txt`, and rebuild.

```bash
# 1. clone into the magnitudes plugin directory of a SeisComP source checkout
cd <seiscomp-src>/src/base/main/plugins/magnitudes
git clone https://github.com/comoglu/seiscomp-mwspec.git mwspec

# 2. register the subdirectory in the parent CMakeLists.txt
#    add a line:   SUBDIRS(mwspec)
$EDITOR <seiscomp-src>/src/base/main/plugins/magnitudes/CMakeLists.txt

# 3. configure + build + install from your SeisComP build directory
cd <seiscomp-build> && cmake . && make mwspec && make install   # -> share/plugins/mwspec.so
ctest -R test_mwspec_brune                                      # unit tests
```
The directory **must** be named `mwspec` inside the source tree. `make install`
places `mwspec.so` in `share/plugins/` and the scconfig docs in
`etc/descriptions/`.

## Configure
```
plugins = ${plugins}, mwspec
amplitudes = ${amplitudes}, Mw(spec)         # scamp
magnitudes = ${magnitudes}, Mw(spec)         # scmag
# shared velocity/Q/density model (Modules -> global in scconfig):
magnitudes.Mw(spec).model = "3.0 5.8 3.2 500 0.7 400 0.7 2.6", \
                            "15.0 6.8 3.9 500 0.7 400 0.7 2.9"
magnitudes.Mw(spec).phase = P
```
Per-station measurement/moment parameters are binding profiles (scconfig →
Bindings → Amplitudes/Magnitudes → `Mw(spec)`); see `descriptions/global_mwspec.xml`.

`phase` must be set in the module configuration (as above), not per binding:
scamp/scolv choose the streams (Z, or N+E) before a binding is read. A binding
phase that disagrees with it makes the processor's setup fail with an error.

## Fit diagnostics (comments)
The values behind each measurement are attached as comments, so they reach the
database/QuakeML and can be shown in scolv or by scripts. Comment ids:

| Object | id | Unit | Meaning |
|---|---|---|---|
| Amplitude | `Om0` | nm·s | spectral flat level (= amplitude value for P) |
| Amplitude | `fc` | Hz | corner frequency used for the magnitude |
| Amplitude | `fmin`, `fmax` | Hz | frequency band that was fitted |
| Amplitude | `fitResidual` | — | Brune-fit misfit (gate: `maxResidual`) |
| Amplitude | `deltaKappa` | s | fitted Δκ (only when the dkappa search is on) |
| Amplitude | `travelTime` | s | travel time used for the Q correction (Q mode only) |
| StationMagnitude | `M0` | N·m | seismic moment |
| StationMagnitude | `fc` | Hz | corner frequency |
| StationMagnitude | `sourceRadius` | m | Brune radius 0.37·c/fc |
| StationMagnitude | `stressDrop` | MPa | Brune stress drop 7/16·M0/r³ |
| StationMagnitude | `geoDistance` | km | geometric-spreading distance R (Q mode only) |

For S the per-component values carry a `.N` / `.E` suffix (`Om0.N`, `fc.E`, …);
the unsuffixed `fc` is the combined value the magnitude uses. The amplitude
`methodID` is `Brune/<phase>/<Q|table>`, and `creationInfo.version` is the
plugin version.

## Preparing for reliable Mw
Spectral Mw needs the right inputs (a regionally-calibrated Q, curated stations,
optional per-station corrections). See **`PREPARING_FOR_RELIABLE_MW.md`** for the
full checklist and the exact config keys (`Q0`/`Qalpha`/`vp`/`vs`/`density`,
`minimumDistance`/`maximumDistance`, `multiplier`/`offset`, QC gates).

## Status & notes
- Production-tested on real regional events (gives sensible station magnitudes).
- **Calibration**: the absolute level depends on FFT/taper/deconvolution
  conventions; a one-time offset can be applied via `amplitudes.Mw(spec).calibration`
  (log10 additive) against a reference (Seisan / GCMT / GA).
- **S-wave** combines both horizontal components. Each horizontal is fitted
  independently and the two Omega0 are combined per `amplitudes.Mw(spec).combiner`
  (default `vector_sum` = sqrt(N²+E²), the total horizontal S motion). P uses the
  vertical only. The registered processor is a component combiner that runs one
  worker (P) or two (S); see `combiner.cpp`.
- **S window**: amplitudes are triggered on the P pick, so the S signal window
  (`signalPreTime`/`signalDuration`) currently starts at the P onset, and the Q
  correction uses the P travel time. Choose `signalDuration` to cover the S
  arrival at your distances.
- Important implementation detail: SeisComP's `deconvolveFFT` removes only the
  normalised response *shape* — the processor divides out the sensitivity (gain)
  itself, as `ML`/`MN`/`A5_2` do.

## Author / credit / license
Author: **Mustafa Comoglu** (Geoscience Australia).

Ports the Seisan **SPEC / AUTOMAG** spectral-Mw algorithm by Lars Ottemöller.
Distributed under the **GNU Affero General Public License v3.0** — see
[`LICENSE`](LICENSE) and the per-file source headers. SeisComP® is a trademark
of gempa GmbH / GFZ; this is an independent, unaffiliated plugin.
