# mwspec — spectral moment magnitude `Mw(spec)` for SeisComP

A native SeisComP magnitude plugin that computes a **per-station moment
magnitude** from a Brune ω²-fit of the displacement spectrum. It is a faithful
port of the Seisan **SPEC / AUTOMAG** spectral-magnitude algorithm
(L. Ottemöller) and fills the one magnitude type Seisan provides that SeisComP
did not: spectral `Mw` computed independently at each station (rather than
scaled from mB/Mwp or derived from a full moment-tensor inversion).

Registers an amplitude processor and a magnitude processor of type
**`Mw(spec)`** with `scamp` / `scmag` / `scolv`.

> **Status:** v0.4.1. P (vertical) and S (N+E vector-sum) phases; configurable
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
- Important implementation detail: SeisComP's `deconvolveFFT` removes only the
  normalised response *shape* — the processor divides out the sensitivity (gain)
  itself, as `ML`/`MN`/`A5_2` do.

## Credit / license
Ports the Seisan **SPEC / AUTOMAG** spectral-Mw algorithm by Lars Ottemöller.
Distributed under the **GNU Affero General Public License v3.0** — see
[`LICENSE`](LICENSE) and the per-file source headers. SeisComP® is a trademark
of gempa GmbH / GFZ; this is an independent, unaffiliated plugin.
