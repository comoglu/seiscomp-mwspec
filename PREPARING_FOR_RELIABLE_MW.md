# Preparing for a reliable Mw(spec)

Spectral Mw is only as good as the inputs you give it. This is the checklist of
what to make ready, in order of impact, with the exact configuration keys. (The
algorithm is a faithful port of Seisan SPEC/AUTOMAG; Seisan gets clean results
because operators prepare exactly these things — a calibrated Q, curated
stations, and per‑station corrections — not because of a different method.)

All "shared model" keys live in **module config** (e.g. `global.cfg`, under
`Modules → global` in scconfig). Measurement/moment/QC keys live in the
**Amplitudes / Magnitudes bindings** for `Mw(spec)`.

---

## 1. Attenuation Q — the single biggest lever  ⚠ required for reliable Mw
Without an attenuation correction the spectrum is over‑attenuated at high
frequency and Mw is distance/magnitude‑biased. Calibrate Q to your region (the
analogue of Seisan **CODAQ / QLG**, or use a published regional Q) and set:

```
magnitudes.Mw(spec).Q0     = 450      # Q(f) = Q0 * f^Qalpha   (region-dependent!)
magnitudes.Mw(spec).Qalpha = 0.7
```
The default `Q0 = 0` means **no** anelastic correction — adequate only at very
short distance. For a depth‑dependent crust use the full layered `model` instead
(see §2), which carries Q per layer.

## 2. Source velocity & density
Used in M0 = 4π·ρ·c³·R·Ω₀ / (radiation·freeSurface). Use regional crustal values:
```
magnitudes.Mw(spec).vp      = 6.1     # km/s   (P; use vs for S phase)
magnitudes.Mw(spec).vs      = 3.5
magnitudes.Mw(spec).density = 2.7     # g/cm^3
```
Or, for depth dependence, the layered model (overrides the scalars above):
```
magnitudes.Mw(spec).model = "3 5.8 3.2 450 0.7 400 0.7 2.6", "15 6.8 3.9 450 0.7 400 0.7 2.9", ...
#   depth vp vs Qp0 Qp_alpha Qs0 Qs_alpha density   (linearly interpolated to source depth)
magnitudes.Mw(spec).phase = P         # P (vertical) or S (horizontal)
```

## 3. Station selection — exclude what spectral Mw can't handle
- **Drop teleseismic stations** (the simple geometric spreading is wrong beyond
  ~10–15°). Cap the distance on the magnitude side:
  ```
  magnitudes.Mw(spec).minimumDistance = 0
  magnitudes.Mw(spec).maximumDistance = 12     # degrees
  ```
- **Don't use low‑quality sensors** (e.g. RaspberryShake) for Mw — exclude them
  from the bindings.
- Prefer **regional broadband** with good low‑frequency response (you need signal
  below the corner frequency to see the Ω₀ plateau).

## 4. Per‑station corrections (optional) — remove residual site bias
The SeisComP base applies `corrected = multiplier·Mw + offset` per station. Set
it in a station's binding:
```
magnitudes.Mw(spec).offset     = -0.20   # = -(that station's site term)
magnitudes.Mw(spec).multiplier = 1.0
```
Generate the offsets from a calibration with the companion tool:
```
python3 station_offsets.py report_clean/station_corrections.csv --min-obs 8 --out offsets.csv
```
(Positive site term ⇒ station reads high ⇒ negative offset.)

## 5. Quality gates & averaging
Reject bad measurements and average robustly:
```
amplitudes.Mw(spec).minSNR      = 2.0    # spectral S/N
amplitudes.Mw(spec).maxResidual = 1.0    # Brune-fit misfit
# scmag: take the median over stations and require a minimum count
magnitudes.average            = median
magnitudes.minimumMagnitudes  = 5        # don't trust an event with too few stations
```

## 6. Window & band
```
amplitudes.Mw(spec).signalPreTime  = 1.0
amplitudes.Mw(spec).signalDuration = 25.0   # long enough to resolve low fc (large events)
amplitudes.Mw(spec).fmin = 0.0              # 0 = automatic S/N band selection
amplitudes.Mw(spec).fmax = 0.0
```

## 7. Absolute calibration (one-time)
After the above, check the median (Mw(spec) − reference Mw) on events that have a
**well‑constrained** moment tensor and, if needed, apply a constant:
```
amplitudes.Mw(spec).calibration = +0.10     # log10 additive on Omega0 (~ +0.07 Mw per +0.1)
```

---

## Reliable operating range (typical, after the above)
- **Magnitude:** ~M3 to ~M6. Below ~M3 the Ω₀ plateau is often below noise;
  above ~M6 the corner frequency drops below the usable band (saturation).
- **Distance:** regional (≲ ~12°) unless you extend the spreading/Q model.
- **Stations:** ≥ ~5 good station magnitudes for a stable network Mw.

## TL;DR — what you must prepare
1. a **regional Q(f)** (Q0, Qalpha) — most important;
2. **regional velocity/density** (or a layered model);
3. a **curated station set** (regional broadband; no teleseismic, no toy sensors);
4. optionally **per‑station offsets** from a calibration;
5. **QC + median** and a sensible distance/magnitude range.
