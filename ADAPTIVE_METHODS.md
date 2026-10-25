# Adaptive Measurement-Fusion Reconstruction

**Primary algorithmic contribution of this extension.**

This document describes two *optional* additions to the reconstruction
pipeline of Wu, Gan & Chen, *"Wide Field-of-View MMW SISO-SAR Image
Reconstruction Based on Curved Linear Array"*, IEEE SPL 2025:

1. **Adaptive spatial weighting**
2. **Adaptive HPBW gating**

Both are **off by default**. With

```matlab
cfg.weighting.mode = 'original';   % default
cfg.gating.mode    = 'fixed';      % default
```

`measurementFusionReconstruction.m` is **numerically identical** to the
original paper implementation (verified byte-for-byte by
`validate_project.m`, check *"Original mode reproduces previous output"*,
max relative difference `0.0e+00`).

---

## 1. What belongs to the paper

| Quantity | Paper eq. | Implementation |
|---|---|---|
| Uniform sub-aperture selection | Eq. (4)–(7) | `createUniformSubAperture.m` |
| HPBW gate `δ = 1{α < B_θ/2}` with **fixed** `B_θ = 60°` | Eq. (10),(12) | `calculateHPBWGate.m` |
| Back-projection integral | Eq. (8) | `backProjectionCore.m` |
| Spatial-grid weight `w = 1 / |{samples illuminating the voxel}|` | Eq. (9) | `calculateSpatialWeights.m` |
| Final image `σ̂' = w ⊙ (gated BPA)` | Eq. (11),(13) | `measurementFusionReconstruction.m` |

Everything below is **our contribution** and is **not** claimed to come
from the paper.

---

## 2. Adaptive spatial weighting  (`computeAdaptiveWeight.m`)

### 2.1 Motivation / limitation of the original

The paper's Eq. (9) weight is *geometry only*: it depends solely on how
many uniform sub-aperture samples see a voxel inside the beam. Two
voxels with the same illumination count get the same weight even if one
holds a strong, phase-coherent scatterer and the other holds only
low-amplitude clutter / sidelobe energy. The weight therefore cannot
concentrate the image on the parts that were actually reconstructed with
high confidence.

### 2.2 Formula

Let `rawImg` be the HPBW-gated back-projection (Eq. 8 + Eq. 12) and
`count` the per-voxel illumination count.

```
geometryWeight = 1 ./ max(count, 1)                       % == Eq. (9)
signalNorm     = |rawImg| ./ (max(|rawImg|) + epsilon)    % in [0, ~1]
confidence     = (signalNorm + epsilon) .^ alpha
adaptiveFactor = clip(confidence, minFactor, maxFactor)
adaptiveFactor = adaptiveFactor ./ mean(adaptiveFactor | count>0)   % unit mean
weight         = geometryWeight .* adaptiveFactor
weight(count==0) = 0
```

### 2.3 Parameters (`config.m`, section *ADAPTIVE RECONSTRUCTION EXTENSIONS*)

| Parameter | Default | Meaning |
|---|---|---|
| `cfg.weighting.mode` | `'original'` | `'original'` → exactly Eq. (9); `'adaptive'` → formula above |
| `cfg.weighting.alpha` | `0.5` | confidence exponent. `alpha = 0` ⇒ `confidence ≡ 1` ⇒ reduces to Eq. (9) |
| `cfg.weighting.minFactor` | `0.5` | lower clip — stops real-but-dim scatterers being suppressed |
| `cfg.weighting.maxFactor` | `2.0` | upper clip — **stops weak noisy voxels being boosted without bound** |
| `cfg.weighting.epsilon` | `1e-3` | numerical floor for the normalised-magnitude ratio |

### 2.4 Numerical stability

* `epsilon` floor ⇒ a near-zero voxel gives `confidence ≈ epsilon^alpha`,
  never `0`, `Inf` or `NaN`.
* `clip(minFactor, maxFactor)` bounds the multiplicative correction on
  both sides. Keep `maxFactor` small (≈ 2) so the method cannot turn
  into a noise amplifier.
* The **unit-mean normalisation** over illuminated voxels keeps the
  global image scale comparable to the original weighting, so dB
  displays and energy-fraction metrics remain directly comparable
  between modes.
* Zero-illumination voxels are forced to `weight = 0` (no measurement to
  weight), identical to the original.

`validate_project.m` checks: weights finite, non-negative, no NaN/Inf,
raw confidence factor within `[minFactor, maxFactor]`, and `weight = 0`
wherever `count = 0`.

### 2.5 This is a heuristic

The formula is a reasonable, bounded, stable choice — not a derived
optimum and **not from the paper**. A maximum-likelihood or
iterative-reweighting formulation would be a natural next step. The
`maxFactor` clip plus unit-mean normalisation were chosen specifically
so the method degrades gracefully rather than amplifying noise; the SNR
experiment (`experimentNoiseRobustness.m`) tests exactly that.

---

## 3. Adaptive HPBW gating  (`computeAdaptiveHPBW.m`)

### 3.1 Motivation

The paper uses one fixed gate half-angle for every sample. On the curved
array the **central elevation rings** are heavily over-sampled (they
produce the central circular clutter the paper describes) while the
**peripheral rings** (`φ → φ_c`) are weakly sampled and lose coverage.

### 3.2 Formula (geometry only — no circular logic)

Per uniform-sub-aperture sample `s` with elevation `φ(s) ∈ [0, φ_c]`:

```
normElev(s)     = φ(s) / max(φ)                                  in [0, 1]
effHPBW_deg(s)  = HPBW_deg + adaptStrength * HPBW_deg * (normElev(s) - 0.5)
effHPBW_deg(s)  = clip(effHPBW_deg(s), minHPBW_deg, maxHPBW_deg)
```

Central rings (`normElev → 0`) gate **tighter** (reject more off-axis
clutter); peripheral rings (`normElev → 1`) gate **wider** (retain more
of their sparse contributions). Every value is a physically meaningful
beamwidth in degrees, hard-bounded by the configured limits.

The **same** per-sample vector feeds both `calculateSpatialWeights`
(so `count` stays consistent) and the gated back-projection
(`backProjectionCore(..., gateHPBW_deg)`).

Because `effHPBW_deg` depends only on the sampling geometry, a single
reconstruction pass is sufficient — there is no dependence on the
reconstructed signal, hence no circularity. (A two-pass
confidence-refined gate was considered and deliberately **not** adopted,
to keep runtime and interpretability simple.)

### 3.3 Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `cfg.gating.mode` | `'fixed'` | `'fixed'` → nominal `cfg.HPBW_deg` gate (paper); `'adaptive'` → formula above |
| `cfg.gating.minHPBW_deg` | `45` | hard lower bound on the effective per-sample HPBW |
| `cfg.gating.maxHPBW_deg` | `75` | hard upper bound on the effective per-sample HPBW |
| `cfg.gating.adaptStrength` | `0.30` | fraction of nominal HPBW swung between central & peripheral rings. `0` ⇒ reduces to `'fixed'` |

With the defaults the effective HPBW spans ≈ `[51°, 69°]` (verified).

`validate_project.m` checks: `effHPBW ∈ [minHPBW_deg, maxHPBW_deg]`, and
that a constant per-sample HPBW vector produces a gate **bit-identical**
to the original scalar gate.

---

## 4. Supported mode combinations

| Mode | `cfg.weighting.mode` | `cfg.gating.mode` | Note |
|---|---|---|---|
| **A** | `original` | `fixed` | the paper — **default** |
| **B** | `adaptive` | `fixed` | adaptive weighting only |
| **C** | `original` | `adaptive` | adaptive gating only |
| **D** | `adaptive` | `adaptive` | both |

All four are handled by config parameters inside the single
`measurementFusionReconstruction.m` — no duplicate functions.

---

## 5. Diagnostics

`[img, diag] = measurementFusionReconstruction(...)` returns, in
addition to the original `diag.subAp / weight / count / rawGatedImage`:

| Field | Contents |
|---|---|
| `diag.geometryWeight` | the Eq. (9) weight (always) |
| `diag.signalConfidence` | `signalNorm` (adaptive weighting only; `[]` otherwise) |
| `diag.adaptiveFactor` | bounded, unit-mean confidence factor (`[]` otherwise) |
| `diag.effectiveHPBW` | **summary** struct: `mode`, `nominal_deg`, `min_deg`, `max_deg`, `mean_deg` (+ `perSample_deg` only if ≤ 20000 samples) |
| `diag.weightingMode`, `diag.gatingMode` | the resolved modes |
| `diag.runtimeSeconds`, `diag.executionMode`, `diag.chunkSizeUsed` | see PERFORMANCE_OPTIMIZATION.md |

Large per-voxel/per-sample arrays are **not** duplicated into `diag`
beyond what the original already stored.

---

## 6. Experiments

| Script | What it does |
|---|---|
| `experiments/experimentAdaptiveMethods.m` | Modes {Classical BPA, Original Measurement Fusion, Adaptive Weighting, Adaptive HPBW, Adaptive Weighting + Adaptive HPBW} on: (A) planar point array — full metric table, (B) resolution vs elevation angle, (C) human point cloud. |
| `experiments/experimentNoiseRobustness.m` | SNR sweep `[Inf, 30, 20, 10] dB` (complex AWGN) comparing Original Measurement Fusion / Adaptive Weighting / Adaptive Weighting + Adaptive HPBW — shows whether adaptive weighting amplifies noise. |

Metric names used verbatim in plots:
`Classical BPA`, `Original Measurement Fusion`, `Adaptive Weighting`,
`Adaptive HPBW`, `Adaptive Weighting + Adaptive HPBW`.

Every metric not in the paper is tagged
**"Additional metric - not reported in original paper"** in the console
output.

### Quick-mode results (indicative, `config('quick')`, planar point array)

| Method | artifact energy frac. | peak CoV | contrast dB |
|---|---|---|---|
| Classical BPA | 0.228 | 0.134 | 7.83 |
| Original Measurement Fusion | 0.197 | 0.094 | 7.70 |
| Adaptive Weighting | **0.080** | 0.138 | **10.81** |
| Adaptive HPBW | 0.207 | 0.103 | 7.50 |
| Adaptive Weighting + Adaptive HPBW | 0.086 | 0.145 | 10.61 |

**Honest reading:** adaptive weighting roughly halves the
out-of-support (artifact) energy and adds ≈ 3 dB of contrast, but on
this coarse quick-mode grid it slightly *worsens* the peak-amplitude
coefficient of variation (0.094 → 0.138). Adaptive HPBW alone is close
to neutral here. These trade-offs are printed and saved unchanged by the
experiment — nothing favourable-only is reported. Re-run with
`config('full')` for paper-scale sampling.

---

## 7. Limitations

* The adaptive weight formula is heuristic (§2.5).
* The adaptive gate uses elevation only; azimuth-dependent sampling
  density is not yet used.
* Quick mode is deliberately under-sampled; some adaptive benefits only
  emerge at `config('full')` density (same caveat the original project
  documents for the point-array metrics).
* No claim is made that these methods appear in, or are endorsed by, the
  original paper.
