# Runtime Optimization & Benchmarking

**Secondary (implementation) contribution of this extension.**

> **What is NOT claimed as novel:** FFT-accelerated back-projection,
> aperture chunking, `parfor` execution and `gpuArray` execution
> **already existed** in `reconstruction/backProjectionCore.m` before
> this extension. We do **not** claim to have invented FFT acceleration,
> chunking, or GPU computing, and none of these is a new radar
> algorithm.
>
> **What IS our contribution:** the unified execution-mode selector, the
> memory-aware *automatic* chunk sizing, the automatic
> hardware/toolbox-aware mode selection (`fft_auto`), the execution
> diagnostics, and the quantitative, reproducible benchmark
> (`experiments/experimentRuntimeBenchmark.m`).

---

## 1. Direct vs FFT back-projection

`backProjectionCore.m` evaluates the discrete Eq. (8) sum

```
img(voxel) = Σ_samples Σ_freq  s(sample,f) · exp(+2j·k(f)·RR(voxel,sample))
```

**Direct** (`cfg.recon.method = 'direct'`): literal double sum,
`O(Nv · Ns · Nf)`.

**FFT** (`'fft'`, default): for a fixed sample the inner sum over `k` is
the inverse DFT of that sample's spectrum evaluated at delay
`τ = 2·RR/c` (substitute `k = 2πf/c` into `exp(2j·k·RR)`). Each sample's
spectrum is IFFT'd **once** (zero-padded by `cfg.recon.zeroPadFactor`
for a fine range grid) and the range profile is linearly interpolated at
`τ(voxel)` per voxel. Cost `O(Ns · (Nfft·log Nfft + Nv))`.

The two are mathematically equivalent up to the constant FFT scale
factor `Nfft` (`MATLAB ifft` normalisation). Both `validate_project.m`
(check *"FFT vs direct Eq.(8) equivalence"*, tolerance `< 0.05` relative
with `zeroPadFactor = 32`) and the benchmark verify this numerically;
the benchmark peak-normalises both images before comparing so the
`Nfft` factor does not appear as an error.

### Complexity summary

| Method | Time | Extra memory |
|---|---|---|
| direct | `Nv·Ns·Nf` | `O(Nv)` |
| fft | `Ns·(Nfft·log Nfft + Nv)`, `Nfft = zeroPadFactor·Nf` | `O(chunk·Nfft + Nv)` |

Indicative (`config('quick')`, `Ns=5760, Nf=301, Nv=343`):
direct ≈ 6.9 s, fft ≈ 0.30 s → **≈ 22× speed-up**, peak-normalised
relative image error ≈ `1.5e-2` (interpolation error at
`zeroPadFactor = 4`; drops well below `1e-3` at `zeroPadFactor ≥ 16`).

---

## 2. Chunking

The aperture is processed in blocks of `cfg.recon.chunkSize` samples to
bound peak memory (the per-chunk IFFT array `G` is
`chunkSize × Nfft` complex). Chunking **never changes the result** — the
coherent sum is associative — only memory and speed.

### 2.1 Fixed chunking (original)

`cfg.performance.chunkMode = 'fixed'` (default). Uses
`cfg.performance.chunkSize` if `> 0`, otherwise the original
`cfg.recon.chunkSize` (3000 full / 1500 quick).

### 2.2 Automatic chunking (our contribution — `estimateChunkSize.m`)

`cfg.performance.chunkMode = 'auto'`. Picks the largest chunk whose
working set fits a memory budget:

```
16 · (chunkSize·Nfft + 8·Nv) · 1.5   ≤   targetBytes
```

* `targetBytes` = `cfg.performance.targetMemoryMB · 1e6 / 1.5`
  (512 MB full / 128 MB quick by default).
* On Windows, `memory` is queried and the target is additionally capped
  at 25 % of `MemAvailableAllArrays`. On platforms without `memory`
  the query is skipped silently and the configured budget is used
  (it always exists, so there is a safe fallback).
* Result clamped to `[1, Ns]`. Fewer, larger chunks ⇒ less loop
  overhead; on the quick point-target problem auto picks `chunkSize ≈
  4427` vs the fixed `1500`.

`validate_project.m` checks: auto chunk size is a valid integer in
`[1, Ns]`, and the reconstruction with auto chunking matches fixed
chunking to `< 1e-10` relative.

---

## 3. Execution modes  (`applyPerformanceMode.m`)

`cfg.performance.mode` is a convenience layer that **resolves onto** the
pre-existing low-level knobs `cfg.recon.method`, `cfg.useGPU`,
`cfg.useParallel`:

| `cfg.performance.mode` | `recon.method` | `useParallel` | `useGPU` |
|---|---|---|---|
| `'direct'` | `direct` | false | false |
| `'fft_serial'` **(default)** | `fft` | false | false |
| `'fft_parallel'` | `fft` | **true** | false |
| `'fft_gpu'` | `fft` | false | **true** |
| `'fft_auto'` | resolved at runtime → GPU if present, else parallel if PCT present, else serial |

The default `'fft_serial'` resolves to *exactly* the original defaults,
so all existing numerical results are reproduced. If `cfg` has **no**
`.performance` field at all (e.g. a hand-built config), the resolver is
a complete no-op.

`fft_auto` reasoning is reported via `info.autoReason` /
`cfg.performance.verbose`.

---

## 4. Parallel CPU

`cfg.performance.mode = 'fft_parallel'` (needs the Parallel Computing
Toolbox — *installed*, not merely licensed).

* Chunks are distributed across a `parfor` loop.
* **Each worker writes only its own column** of `imgParts` /
  `countParts` (local partial accumulators). Partial images are summed
  after the loop → **no race condition**.
* The pool is **not** opened/closed per call — `parfor` reuses the
  existing pool (`gcp('nocreate')` is only *queried*, never forced).
* If the toolbox is unavailable the code falls back to the serial loop
  with a **single** one-time warning (`backProjectionCore:noParallel`).

`validate_project.m` check *"Parallel CPU vs serial numerically close"*
(`< 1e-9` relative) — **skipped and not counted** if the toolbox is
absent (as on the development machine).

---

## 5. GPU

`cfg.performance.mode = 'fft_gpu'` (needs the Parallel Computing Toolbox
and a supported GPU).

* Each chunk's heavy vectorised arithmetic (`RR`, `τ`, gate angle,
  range-profile interpolation) runs on `gpuArray`s.
* Data is moved host→device **once per chunk**, not inside the inner
  per-sample loop; only the per-chunk partial image/count is
  `gather`ed back.
* If no supported GPU is found: **one** warning
  (`backProjectionCore:noGPU`), then graceful CPU fallback — never a
  crash.

`validate_project.m` check *"GPU vs serial numerically close"*
(`< 1e-5` relative) — **skipped and not counted** if no GPU.

---

## 6. Execution diagnostics

`[img, count, info] = backProjectionCore(...)` (3rd output, backward
compatible) and `diag` from `measurementFusionReconstruction` expose:

| Field | Meaning |
|---|---|
| `info.requestedMode` | `cfg.performance.mode` as asked |
| `info.resolvedMode` | after `fft_auto` resolution |
| `info.executionMode` | what actually ran (after GPU/PCT fallback) |
| `info.numWorkers` | parallel pool size (1 if serial) |
| `info.chunkMode` / `info.chunkSize` / `info.numChunks` | chunking actually used |
| `diag.runtimeSeconds` | wall-clock time of the fusion reconstruction |

Set `cfg.performance.verbose = true` to print the line
`[backProjectionCore] requested=… resolved=… actual=… workers=… chunkMode=… chunkSize=… chunks=…`.

---

## 7. Benchmark methodology  (`experimentRuntimeBenchmark.m`)

Small, **direct-feasible** problem (single point target, 3-D window,
`config('quick')` sizes). For each method:

1. one warm-up call (JIT / pool spin-up excluded from timing),
2. one timed call,
3. peak-normalise the image, compare to the **direct** reference
   (`‖imgN − refN‖ / ‖refN‖`),
4. peak-location error vs the true target,
5. `estimateChunkSize` memory estimate, actual execution mode, chunk
   size.

Recorded per method: `runtimeSeconds`, `speedupVsDirect`,
`speedupVsFFTSerial`, `peakLocationError`,
`relativeImageErrorVsReference`, `memoryEstimateMB`,
`actualExecutionMode`, `chunkSizeUsed`.

**Unavailable modes are written as `unavailable` in the CSV and omitted
from the plots — never given invented numbers**, and no speed-up is
reported for them.

Outputs (under `cfg.resultsDir`):
`runtime_benchmark.mat`, `runtime_benchmark.csv`,
`runtime_comparison.png`, `speedup_comparison.png`,
`reconstruction_error_comparison.png`.

---

## 8. Hardware / toolbox dependencies

| Feature | Requires | Fallback |
|---|---|---|
| `direct`, `fft_serial`, auto chunking | base MATLAB only | — |
| `fft_parallel` | Parallel Computing Toolbox (installed) | serial + 1 warning |
| `fft_gpu` | Parallel Computing Toolbox + supported GPU | CPU + 1 warning |
| Windows memory probe in auto chunking | `memory` builtin (Windows) | configured `targetMemoryMB` |

CPU-only, base-MATLAB operation is fully functional and is the default
and tested path.

---

## 9. Limitations

* The auto chunk-size model accounts for the dominant arrays only; a
  pathological `Nv` can still exceed the budget — lower
  `targetMemoryMB` or set an explicit `cfg.performance.chunkSize`.
* GPU / parallel paths could not be executed on the development machine
  (no PCT / GPU); their correctness relies on the unchanged pre-existing
  code paths plus the skipped-not-faked validation checks.
* The benchmark uses one small scene; absolute speed-ups scale with
  problem size and hardware.
