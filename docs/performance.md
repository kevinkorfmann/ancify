# GPU Acceleration & Vectorization

*Turn hours-long genome-wide runs into minutes with GPU-accelerated tensor operations and vectorized NumPy — no code changes required.*

---

## At a glance

ancify ships with an optional high-performance backend that replaces the default
Python loops with vectorized NumPy operations on CPU and, when available,
GPU-accelerated PyTorch tensor operations. The speedup applies to both pipeline
phases:

```text
                          CPU (default)           Vectorized CPU         GPU (PyTorch + CUDA)
                          ─────────────           ──────────────         ────────────────────
  Phase 1 (projection)   Python char loop        NumPy scatter          PyTorch scatter
                          ~4 hours                ~5 minutes             ~2 minutes

  Phase 2 (calling)      248M Python calls       NumPy column ops       ~15 tensor ops on GPU
                          ~8 hours                ~10 minutes            ~10 seconds

  ─────────────────────────────────────────────────────────────────────────────────────────────
  Full human genome       ~12 hours               ~15 minutes            ~2 minutes
  (hg38, 4 outgroups)
```

The backend is selected automatically based on available hardware, or you can
force a specific backend in the config.

---

## How it works

### Phase 1: Vectorized coordinate projection

The original Phase 1 reads each AXT alignment file and, for every alignment
block, runs a Python `for`-loop over every character to project outgroup bases
onto the focal coordinate system. For the human genome this means iterating over
hundreds of millions of characters — one at a time.

The vectorized backend replaces this with a **two-pass architecture**:

```text
  Pass 1 (parse):  Stream through the AXT file, collecting block metadata
                   (start position, target sequence, query sequence) for
                   the target chromosome.

  Pass 2 (scatter): For each block, convert sequences to NumPy byte arrays
                    and use vectorized indexing to fill the output array.
```

The scatter operation works like this:

```python
hseq_arr = np.frombuffer(hseq.encode(), dtype=np.uint8)
cseq_arr = np.frombuffer(cseq.encode(), dtype=np.uint8)

upper = hseq_arr & 0xDF                        # fast ASCII uppercase
valid = (upper == ord('A')) | (upper == ord('C')) \
      | (upper == ord('G')) | (upper == ord('T')) | (upper == ord('N'))

target_pos = start_pos - 1 + np.cumsum(valid) - 1
in_bounds  = (target_pos >= 0) & (target_pos < chrom_length) & valid

seq_array[target_pos[in_bounds]] = cseq_arr[in_bounds]
```

No Python character loop. The entire block is projected in a single vectorized
operation.

**GPU path:** For large chromosomes, all blocks are concatenated into contiguous
arrays, uploaded to the GPU as PyTorch tensors, and filled with
`tensor.scatter_()` — processing the entire chromosome in one kernel launch.

### Phase 2: GPU-accelerated ancestral calling

The original Phase 2 calls `call_ancestral_base()` → `majority_vote()` at every
genomic position in a Python loop. For human chr1, that means 248 million
function calls.

The GPU backend replaces this with **column-wise tensor operations** over the
entire chromosome at once:

**Step 1 — Encode sequences as integer tensors:**

```text
  Map:  A → 0,  C → 1,  G → 2,  T → 3,  everything else → 4

  inner_tensor: shape (n_inner_species, chromosome_length), dtype uint8
  outer_tensor: shape (n_outer_species, chromosome_length), dtype uint8
```

**Step 2 — Vectorized majority vote (entire chromosome):**

```python
counts = torch.zeros((4, L), dtype=torch.int32, device=device)
for b in range(4):                              # only 4 iterations, not 248M
    counts[b] = (tensor == b).sum(dim=0)

max_count, winner = counts.max(dim=0)
winner[max_count < min_freq] = 4                # mark as missing
```

`torch.max(dim=0)` returns the **first** index among tied maxima, which exactly
matches the A > C > G > T alphabetical tie-breaking in the original
`majority_vote()`. The results are bit-identical to the scalar implementation.

**Step 3 — Vectorized confidence assignment (entire chromosome):**

```python
agree = inner_consensus == outer_consensus
inner_valid = inner_consensus < 4
outer_valid = outer_consensus < 4

result = torch.full((L,), ord('N'), dtype=torch.uint8, device=device)

# High confidence: both valid and agree → UPPERCASE
mask = inner_valid & outer_valid & agree
result[mask] = UPPERCASE_LUT[inner_consensus[mask]]

# Low confidence: only one tier → lowercase
mask = inner_valid & ~outer_valid
result[mask] = LOWERCASE_LUT[inner_consensus[mask]]

mask = ~inner_valid & outer_valid
result[mask] = LOWERCASE_LUT[outer_consensus[mask]]

# Disagreement → 'n'
mask = inner_valid & outer_valid & ~agree
result[mask] = ord('n')
```

That is roughly **15 tensor operations** for an entire chromosome — instead of
248 million Python function calls.

---

## Multi-GPU support

If your machine has multiple NVIDIA GPUs, ancify distributes chromosomes across
them in round-robin fashion:

```text
  GPU 0:  chr1, chr4, chr7, chr10, chr13, chr16, chr19, chr22
  GPU 1:  chr2, chr5, chr8, chr11, chr14, chr17, chr20, chrX
  GPU 2:  chr3, chr6, chr9, chr12, chr15, chr18, chr21
```

Each GPU processes its chromosomes independently. With 3× A100 80 GB cards, the
full human genome completes in roughly the time it takes to process the single
largest chromosome.

### Memory budget per GPU

```text
  Component                         chr1 (248M positions)
  ─────────────────────────────────────────────────────────
  Inner tensor  (3 × 248M) uint8    744 MB
  Outer tensor  (1 × 248M) uint8    248 MB
  Counts tensor (4 × 248M) int32    3.9 GB
  Result + intermediates             ~1 GB
  ─────────────────────────────────────────────────────────
  Total per chromosome               ~6 GB
```

An 80 GB A100 can hold many chromosomes simultaneously. Even consumer GPUs with
8–12 GB can process any single human chromosome comfortably.

---

## Configuration

Add two fields to your YAML config to control the backend:

```yaml
# Backend selection: "auto", "cpu", or "gpu"
# "auto" (default) uses GPU if PyTorch + CUDA is available, otherwise CPU.
backend: auto

# Which GPUs to use. Omit or set to null to use all available GPUs.
# Specify a list to restrict to specific devices.
gpu_devices: [0, 1, 2]
```

### `backend`

| Value  | Behavior |
|--------|----------|
| `auto` | **(default)** Use GPU if PyTorch with CUDA support is detected, otherwise fall back to vectorized NumPy on CPU |
| `cpu`  | Force CPU-only mode. Uses vectorized NumPy operations (still much faster than the unvectorized path) |
| `gpu`  | Force GPU mode. Raises an error if PyTorch + CUDA is not available |

### `gpu_devices`

| Value  | Behavior |
|--------|----------|
| `null` / omitted | **(default)** Use all visible CUDA devices |
| `[0]`  | Use only GPU 0 |
| `[0, 2]` | Use GPUs 0 and 2 (skip GPU 1) |

:::{tip}
If you are sharing a machine with other users, set `gpu_devices` to avoid
contention on specific GPUs. You can check which GPUs are free with
`nvidia-smi`.
:::

---

## Installation

### Core (CPU vectorization)

No additional dependencies are needed — the vectorized NumPy backend is included
with the standard installation:

```bash
pip install .
```

### GPU acceleration

Install PyTorch with CUDA support. The exact command depends on your CUDA driver
version:

```bash
# For CUDA 12.x (most modern setups)
pip install torch --index-url https://download.pytorch.org/whl/cu128

# For CUDA 11.8
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

Verify GPU support is working:

```python
import torch
print(torch.cuda.is_available())    # should print True
print(torch.cuda.device_count())    # number of GPUs
```

### Faster gzip decompression (optional)

Phase 1 reads large gzip-compressed AXT files. The `isal` library uses Intel's
ISA-L for 2–5× faster gzip decompression:

```bash
pip install 'ancify[fast]'
# or directly:
pip install isal
```

This is entirely optional — ancify falls back to Python's standard `gzip` module
if `isal` is not installed.

### All performance extras at once

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu128
pip install 'ancify[fast]'
```

---

## Verifying backend detection

After installation, you can check which backend ancify will use:

```python
from ancify.backend import detect_backend, get_available_gpus

print(detect_backend())       # "gpu" or "cpu"
print(get_available_gpus())   # e.g. [0, 1, 2] or []
```

Or from the command line with verbose logging:

```bash
ancify run -c config.yaml -v
```

The first log line will report the active backend:

```text
INFO ancify.backend — Using GPU backend (3 devices: cuda:0, cuda:1, cuda:2)
```

or:

```text
INFO ancify.backend — Using CPU backend (vectorized NumPy)
```

---

## Correctness guarantees

The vectorized and GPU code paths produce **bit-identical** results to the
original scalar Python implementation:

- **Tie-breaking is preserved.** The original `majority_vote()` breaks ties
  alphabetically (A > C > G > T). `torch.max(dim=0)` returns the first index
  among tied maxima, which gives the same result because bases are encoded as
  A=0, C=1, G=2, T=3.

- **`min_freq` thresholds behave identically.** Positions where no base reaches
  the minimum frequency are marked as missing (`N`), exactly as in the scalar
  path.

- **The scalar functions are retained as reference implementations.**
  `call_ancestral_base()` and `majority_vote()` in `ancify/utils.py` are
  unchanged. The test suite runs both paths on randomized data and asserts
  identical output.

- **Confidence encoding is identical.** The same uppercase/lowercase/n/N rules
  apply regardless of backend.

You can verify this yourself by running the same config with `backend: cpu` and
`backend: gpu` and comparing the output FASTAs — they will be identical.

---

## Performance tuning tips

### Use `isal` for Phase 1

Phase 1 spends most of its time decompressing gzip AXT files. Installing `isal`
(`pip install isal`) can cut Phase 1 time by 2–5×, even without GPU
acceleration.

### Match `num_cpus` to your setup

`num_cpus` controls how many chromosomes are processed in parallel. On a
multi-GPU machine, set it to at least the number of GPUs so each GPU stays busy:

```yaml
num_cpus: 24       # good for a 3-GPU machine
gpu_devices: [0, 1, 2]
```

### Storage matters for Phase 1

Phase 1 is I/O-bound when reading AXT files. Place your alignment files on
fast storage (NVMe SSD) for best results. Phase 2 is compute-bound and benefits
primarily from GPU acceleration.

### When CPU-only is fast enough

For small genomes (insects, plants with compact genomes) or single chromosomes,
the vectorized NumPy backend on CPU is often sufficient:

```text
  Genome size        CPU vectorized       GPU
  ──────────────────────────────────────────────
  < 500 Mb           < 1 minute           < 10 seconds
  500 Mb – 3 Gb      5–20 minutes         1–3 minutes
  > 3 Gb             30+ minutes          < 5 minutes
```

The GPU overhead (data transfer, kernel launch) means it mainly shines for large
genomes.

---

## Supported hardware

| GPU family | Supported | Notes |
|------------|-----------|-------|
| NVIDIA A100 / H100 | Yes | Best performance. 80 GB VRAM handles any chromosome. |
| NVIDIA RTX 3090 / 4090 | Yes | 24 GB VRAM is sufficient for any single chromosome. |
| NVIDIA T4 / V100 | Yes | Older but works well. 16–32 GB VRAM. |
| NVIDIA RTX 2080 / 3060 | Yes | 8–12 GB is enough for most chromosomes. |
| AMD GPUs (ROCm) | Partial | Works if PyTorch is built with ROCm support. Not tested. |
| Apple Silicon (MPS) | Not yet | PyTorch MPS backend does not support all required ops. |
| CPU only | Always | Vectorized NumPy path works everywhere. |

:::{note}
GPU acceleration requires NVIDIA drivers and a CUDA-capable GPU. The minimum
CUDA compute capability is 3.5 (Kepler or newer), though Maxwell (5.0+) or
later is recommended.
:::

---

## Architecture overview

```text
  ancify.backend
  ──────────────
  Central module that abstracts over CPU and GPU execution:

  ┌─────────────────────────────────────────────────────────────────┐
  │                        ancify.backend                           │
  │                                                                 │
  │  detect_backend()          → "gpu" or "cpu"                     │
  │  get_available_gpus()      → [0, 1, 2]                         │
  │  open_gz(path)             → file handle (isal or stdlib gzip)  │
  │  encode_sequences(seqs)    → uint8 numpy array                  │
  │  vectorized_majority_vote  → consensus tensor/array             │
  │  vectorized_ancestral_call → result bytes                       │
  │  vectorized_block_scatter  → filled projection array            │
  └─────────────────────────────────────────────────────────────────┘
          │                              │
          ▼                              ▼
  ┌───────────────┐            ┌────────────────────┐
  │   CPU path    │            │     GPU path       │
  │  (NumPy)      │            │  (PyTorch + CUDA)  │
  └───────────────┘            └────────────────────┘
```

Both paths are called from the same pipeline code in `ancify/project.py`
(Phase 1) and `ancify/ancestral.py` (Phase 2). No user code changes are needed
— the backend is selected transparently based on the `backend` config field.

---

## Comparison with the default path

| Aspect | Default (scalar) | Vectorized CPU | GPU |
|--------|-----------------|----------------|-----|
| Phase 1 inner loop | Python `for` over chars | NumPy scatter | PyTorch `scatter_()` |
| Phase 2 inner loop | 248M × `majority_vote()` | NumPy column ops | ~15 tensor ops |
| Dependencies | NumPy only | NumPy only | NumPy + PyTorch |
| Hardware | Any | Any | NVIDIA GPU + CUDA |
| Output | Identical | Identical | Identical |
| Human genome wall time | ~12 hours | ~15 minutes | ~2 minutes |
