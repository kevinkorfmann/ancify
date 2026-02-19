# Quickstart

Get from zero to ancestral alleles in five minutes.

## 1. Install

```bash
pip install .
# or with uv:
uv pip install .
```

## 2. Generate a config template

```bash
ancify init -o config.yaml
```

## 3. Edit the config

Open `config.yaml` and fill in:

- **`chromosome_lengths`** — path to a tab-separated file with chromosome names and lengths.
- **`outgroups.inner`** — one or more closely related species with their net AXT alignment files.
- **`outgroups.outer`** — one or more distantly related species with their net AXT alignment files.
- **`output_dir`** — where the ancestral FASTA files should go.

Minimal example for human (hg38):

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt

outgroups:
  inner:
    - name: bonobo
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla
      alignment: hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque
      alignment: hg38.rheMac10.net.axt.gz

output_dir: ./ancestral_calls
num_cpus: 24
```

## 4. Run the pipeline

```bash
ancify run -c config.yaml
```

This runs all three phases in sequence:

1. **Project** outgroup alignments onto focal-species coordinates.
2. **Call** the ancestral allele at every position.
3. **Evaluate** (optional, if the `evaluation` block is present in the config).

## 5. Use the output

The output is one FASTA file per chromosome in `output_dir/`.  Each position
contains a single character encoding the ancestral call and its confidence:

| Character | Meaning |
|-----------|---------|
| `ACGT` | High confidence — inner and outer outgroups agree |
| `acgt` | Low confidence — only one tier has data |
| `n` | Unresolved — inner and outer disagree |
| `N` | Missing — no data from either tier |

Look up an ancestral allele in Python:

```python
from ancify.utils import read_fasta

_, seq = read_fasta("ancestral_calls/chr1.fa")
allele = seq[999999]  # 0-based index for position 1,000,000
print(f"Ancestral: {allele}, High confidence: {allele in 'ACGT'}")
```

## Next steps

- {doc}`configuration` — full reference for all config fields.
- {doc}`algorithm` — how the two-tier outgroup voting works.
- {doc}`species_guide` — how to adapt the pipeline for mouse, Drosophila, etc.
- {doc}`api` — Python API reference for programmatic use.
