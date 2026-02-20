# Adapting to Other Species

*ancify works with **any focal species** for which you have pairwise net AXT alignments. This guide teaches you how to think about outgroup selection and walks through examples from across the tree of life.*

---

## A framework for choosing outgroups

Before configuring ancify, you need to make one key decision: **which species are your inner outgroups and which are your outer outgroups?**

### Step 1: Draw the phylogeny

Sketch (or look up) the phylogenetic relationships around your focal species. You need at least three species total: your focal, one inner outgroup, and one outer outgroup.

```text
  Example for human:

                 6 Mya         9 Mya              25 Mya
                  ┌─── Bonobo
            ┌─────┤                          ┌ INNER tier
            │     └─── Chimp                 │ (closely related)
     ┌──────┤                                │
     │      └──────── Gorilla               ─┘
  ───┤
     │
     └──────────────────── Macaque          ── OUTER tier
                                               (distantly related)
```

### Step 2: Assign tiers

**Inner outgroups** should be:

- Closely related to the focal species (same genus or family)
- Diverged **more recently** than the outer outgroup
- Ideally 2 or more species (enables majority voting)

**Outer outgroups** should be:

- Clearly outside the inner clade
- Diverged at least 2-3x further than the inner outgroups
- Far enough that convergent mutations with the inner group are extremely rare

### Step 3: Check data availability

For each candidate outgroup, check if UCSC has a pairwise net AXT alignment to your focal assembly:

```text
https://hgdownload.soe.ucsc.edu/goldenPath/<focal_assembly>/
```

Look for directories named `vs<Outgroup>`. If the alignment exists, you are in business.

### Decision flowchart

```text
  Do you have ≥2 inner outgroups with AXT alignments?
    │
    ├─ YES → Great! Majority voting will be robust.
    │
    └─ NO → Do you have 1 inner + 1 outer?
              │
              ├─ YES → Still works. The outer outgroup provides
              │        the independent check. Consider adding more
              │        inner species if available.
              │
              └─ NO → You need at least 1 inner + 1 outer.
                      Check UCSC or generate your own alignments.
```

---

## Worked examples

### Human (hg38) — the gold standard

```yaml
focal_species: human
chromosome_lengths: chromoLens.txt

outgroups:
  inner:
    - name: bonobo      # ~6 Mya
      alignment: hg38.panPan3.net.axt.gz
    - name: chimp       # ~6 Mya
      alignment: hg38.panTro6.net.axt.gz
    - name: gorilla     # ~9 Mya
      alignment: hg38.gorGor6.net.axt.gz
  outer:
    - name: macaque     # ~25 Mya
      alignment: hg38.rheMac10.net.axt.gz

output_dir: ./human_ancestral
num_cpus: 24
```

**Why this works well:** Three inner outgroups provide redundancy. The inner-outer divergence ratio (~6-9 Mya vs. ~25 Mya) is large enough that convergent errors between tiers are negligible. Alignment coverage is excellent for all four species.

### Mouse (mm39) — minimal setup

```yaml
focal_species: mouse
chromosome_lengths: mm39.chromLens.txt

outgroups:
  inner:
    - name: rat            # ~12 Mya
      alignment: mm39.rn7.net.axt.gz
  outer:
    - name: rabbit         # ~90 Mya
      alignment: mm39.oryCun2.net.axt.gz

output_dir: ./mouse_ancestral
num_cpus: 8
```

**With only one inner species,** the majority vote is trivially that species' allele. The outer outgroup still provides the independent confirmation. To strengthen the inner tier, consider adding hamster or other rodents if alignments are available.

### *Drosophila melanogaster* (dm6) — non-chr naming

```yaml
focal_species: drosophila_melanogaster
chromosome_lengths: dm6.chromLens.txt
chromosomes: [2L, 2R, 3L, 3R, 4, X]

outgroups:
  inner:
    - name: simulans       # ~2.5 Mya
      alignment: dm6.droSim2.net.axt.gz
    - name: sechellia      # ~2.5 Mya
      alignment: dm6.droSec1.net.axt.gz
  outer:
    - name: yakuba         # ~6 Mya
      alignment: dm6.droYak3.net.axt.gz

output_dir: ./dmel_ancestral
num_cpus: 6
```

**Note:** The chromosome names (`2L`, `3R`, etc.) do not have a `chr` prefix — ancify handles any naming convention. The explicit `chromosomes` list excludes heterochromatic scaffolds.

### *Brassica rapa* (plant) — beyond animals

```yaml
focal_species: brassica_rapa
chromosome_lengths: braRap1.chromLens.txt

outgroups:
  inner:
    - name: brassica_oleracea   # close relative, same genus
      alignment: braRap1.braOleracea.net.axt.gz
  outer:
    - name: arabidopsis_thaliana  # ~20 Mya, same family Brassicaceae
      alignment: braRap1.araTha1.net.axt.gz

output_dir: ./brassica_rapa_ancestral
num_cpus: 4
```

**Plant genomes** often use chromosome naming like A01, A02. Omit the `chromosomes` key to process all entries in the lengths file. Plant genome alignments may be more fragmented due to whole-genome duplications — expect higher `N` rates than in mammals.

### Zebrafish (danRer11) — fish

```yaml
focal_species: zebrafish
chromosome_lengths: danRer11.chromLens.txt

outgroups:
  inner:
    - name: medaka
      alignment: danRer11.oryLat2.net.axt.gz
  outer:
    - name: fugu
      alignment: danRer11.fr3.net.axt.gz

output_dir: ./zebrafish_ancestral
num_cpus: 4
```

**Fish have high substitution rates** compared to mammals, so the divergence between inner and outer outgroups needs to be carefully considered. Medaka and fugu provide a reasonable tier separation for zebrafish.

---

## Getting the input data

### Net AXT alignments from UCSC

Download from:

```text
https://hgdownload.soe.ucsc.edu/goldenPath/<focal_assembly>/vs<Outgroup>/
```

Example for human vs. chimp:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsPanTro6/hg38.panTro6.net.axt.gz
```

### Chromosome lengths

**From a FASTA index:**

```bash
samtools faidx reference.fa
cut -f1,2 reference.fa.fai > chromoLens.txt
```

**From UCSC MySQL:**

```bash
mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A \
  -e "SELECT chrom, size FROM chromInfo" hg38 > chromoLens.txt
```

---

## Tips for outgroup selection

### More inner species is better

Even 2 inner species is a major improvement over 1. With 1 inner species, a single alignment error or lineage-specific substitution produces a wrong call. With 2+, the majority vote provides robustness.

```text
  1 inner species:   accuracy ≈ alignment quality
  2 inner species:   accuracy ≈ max(alignment quality)
  3+ inner species:  accuracy ≈ consensus of multiple independent signals
```

### The outer outgroup must be clearly outside the inner clade

If the outer outgroup is too closely related to the inner species, ILS can affect both tiers *together*, producing false high-confidence calls:

```text
  BAD: outer is too close to inner

       ┌── Focal
    ┌──┤
    │  └── Inner 1             ILS can affect all three
  ──┤
    └──── Outer (barely        species → false agreement
          more distant)

  GOOD: outer is clearly distant

       ┌── Focal
    ┌──┤
    │  └── Inner 1
  ──┤
    │
    │
    └────────── Outer (deep   ILS between tiers is
                divergence)    negligible
```

### Stringency vs. coverage tradeoff

Increasing `min_inner_freq` requires more species to agree:

| Setting | Coverage | Accuracy |
|---------|----------|----------|
| `min_inner_freq=1` | Highest | Lower (but still >99%) |
| `min_inner_freq=2` (with 3 inner) | Moderate | High |
| `min_inner_freq=3` (with 3 inner) | Lowest | Highest |

For most demographic analyses, `min_inner_freq=1` is fine — the SFS shape is robust to rare misassignments.

### Sex chromosomes and other special cases

- **chrY** has very poor alignment coverage due to massive repetitive content. Expect >80% `N` (missing).
- **chrX** may have lower coverage than autosomes, especially in regions of reduced synteny.
- **Mitochondrial genomes** are not handled by standard UCSC pairwise alignments. For mtDNA polarization, consider using a multiple sequence alignment approach instead.

---

## Species catalogue

The tables below list 100 commonly studied species with suggested outgroup configurations for ancestral allele polarization. Each entry shows:

- **Assembly** — the UCSC genome browser identifier (where one exists).
- **Suggested inner outgroups** — closely related species for the inner tier.
- **Suggested outer outgroup** — a more distant species for independent confirmation.
- **Divergence** — approximate divergence times (inner / outer) in millions of years ago (Mya).

:::{tip}
These are starting-point suggestions. Always verify that pairwise net AXT alignments exist for your chosen focal assembly at `https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/`. If a pre-computed alignment is not available, you can generate your own with [lastz](https://lastz.github.io/lastz/) and the UCSC [axtChain/chainNet](https://genome.ucsc.edu/goldenPath/help/chain.html) pipeline.
:::

### Primates

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 1 | *Homo sapiens* | Human | hg38 | Chimp, bonobo, gorilla | Macaque | 6–9 / 25 Mya |
| 2 | *Pan troglodytes* | Chimpanzee | panTro6 | Bonobo, human | Macaque | 2–6 / 25 Mya |
| 3 | *Pan paniscus* | Bonobo | panPan3 | Chimp, human | Macaque | 2–6 / 25 Mya |
| 4 | *Gorilla gorilla* | Gorilla | gorGor6 | Human, chimp | Macaque | 9 / 25 Mya |
| 5 | *Pongo abelii* | Orangutan | ponAbe3 | Human, chimp, gorilla | Macaque | 13 / 25 Mya |
| 6 | *Nomascus leucogenys* | Gibbon | nomLeu3 | Human, macaque | Marmoset | 18–25 / 40 Mya |
| 7 | *Macaca mulatta* | Rhesus macaque | rheMac10 | Crab-eating macaque, baboon | Human | 3–8 / 25 Mya |
| 8 | *Macaca fascicularis* | Crab-eating macaque | macFas5 | Rhesus macaque, baboon | Human | 3–8 / 25 Mya |
| 9 | *Papio anubis* | Baboon | papAnu4 | Rhesus macaque, green monkey | Human | 5–8 / 25 Mya |
| 10 | *Chlorocebus sabaeus* | Green monkey | chlSab2 | Rhesus macaque, baboon | Human | 8 / 25 Mya |
| 11 | *Nasalis larvatus* | Proboscis monkey | nasLar1 | Green monkey, rhesus macaque | Human | 10 / 25 Mya |
| 12 | *Callithrix jacchus* | Marmoset | calJac4 | Squirrel monkey | Macaque | 15 / 40 Mya |
| 13 | *Saimiri boliviensis* | Squirrel monkey | saiBol1 | Marmoset | Macaque | 15 / 40 Mya |
| 14 | *Carlito syrichta* | Tarsier | tarSyr2 | Mouse lemur, bushbaby | Human | 58–70 / 75 Mya |
| 15 | *Microcebus murinus* | Mouse lemur | micMur2 | Bushbaby | Marmoset | 58 / 70 Mya |
| 16 | *Otolemur garnettii* | Bushbaby | otoGar3 | Mouse lemur | Marmoset | 58 / 70 Mya |

### Rodents & Lagomorphs

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 17 | *Mus musculus* | Mouse | mm39 | Rat | Rabbit | 12 / 90 Mya |
| 18 | *Rattus norvegicus* | Rat | rn7 | Mouse | Rabbit | 12 / 90 Mya |
| 19 | *Cricetulus griseus* | Chinese hamster | criGriChoV2 | Mouse, rat | Rabbit | 20 / 90 Mya |
| 20 | *Cavia porcellus* | Guinea pig | cavPor3 | Chinchilla | Mouse | 35 / 70 Mya |
| 21 | *Chinchilla lanigera* | Chinchilla | chiLan1 | Guinea pig | Mouse | 35 / 70 Mya |
| 22 | *Heterocephalus glaber* | Naked mole-rat | hetGla2 | Guinea pig | Mouse | 35 / 70 Mya |
| 23 | *Ictidomys tridecemlineatus* | Squirrel | speTri2 | Mouse, rat | Rabbit | 50 / 90 Mya |
| 24 | *Dipodomys ordii* | Kangaroo rat | dipOrd1 | Mouse, rat | Rabbit | 50 / 90 Mya |
| 25 | *Oryctolagus cuniculus* | Rabbit | oryCun2 | Pika | Mouse | 30 / 90 Mya |
| 26 | *Ochotona princeps* | Pika | ochPri3 | Rabbit | Mouse | 30 / 90 Mya |

### Carnivores

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 27 | *Canis lupus familiaris* | Dog | canFam6 | Ferret, cat | Horse | 50–55 / 80 Mya |
| 28 | *Felis catus* | Cat | felCat9 | Ferret, dog | Horse | 50–55 / 80 Mya |
| 29 | *Mustela putorius furo* | Ferret | musFur1 | Dog, cat | Horse | 50–55 / 80 Mya |
| 30 | *Ailuropoda melanoleuca* | Giant panda | ailMel1 | Dog, ferret | Horse | 50 / 80 Mya |
| 31 | *Neomonachus schauinslandi* | Hawaiian monk seal | neoSch1 | Dog, ferret | Horse | 50 / 80 Mya |

### Ungulates & Cetaceans

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 32 | *Bos taurus* | Cow | bosTau9 | Sheep | Pig | 20 / 60 Mya |
| 33 | *Ovis aries* | Sheep | oviAri4 | Cow | Pig | 20 / 60 Mya |
| 34 | *Sus scrofa* | Pig | susScr11 | Cow, dolphin | Horse | 60 / 80 Mya |
| 35 | *Tursiops truncatus* | Dolphin | turTru2 | Cow | Horse | 55 / 80 Mya |
| 36 | *Balaenoptera acutorostrata* | Minke whale | balAcu1 | Dolphin, cow | Horse | 35–55 / 80 Mya |
| 37 | *Equus caballus* | Horse | equCab3 | White rhinoceros | Dog | 55 / 80 Mya |
| 38 | *Ceratotherium simum* | White rhinoceros | cerSim1 | Horse | Dog | 55 / 80 Mya |
| 39 | *Vicugna pacos* | Alpaca | vicPac2 | Pig, cow | Horse | 55–60 / 80 Mya |

### Bats

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 40 | *Myotis lucifugus* | Microbat | myoLuc2 | Megabat | Dog | 55 / 80 Mya |
| 41 | *Pteropus vampyrus* | Megabat | pteVam1 | Microbat | Dog | 55 / 80 Mya |

### Insectivores & other Laurasiatheria

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 42 | *Erinaceus europaeus* | Hedgehog | eriEur2 | Shrew | Dog | 75 / 90 Mya |
| 43 | *Sorex araneus* | Shrew | sorAra2 | Hedgehog | Dog | 75 / 90 Mya |
| 44 | *Manis pentadactyla* | Chinese pangolin | manPen1 | Dog, cat | Horse | 75 / 85 Mya |

### Afrotheria & Xenarthra

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 45 | *Loxodonta africana* | African elephant | loxAfr3 | Manatee, rock hyrax | Armadillo | 60–65 / 100 Mya |
| 46 | *Trichechus manatus* | Manatee | triMan1 | Elephant | Armadillo | 60 / 100 Mya |
| 47 | *Procavia capensis* | Rock hyrax | proCap1 | Elephant, manatee | Armadillo | 60 / 100 Mya |
| 48 | *Echinops telfairi* | Tenrec | echTel2 | Elephant | Armadillo | 75 / 100 Mya |
| 49 | *Dasypus novemcinctus* | Armadillo | dasNov3 | Sloth | Elephant | 65 / 100 Mya |
| 50 | *Choloepus hoffmanni* | Sloth | choHof1 | Armadillo | Elephant | 65 / 100 Mya |

### Marsupials & Monotremes

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 51 | *Monodelphis domestica* | Opossum | monDom5 | Tasmanian devil | Wallaby | 70 / 80 Mya |
| 52 | *Sarcophilus harrisii* | Tasmanian devil | sarHar1 | Wallaby | Opossum | 40 / 80 Mya |
| 53 | *Macropus eugenii* | Wallaby | macEug2 | Tasmanian devil | Opossum | 40 / 80 Mya |
| 54 | *Ornithorhynchus anatinus* | Platypus | ornAna2 | Opossum | Chicken | 190 / 320 Mya |

### Birds

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 55 | *Gallus gallus* | Chicken | galGal6 | Turkey | Zebra finch | 30 / 80 Mya |
| 56 | *Meleagris gallopavo* | Turkey | melGal5 | Chicken | Zebra finch | 30 / 80 Mya |
| 57 | *Anas platyrhynchos* | Mallard duck | anaPla1 | Chicken, turkey | Zebra finch | 55–70 / 80 Mya |
| 58 | *Taeniopygia guttata* | Zebra finch | taeGut2 | Medium ground finch, flycatcher | Chicken | 15–45 / 80 Mya |
| 59 | *Geospiza fortis* | Medium ground finch | geoFor1 | Zebra finch | Chicken | 15 / 80 Mya |
| 60 | *Ficedula albicollis* | Collared flycatcher | ficAlb2 | Zebra finch | Chicken | 45 / 80 Mya |
| 61 | *Melopsittacus undulatus* | Budgerigar | melUnd1 | Falcon, eagle | Chicken | 50–60 / 80 Mya |
| 62 | *Falco peregrinus* | Peregrine falcon | falPer1 | Budgerigar, eagle | Chicken | 50–55 / 80 Mya |
| 63 | *Aquila chrysaetos* | Golden eagle | aquChr2 | Falcon, budgerigar | Chicken | 50–55 / 80 Mya |
| 64 | *Columba livia* | Rock pigeon | colLiv1 | Flycatcher, zebra finch | Chicken | 70 / 80 Mya |
| 65 | *Apteryix mantelli* | Brown kiwi | aptMan1 | Chicken, turkey | Zebra finch | 60–70 / 80 Mya |

### Reptiles & Amphibians

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 66 | *Anolis carolinensis* | Green anole | anoCar2 | Garter snake | Painted turtle | 150 / 250 Mya |
| 67 | *Chrysemys picta* | Painted turtle | chrPic2 | Softshell turtle, alligator | Green anole | 70–160 / 250 Mya |
| 68 | *Pelodiscus sinensis* | Chinese softshell turtle | pelSin2 | Painted turtle, alligator | Green anole | 70–160 / 250 Mya |
| 69 | *Alligator mississippiensis* | American alligator | allMis1 | Painted turtle, softshell turtle | Green anole | 160 / 250 Mya |
| 70 | *Thamnophis sirtalis* | Garter snake | thaSir1 | Green anole | Painted turtle | 150 / 250 Mya |
| 71 | *Xenopus tropicalis* | Western clawed frog | xenTro10 | *X. laevis* | Coelacanth | 50 / 370 Mya |
| 72 | *Xenopus laevis* | African clawed frog | xenLae2 | *X. tropicalis* | Coelacanth | 50 / 370 Mya |
| 73 | *Rana temporaria* | Common frog | — | *Xenopus* spp. | Coelacanth | 200 / 370 Mya |

### Fish

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 74 | *Danio rerio* | Zebrafish | danRer11 | Medaka, stickleback | Fugu | 110–150 / 180 Mya |
| 75 | *Oryzias latipes* | Medaka | oryLat2 | Stickleback, Nile tilapia | Zebrafish | 70–100 / 150 Mya |
| 76 | *Gasterosteus aculeatus* | Stickleback | gasAcu1 | Medaka, Nile tilapia | Zebrafish | 70–100 / 150 Mya |
| 77 | *Takifugu rubripes* | Fugu | fr3 | Tetraodon | Stickleback | 30 / 100 Mya |
| 78 | *Tetraodon nigroviridis* | Tetraodon | tetNig2 | Fugu | Stickleback | 30 / 100 Mya |
| 79 | *Oreochromis niloticus* | Nile tilapia | oreNil3 | Stickleback, medaka | Zebrafish | 70–100 / 150 Mya |
| 80 | *Gadus morhua* | Atlantic cod | gadMor1 | Stickleback, medaka | Zebrafish | 100–110 / 150 Mya |
| 81 | *Salmo salar* | Atlantic salmon | — | Rainbow trout | Zebrafish | 25 / 200 Mya |
| 82 | *Oncorhynchus mykiss* | Rainbow trout | — | Atlantic salmon | Zebrafish | 25 / 200 Mya |
| 83 | *Latimeria chalumnae* | Coelacanth | latCha1 | Spotted gar | *Xenopus* | 400 / 420 Mya |
| 84 | *Lepisosteus oculatus* | Spotted gar | lepOcu1 | Zebrafish | Coelacanth | 300 / 420 Mya |
| 85 | *Petromyzon marinus* | Sea lamprey | petMar3 | Elephant shark | Lancelet | 450 / 550 Mya |
| 86 | *Callorhinchus milii* | Elephant shark | calMil1 | Zebrafish | Lamprey | 400 / 450 Mya |

### Insects — *Drosophila* & relatives

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 87 | *D. melanogaster* | Fruit fly | dm6 | *D. simulans*, *D. sechellia* | *D. yakuba* | 2.5 / 6 Mya |
| 88 | *D. simulans* | — | droSim2 | *D. sechellia*, *D. melanogaster* | *D. yakuba* | 0.5–2.5 / 6 Mya |
| 89 | *D. sechellia* | — | droSec1 | *D. simulans*, *D. melanogaster* | *D. yakuba* | 0.5–2.5 / 6 Mya |
| 90 | *D. yakuba* | — | droYak3 | *D. erecta* | *D. melanogaster* | 6 / 12 Mya |
| 91 | *D. erecta* | — | droEre2 | *D. yakuba* | *D. melanogaster* | 6 / 12 Mya |
| 92 | *D. ananassae* | — | droAna3 | *D. melanogaster* | *D. pseudoobscura* | 25 / 40 Mya |
| 93 | *D. pseudoobscura* | — | dp4 | *D. persimilis* | *D. melanogaster* | 2 / 25 Mya |
| 94 | *D. virilis* | — | droVir3 | *D. mojavensis*, *D. grimshawi* | *D. melanogaster* | 20 / 40 Mya |

### Insects — other orders

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 95 | *Anopheles gambiae* | Malaria mosquito | anoGam3 | *Aedes aegypti* | *D. melanogaster* | 150 / 250 Mya |
| 96 | *Apis mellifera* | Honeybee | apiMel4 | *Bombus* spp. | *Nasonia* (jewel wasp) | 70 / 150 Mya |
| 97 | *Tribolium castaneum* | Red flour beetle | triCas2 | *D. melanogaster* | Honeybee | 250 / 300 Mya |
| 98 | *Bombyx mori* | Silkworm | — | *Manduca sexta* (tobacco hornworm) | *D. melanogaster* | 100 / 250 Mya |

### Worms & other invertebrates

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 99 | *C. elegans* | Nematode | ce11 | *C. briggsae*, *C. remanei* | *C. japonica* | 30 / 60 Mya |
| 100 | *C. briggsae* | — | cb4 | *C. remanei*, *C. elegans* | *C. japonica* | 20–30 / 60 Mya |
| 101 | *Strongylocentrotus purpuratus* | Purple sea urchin | strPur2 | Lancelet | *C. elegans* | 520 / 650 Mya |
| 102 | *Branchiostoma floridae* | Lancelet | braFlo1 | Sea urchin | Lamprey | 520 / 550 Mya |
| 103 | *Ciona intestinalis* | Sea squirt | ci3 | Lancelet | Lamprey | 520 / 550 Mya |
| 104 | *Aplysia californica* | Sea hare | aplCal1 | *C. elegans* | Sea urchin | 600 / 650 Mya |

### Plants

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 105 | *Arabidopsis thaliana* | Thale cress | araTha1† | *A. lyrata* | *Brassica rapa* | 5 / 20 Mya |
| 106 | *Brassica rapa* | Turnip / Chinese cabbage | braRap1† | *B. oleracea* | *Arabidopsis thaliana* | 4 / 20 Mya |
| 107 | *Oryza sativa* | Rice | — | *O. rufipogon*, *O. glaberrima* | *Brachypodium distachyon* | 1–2 / 50 Mya |
| 108 | *Zea mays* | Maize | — | *Sorghum bicolor* | *Oryza sativa* | 12 / 50 Mya |
| 109 | *Solanum lycopersicum* | Tomato | — | *S. tuberosum* (potato) | *Coffea canephora* | 8 / 80 Mya |
| 110 | *Glycine max* | Soybean | — | *Phaseolus vulgaris* (common bean) | *Arabidopsis thaliana* | 20 / 100 Mya |
| 111 | *Vitis vinifera* | Grape | — | *Coffea canephora* | *Arabidopsis thaliana* | 80 / 110 Mya |
| 112 | *Populus trichocarpa* | Poplar | — | *Salix* spp. (willow) | *Arabidopsis thaliana* | 8 / 100 Mya |

### Fungi

| # | Focal species | Common name | Assembly | Suggested inner outgroups | Suggested outer outgroup | Div. (inner / outer) |
|---|---------------|-------------|----------|---------------------------|--------------------------|----------------------|
| 113 | *Saccharomyces cerevisiae* | Baker's yeast | sacCer3 | *S. paradoxus* | *S. pombe* | 5 / 600 Mya |
| 114 | *Schizosaccharomyces pombe* | Fission yeast | — | *S. japonicus* | *Neurospora crassa* | 150 / 500 Mya |
| 115 | *Neurospora crassa* | Red bread mold | — | *N. tetrasperma* | *S. cerevisiae* | 30 / 500 Mya |

:::{note}
Assemblies marked with **†** are not part of the core UCSC Genome Browser but may be available through UCSC's GenArk or via Ensembl/Phytozome. Entries with **—** in the assembly column have public genome assemblies but typically no UCSC pairwise alignments; you will need to generate AXT alignments yourself (see below).
:::

### Generating your own alignments

For species not covered by UCSC's pre-computed pairwise alignments, you can
produce `net.axt.gz` files from any pair of genome assemblies:

```bash
# 1. Run lastz pairwise alignment
lastz target.2bit query.2bit \
    --format=axt --ambiguous=iupac \
    > raw.axt

# 2. Chain and net the alignment (UCSC Kent tools)
axtChain raw.axt target.2bit query.2bit chain.txt
chainSort chain.txt sorted.chain
chainNet sorted.chain target.chrom.sizes query.chrom.sizes target.net query.net
netToAxt target.net sorted.chain target.2bit query.2bit stdout \
    | axtSort stdin net.axt

# 3. Compress
gzip net.axt
```

The resulting `net.axt.gz` can be used directly as an alignment file in your
ancify config. See the [lastz documentation](https://lastz.github.io/lastz/)
and the [UCSC Kent utilities](https://hgdownload.soe.ucsc.edu/admin/exe/) for
installation and detailed options.
