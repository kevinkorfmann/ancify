# ancify terminal demo — render with: vhs ancify_demo.tape
# Install vhs: https://github.com/charmbracelet/vhs

Output docs/_static/ancify_demo.gif
Output docs/_static/ancify_demo.webm

Set Theme "Tokyo Night"
Set FontFamily "Fira Code"
Set FontSize 16
Set Width 900
Set Height 550
Set Padding 20
Set TypingSpeed 40ms
Set CursorBlink false

Sleep 500ms

Type "pip install ."
Enter
Sleep 1.5s
Type "ancify init -o config.yaml"
Enter
Sleep 800ms

Type "# edit config.yaml: species, alignments, output_dir"
Enter
Sleep 1.5s

Type "ancify run -c config.yaml"
Enter
Sleep 500ms

# Simulated output
Type "Phase 1: projecting 4 outgroups × 23 chromosomes..."
Enter
Sleep 2s
Type "  ✓ bonobo    ✓ chimp    ✓ gorilla    ✓ macaque"
Enter
Sleep 1s
Type "Phase 2: calling ancestral alleles..."
Enter
Sleep 1.5s
Type "Done. Output: ./ancestral_calls/"
Enter
Sleep 2s

Type "head -c 80 ancestral_calls/chr22.fa"
Enter
Sleep 500ms
Type ">chr22"
Enter
Type "ACGTACGTacgtACGTNNNNACGTacgtACGTnACGTACGTacgtACGT"
Enter
Sleep 2s

Type "# ACGT = high confidence   acgt = low   n = unresolved   N = missing"
Enter
Sleep 3s

Type "# ancestral alleles, decoded."
Enter
Sleep 3s
