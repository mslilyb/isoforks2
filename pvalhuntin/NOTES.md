# Steps taken

1. Download .gff3 and .fa files from WormBase. Using bioproject PRNA0013758.
    - Filenames:
        - c_elegans.PRJNA13758.WS282.genomic.fa.gz - sequence
        - c_elegans.PRJNA13758.WS282.annotations.gff3.gz - features

2. Strip gff3 to only entries including WormBase genes and RNASeq data.
Command used:
```bash
gunzip -c ~/Downloads/c_elegans.PRJNA13758.WS282.annotations.gff3.gz | grep -E 'WormBase|RNASeq' > ws282.gff3
```
3. Counted number of genes total in ws282 .gff3
Command used:
```bash
cat ws282.gff3 | grep -E 'WormBase[[:blank:]]{1,}gene' | uniq | wc -l
```
4. Verified counts of unique genes from stripped gff3 was same as full file.
Command used:
```bash
gunzip -c ~/Downloads/c_elegans.PRJNA13758.WS282.annotations.gff3.gz | grep -E 'WormBase[[:blank:]]{1,}gene' | uniq | wc -l
```
Results: 47651 genes found in both.

5. Adapted `icountvar.py` from `datacore2024/genome_celegans`.
    - original script found multi intron, single isoform genes in genome.
    - intermediate script finds all gene ids and their intron counts.
    
6. Used intermediate script to determine how many genes fit our criteria across entire genome
    - detailed output is in textfile `rawcounts`
Command used:
```bash
cat rawcounts | uniq | wc -l
```
Results: 14374 genes, 30.2% of genes in the .gff3 file.

7. Final script calculates taxicab distance of intron counts from the mean.
    - Methods:
        - Find average counts across all introns in gene
        - Convert counts into ratios by dividing by the count mean
        - find distance from 'ground truth', i.e. count ratios - 1.
            - 1 is used as the ground truth because a single isoform should have only those introns, and there should be equal counts of all (i.e. a probability/frequency of 1.00)

Results:
1pct average distance: 1.1341249130710245
Whole Genome average distance: 1.1096404743864219
