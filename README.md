Process explanation:

Part 1 - Get enhancer sequences and target gene names:

1. Read elite GeneHancer bed file, line by line
2. For each line, extract coordinates of the enhancer and the name of the target gene (in Ensembl)
3. Look up the Hugo name of the target gene and its biotype/description in genes.ENSG.tbl
4. Extract the enhancer sequence from referenceGenome matching the coordinates
5. Output this intermediary file
#
Part 2 - Get TF binding site motifs and apply them to enhancers
1. Read in Hocomoco PFM matrices one at a time
2. For each matrix, generate a Biopython Motif
3. For each enhancer region from 1) find which BS Motifs match there (on + and - strands)
4. Calculate how many binding sites on average match within 1 enhancer region - to be informed on the supernode
#
Part 3 - Use co-expression data to find TFBS clusters 
1. Use GTEX co-expression data to figure out how TFs that match within the same enhancer region
   regulate the transcription of the target gene - to be informed on the supernode
#
Files expected to be present:
1. **GRCh38.primary_assembly.genome.fa** - reference genome ch38
2. **genes.ENSG.tbl** - gene names in Ensembl and Hugo forms, along with biotype description
3. **elite_ensg_enhan_fused_ensembl_prom_500b.hg38.bed** - elite enhancer coordinates from GeneHancer
4. **elite_enhancer_sequences.fa** - elite enhancers sequences (generated)
(You can generate it by running the following command:
`bedtools getfasta -fi GRCh38.primary_assembly.genome.fa
  -bed elite_ensg_enhan_fused_ensembl_prom_500b.hg38.bed -fo elite_enhancer_sequences.fa`