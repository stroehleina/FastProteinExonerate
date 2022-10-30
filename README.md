# FastProteinExonerate
Match protein sequences to a genome and predict genes in the matching genome regions. Using `pblat` and running `exonerate` only on the smaller matching section of the genome improves the speed of the generally slower process of running `exonerate` genome-wide.

Takes a protein and DNA input file (usually scaffolds), `PBLAT`s the proteins against the DNA sequences, then predicts genes in the matched sections (+-500nt) using `exconerate protein2genome`

# Usage

`FastProteinExonerate_v220221.sh <protein file> <DNA file> <n cores> <maxIntron>`

# Dependencies

The script attempts to create a `conda` environment `proteinexonerate` which will install the following dependencies.

1. `pblat`
2. `bedtools`
3. `exonerate`
4. `gffread`
