# FastProteinExonerate
Match protein sequences to a genome and predict genes in the matching genome regions. Using `pblat` and running `exonerate` only on the smaller matching section of the genome improves the speed of the generally slower process of running `exonerate` genome-wide.

Takes a protein and DNA input file (usually scaffolds), `plat`s the proteins against the DNA sequences, then predicts genes in the matched sections (+-500nt) using `exconerate protein2genome`

# Usage

`FastProteinExonerate_v220221.sh <protein file> <DNA file> <n cores> <maxIntron>`

By default, the script tries to find a conda executable in

`CONDASH=/data/miniconda3/etc/profile.d/conda.sh`

If your `conda.sh` is in a different location, edit the path to `$CONDASH` in the script.

# Output
All output files will be in a new folder called `protExon`. If this folder exists it will be overwritten!

1. `cleaned_proteins.fasta` (Basic clean up of input `.fasta` file, a `.fasta` file in one-line format)
2. `protein_out.psl` (Output of `pblat`)
3. `pblat.log` (`pblat` log file)
4. `pblat.err` (`pblat` error file)
5. `best_hits_protein_out.psl` (Filtered `pblat` output file, only take the best hit for each input protein query)
6. `coord.info.tsv` (a `.tsv` file with genome coordinate info +- 500nt of the matched region, one line for each matched protein sequence)
7. `match_coord.bed` (matched coordinates in `.bed` format)
8. `match_sections.fasta` (nucleotide sequences of regions)
9. `run.sh` (The actual script that does all the work. It is created at runtime and will be quite large as it contains sequence data, not recommended to `less`/`more`/`cat` it)
10. `run.log` (Log `STDOUT` file of the run)
11. `run.err` (Error `STDERR` file of the run)
12. `final.gff` (Output: Predicted genes in GFF format)
13. `final.proteins.fa` (Output: translated protein sequences)
14. `final.cds.fa` (Output: coding sequences (CDSs))

# Dependencies

The script attempts to create a `conda` environment `proteinexonerate` which will install the following dependencies. If `proteinexonerate` exists, it will activate the existing environment.

1. `pblat`
2. `bedtools`
3. `exonerate`
4. `gffread`
