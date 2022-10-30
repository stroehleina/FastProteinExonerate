#!/bin/bash
#Takes a protein and DNA input file, pblats the proteins against the DNA sequences, then predicts genes in the matched sections using exconerate protein2genome

CONDASH=/data/miniconda3/etc/profile.d/conda.sh

if [ ! -f $CONDASH ]
  then
    echo "Conda executable does not exist in "$CONDASH" Please adjust set $CONDASH and rerun. Exiting..."
    exit
  else
    source $CONDASH
fi

ENVIRO=$(conda env list | grep -v "#" | awk '{print $1}' | sed '/^$/d' | grep "proteinexonerate")

if [ "$ENVIRO" == "proteinexonerate" ]
  then
    echo "Conda environment \"proteinexonerate\" exist. Activating environment"
    conda activate proteinexonerate
  else
    echo "Creating conda environment \"proteinexonerate\" "
    conda create -n proteinexonerate -c bioconda pblat bedtools exonerate gffread
    conda activate proteinexonerate
fi

if [ "$#" -lt "4" ]
  then
    echo "Usage: FastProteinExonerate_v220221.sh <protein file> <DNA file> <n cores> <maxIntron>"
    exit
fi

if [[ "$1" == "-h" ]]
  then
    echo "Usage: FastProteinExonerate_v220221.sh <protein file> <DNA file> <n cores> <maxIntron>"
    exit
fi


if [[ "$1" == "--help" ]]
  then
    echo "Usage: FastProteinExonerate_v220221.sh <protein file> <DNA file> <n cores> <maxIntron>"
    exit
fi



echo "$(date): Starting run"
echo "$(date): Running pblat and Exonerate for protein data set $1 against DNA sequences $2 using $3 cores and maxIntron length $4"
if [ -d "protExon" ]
  then
    echo "$(date): Removing existing folder protExon"
    rm -rf protExon
fi

mkdir protExon && cd protExon

#clean input protein data set (read in from command line via $1)"
echo "$(date): Cleaning protein input $1"
awk '{print $1}' ../$1 | sed 's/|/_/g' | awk 'BEGIN{getline; print $0; seq="";}{if($0 ~ "^>"){print seq; seq=""; print $0;}else{seq=seq""$0}}END{print seq}' > cleaned_proteins.fasta

#run pblat with number of threads $3
echo "$(date): Run pblat with $3 cores: $1 against $2"
pblat -threads=$3 -q=prot -t=dnax -fine -maxIntron=$4 ../$2 cleaned_proteins.fasta protein_out.psl 1>pblat.log 2>pblat.err 

#Pick top hit only
echo "$(date): Picking best pblat hit for each protein"
awk -F"\t" '{if($0 ~ /^psLayout/){for(i=1; i < 5; i++){getline;}}else{print $0}}' protein_out.psl | sort -t$'\t' -k10,10 -k1,1nr -k15,15nr | sort -u -k10,10 --merge > best_hits_protein_out.psl

#create coordinate file with +- 500nt distances
echo "$(date): Adjusting match coordinates and creating lookup files"
cut -f 10,14-17 best_hits_protein_out.psl | awk -F"\t" '{printf $1"\t"$2"\t"$3; if($4-500 < 0){printf "\t0"}else{printf "\t"$4-500}; if($5+500 > $3-1){print "\t"$3-1}else{print "\t"$5+500};}' | sed 's/|/_/g' > coord.info.tsv
#awk -F"\t" '{if(){}; }' > tmp && mv tmp coord.info.tsv

#cut -f 10,14-17 best_hits_protein_out.psl | awk -F"\t" '{printf $0; if($4-500 < 0){printf "\t0"}else{printf "\t"$4-500}; if($5+500 > $3){print "\t"$3}else{print "\t"$5+500};}' | sed 's/|/_/g' > coord.info.tsv
cut -f2,4,5 coord.info.tsv | sort | uniq > match_coord.bed

#extract scaffold regions using bedtools
echo "$(date): Using bedtools to extract the match regions from $2"
bedtools getfasta -fi ../$2 -bed match_coord.bed | sed 's/:/__/' > match_sections.fasta

#prepare the script that creates input files, runs exonerate, cleans up gff etc.

echo "$(date): Creating run.sh script. This script will be very large and contains sequence data, do not attempt to less/more/cat it"
#IF YOUR conda.sh IS SOMEWHERE ELSE MODIFY THE NEXT LINE ACCORDINGLY
echo "source $CONDASH" > run.sh
echo "conda activate proteinexonerate" >> run.sh
echo "echo \"\$(date): Start run.sh\"" >> run.sh
awk -F"\t" 'BEGIN{
        while(getline < "match_sections.fasta"){
                a=$1; getline < "match_sections.fasta"; seq[a]=$1;
        };

        while(getline < "cleaned_proteins.fasta"){
                b=$1; getline < "cleaned_proteins.fasta"; prot[b]=$1;
        };
}{
        q=$1; t=$2; len=$3; sn=$4; en=$5;
        q_name=q".fasta"
        t_comb=t"__"sn"-"en
        t_name=t_comb".fasta"
        print("echo \">"t_comb"\" > "t_name);
        print("echo \""seq[">"t_comb]"\" >>  "t_name);
        print("echo \">"q"\" > "q_name);
        print("echo \""prot[">"q]"\" >>  "q_name);
        print("exonerate --minintron 20 --model protein2genome -S FALSE --showtargetgff TRUE --showvulgar FALSE --showalignment no --singlepass FALSE --bestn 1 -q "q_name" -t "t_name" > "q".gff")
        #print("rm "q_name" "t_name)
        print("sed '\''s/\\tgene_id 0 /\\tID="q"/'\'' "q".gff | awk -F\"\\t\" '\''{if($3 == \"gene\"){print $0; gsub(/gene/,\"mRNA\",$3); for(i=1; i<NF; i++){printf $i\"\\t\"}; print $NF}else{print $0}}'\'' | awk -F\"\\t\" '\''{if($9 !~ /^ID=/ && $9 !~ /^Parent=/){print $0\"\\tParent="q";\"$9}else{print $0\"\\t\"$9}}'\'' | cut -f1-8,10 > tmp; mv tmp "q".gff")
        print("grep -v -e \"^-- completed exonerate analysis$\" -e \"^#\" -e \"^Command line:\" -e \"^Hostname:\" "q".gff | awk -F\"\\t\" '\''{if($3 != \"similarity\"){$4=$4+"sn"; $5=$5+"sn"; printf \""t"\\t\"; for(i=2; i<NF; i++){printf $i\"\\t\"}; print $NF}}'\'' >> final.gff")
        #print("rm "q".gff")

}' coord.info.tsv >> run.sh

#NORUN #print("sed '\''s/\\tgene_id 0 /\\tID="q"/'\'' "q".gff | awk -F\"\\t\" '\''{if($3 == \"gene\"){print $0; gsub(/gene/,\"mRNA\",$3); gsub(/^ID=/,\"Parent=\",$9); for(i=1; i<NF; i++){printf $i\"\\t\"}; print $NF}else{print $0}}'\'' | awk -F\"\\t\" '\''{if($9 !~ /^ID=/ && $9 !~ /^Parent=/){print $0\"\\tParent="q";\"$9}else{print $0\"\\t\"$9}}'\'' | cut -f1-8,10 > tmp; mv tmp "q".gff")

echo "echo \"\$(date): End run.sh\"" >> run.sh
echo "$(date): run.sh created"
echo "$(date): Size of run sh: $(du -h run.sh | awk '{print $1}')"

echo "$(date): Executing run.sh. Log file: run.log. Error file: err.log"
bash run.sh 1> run.log 2> run.err
awk -F"\t" '{if($3 == "mRNA" || $3 == "cds" || $3 == "exon"){print $0}}' final.gff | sed 's/\tcds\t/\tCDS\t/' | awk -F"\t" '{if($9 ~ /^ID=/){print $0"\t"$9}else{split($9,attr,"[;=]"); count[$3""attr[2]]++; print $0"\tID="attr[2]"_"$3"_"count[$3""attr[2]]"; "$9}}' | cut -f 1-8,10 > tmp; mv tmp final.gff
#awk -F"\t" '{if($3 == "gene" || $3 == "mRNA" || $3 == "cds" || $3 == "exon"){print $0}}' final.gff | sed 's/\tcds\t/\tCDS\t/' > tmp; mv tmp final.gff

echo "$(date): Creating CDS and protein files"
gffread final.gff -g ../$2 -y final.proteins.fa
gffread final.gff -g ../$2 -x final.cds.fa
echo "$(date): Done running run.sh. Analysis complete. Output in final.gff, final.cds.fa, final.proteins.fa"
