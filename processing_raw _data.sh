#!/bin/bash

source /etc/profile.d/modules.sh
module load mothur/v1.41.3

MOTHUR=/opt/mothur/v1.41.3/bin/mothur
USEARCH=usearch11.0.667
VSEARCH=vsearch-2.12.0

PROCS=$(nproc)
MAXDIFFS=2
MAXEE=0.5
MINSIZE=1

export OMP_NUM_THREADS=$PROCS

DB_DIR="/data/skin_microbiome/db/"

RDP_DB="/data/skin_microbiome/db/rdp_gold_v9.fasta"

MOTHUR_PCR=$DB_DIR"silva.123/silva.v123.align"

# RDP v18
FASTA_RDP=$DB_DIR"rdp_16s_v18.fa"
SINTAX_RDP=$DB_DIR"rdp_16s_v18.udb"

echo "start"
date;

echo "Processors " $PROCS
echo "Maxdiffs " $MAXDIFFS
echo "Maxee " $MAXEE
echo "Minsize " $MINSIZE

echo "OMP_NUM_THREADS " $OMP_NUM_THREADS

cd ./00_fastq

############
#
# merge fastq files (forward/reverse read)
#

# use usearch instead of vsearch for this dataset
#ls *1.fastq |sed 's/\./ /g' | uniq | sort -u| awk -v usearch=$USEARCH -v maxdiffs=$MAXDIFFS -v procs=$PROCS '{print ""usearch" -fastq_mergepairs "$1".1.fastq -reverse "$1".2.fastq -fastq_maxdiffs "maxdiffs" -fastqout "$1".merged.fastq -threads "procs""}' > usearch.mergepairs.sh
ls *R1_001.fastq.gz | sed 's/\_/ /g' | uniq | sort -u| awk -v usearch=$VSEARCH -v maxdiffs=$MAXDIFFS -v procs=$PROCS '{print ""usearch" -fastq_mergepairs "$1"_"$2"_L001_R1_001.fastq.gz -reverse "$1"_"$2"_L001_R2_001.fastq.gz -fastq_maxdiffs "maxdiffs" -fastqout "$1".merged.fastq -threads "procs" --fastq_minovlen 100"}' > vsearch.mergepairs.sh

bash vsearch.mergepairs.sh

mkdir ../01_merged
mv *merged* ../01_merged

gzip *.fastq

cd ../01_merged

############
#
# rename reads & filter
#

ls *merged.fastq |sed 's/\./ /g' | uniq | sort -u | awk -v vsearch=$VSEARCH -v maxee=$MAXEE -v procs=$PROCS '{print ""vsearch" -fastq_filter "$1"."$2".fastq -fastq_maxee "maxee" -fastaout "$1".filtered.fasta -relabel "$1". -threads "procs""}' > vsearch.filtered.sh
bash vsearch.filtered.sh

mkdir ../02_filtered
mv *.filtered.fasta ../02_filtered

gzip *.fastq

cd ../02_filtered

############
#
# discover chimeras
#

ls *.filtered.fasta |sed 's/\./ /g' | uniq |sort -u| awk -v vsearch=$VSEARCH -v rdp_db=$RDP_DB -v procs=$PROCS '{print ""vsearch" -uchime_ref "$1".filtered.fasta -db "rdp_db" -nonchimeras "$1".noChime.fasta -strand plus -threads "procs""}' > uchime.db.chimera.sh

bash uchime.db.chimera.sh

mkdir ../03_noChime
mv *.noChime.fasta ../03_noChime

gzip *.fasta

#
#
############

############
#
# 'merge' all subsampled reads for OTU identification
#

cd ../
cat 03_noChime/*.fasta > all.fa
grep -c ">" < all.fa

gzip 03_noChime/*.fasta
gzip all.fa

############
#
# statistics about yield per filtering step
#
# raw data
zgrep -c ^ 00_fastq/*R1*.gz | sed 's/\:/ /g' | awk '{print $1," ",$2/4}' > data.00.raw.pairs.txt
# merged reads
zgrep -c ^ 01_merged/*merged*.gz | sed 's/\:/ /g' | awk '{print $1," ",$2/4}' > data.01.merged.txt
# quality filtered
zgrep -c ">" 02_filtered/*filtered*.gz | sed 's/[\:]/ /g' | awk '{print $1," ",$2}' > data.02.filtered.txt
# chimeras removed
zgrep -c ">" 03_noChime/*noChime*.gz | sed 's/[\:]/ /g' | awk '{print $1," ",$2}' > data.03.noChime.txt

echo -e "File Raw Merged Filter Chimerafree" > stats.data.txt

paste data.00.raw.pairs.txt data.01.merged.txt data.02.filtered.txt data.03.noChime.txt  | awk '{print $1,$2,$4,$6,$8}' >> stats.data.txt

rm data.00.raw.pairs.txt data.01.merged.txt data.02.filtered.txt data.03.noChime.txt
#
#############

############
#
# finding unique sequences (Dereplication)
#
# and sort them by size
#

$VSEARCH --derep_fulllength all.fa.gz --relabel Uniq --sizeout --output uniques.fa --threads $PROCS --log log.vsearch.derep.txt --minuniquesize 2 --fasta_width 0

############
#
# Clustering/OTUs
#

$USEARCH -cluster_otus uniques.fa -minsize $MINSIZE -otus otus.fa -relabel Otu -log log.usearch.cluster_otus.txt

# $USEARCH -unoise3 uniques.fa -zotus zotus.fa -tabbedout unoise3.txt

gzip uniques.fa

############
#
# assign taxonomy using SINTAX (usearch)
#

$USEARCH -makeudb_usearch $FASTA_RDP -output $SINTAX_RDP

$USEARCH -sintax otus.fa -db $SINTAX_RDP -tabbedout otus.sintax.txt -strand both -sintax_cutoff 0.8 -threads $PROCS

_sintax_to_Fasta.pl otus.sintax.txt otus.fa

############
#
# Generating an OTU table
#

$VSEARCH --usearch_global all.fa.gz --db sintax.otus.fa --strand plus --id 0.97 --biomout otu_table.sintax.biom --threads $PROCS -log log.vsearch.usearch_global.otus.txt

echo "end usearch"
date;

##########
#
# Align sequences & FastTreeMP
#

#### MOTHUR #####
#
# align otus
$MOTHUR "#set.logfile(name=log.mothur.txt, append=T); align.seqs(fasta=otus.fa, reference=$MOTHUR_PCR, processors=$PROCS)"
$MOTHUR "#set.logfile(name=log.mothur.txt, append=T); summary.seqs(fasta=otus.align, processors=$PROCS)"

#export OMP_NUM_THREADS=$PROCS
FastTreeMP -gtr -nt -gamma otus.align > ./tree.otus.tree


