#!/bin/bash -
#title          :NifMAP.sh
#description    :NifMAP - a bioinformatics pipline for analysing nifH amplicon data
#author         :Roey Angel
#date           :20210309
#version        :1.2
#usage          :./nifmap.sh <inputDir>
#notes          : Unlike V1.0, this pipeline is now independent and performs quality-filtering and OTU generation.
#               See Angel, R., Nepel, M., Panhölzl, C., Schmidt, H., Herbold, C. W., Eichorst, S. A., et al. (2018). Evaluation of primers targeting the diazotroph functional gene and development of NifMAP – a bioinformatics pipeline for analyzing nifH amplicon data. Front Microbiol 9, 703. doi:10.3389/fmicb.2018.00703.
#dependencies   : 1. HMMER - http://hmmer.org/
#                 2. FrameBot - https://github.com/rdpstaff/Framebot
#                 3. seqmagick - https://fhcrc.github.io/seqmagick/
#                 4. CART - https://wwwzehr.pmc.ucsc.edu/CART_model_public/
#                 5. R - https://www.r-project.org/
#bash_version   :4.3.48(1)-release
#============================================================================

# Optimized pipeline for moss nifH analyses. Dennis Escolástico-Ortiz. 2022

module load usearch/9.2.64_i86linux32
module load hmmer/3.3

module load mafft/7.453
module load python/3.7 seqmagick/0.8.4




eval "$(conda shell.bash hook)" # this is needed to be able to use conda activate

# Define variables
inputDir="./filtered" #$1 # input fasta file containing merged MiSeq reads
HOMEFOLDER=`pwd` # base library
FRAMEBOTPATH="../RDPTools"
WORKFOLDER=${HOMEFOLDER}/nifH_work
RESULTSFOLDER=${HOMEFOLDER}/nifH_results
RESOURCEFOLDER=${HOMEFOLDER}/Resources
INPUTREFERENCEALIGNMENT=cluster.rep.nifH-chL-bchX.0.9_noDot.fasta
INPUTREFERENCETREE=RAxML_bipartitions.MultipleOriginal.tree
id=`echo 0.94 | bc` # OTU clustering, 1.0 = zOTUS, 0.97 = 97%
cores=4 # number of cores to use
MINLEN=200 # merged reads with a smaller size will be removed
DB="/biodata/blast/refseq_protein"


# 8. Classify sequences using CART
module load hmmer/3.3
cat ${RESOURCEFOLDER}/AztVine8_AA.fasta ${RESULTSFOLDER}/nifH_corr_prot.fasta > ${WORKFOLDER}/nifH_corr_prot_4classification.fa
hmmalign --amino ${RESOURCEFOLDER}/Zehr_2014_1812genomes_nifH_AA_noGaps.hmm ${WORKFOLDER}/nifH_corr_prot_4classification.fa > ${WORKFOLDER}/nifH_corr_prot_hmm.sth

module load python/3.7 seqmagick/0.8.4 ### DO NOT MOVE OR ERASE THIS LINE, OTHERWISE WILL NOT WORK
seqmagick convert ${WORKFOLDER}/nifH_corr_prot_hmm.sth  ${WORKFOLDER}/nifH_corr_prot_hmmAln.fasta

python3 ${RESOURCEFOLDER}/NifH_Clusters.py ${WORKFOLDER}/nifH_corr_prot_hmmAln 1
grep '>' ${WORKFOLDER}/nifH_corr_prot_hmmAln_Clusters.fasta |  sed 's/>//' | awk '{print $1, $5, $8}' > ${RESULTSFOLDER}/nifH_corr_prot.zehr.classification

# 9. Classify sequences using BLASTP

# The following two lines allow downloading the taxa ids
#wget http://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz ./
#tar -zxvf taxdb.tar.gz

blastp -query ${RESULTSFOLDER}/nifH_corr_prot.fasta -db $DB -num_threads 4 -max_target_seqs 1 -out ${RESULTSFOLDER}/nifH_blastp.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sblastnames sskingdoms stitle"

# 8. Return results and clean up
#return results
#cp ${RESULTSFOLDER} $HOMEFOLDER/

#cleanup
#rm -r ${RESULTSFOLDER}/tmp.log ${RESULTSFOLDER}/unique.fa ${RESULTSFOLDER}/unoise3.txt
#rm -r ${WORKFOLDER}
#rm taxdb*
