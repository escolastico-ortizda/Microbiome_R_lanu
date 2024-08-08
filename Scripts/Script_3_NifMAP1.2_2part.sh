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

# Code for moss nifH analyses.


module load usearch/9.2.64_i86linux32
module load hmmer/3.3

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

# 4. Generate OTUs
# Dereplicate sequences for OTU table
echo "usearch -fastx_uniques">> ${RESULTSFOLDER}/NifMAP.log
echo "Singletons will not be allowed to generate OTUs">> ${RESULTSFOLDER}/NifMAP.log
usearch9.2.64_i86linux32 -fastx_uniques ${RESULTSFOLDER}/merged_reads.fa -minuniquesize 2 -fastaout ${RESULTSFOLDER}/unique.fa -threads $cores --log ${RESULTSFOLDER}/tmp.log # -sizeout
cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
echo "--">> ${RESULTSFOLDER}/NifMAP.log

# Make zOTUs or OTUs
if [[ (( $id == 1.0 )) ]]; then
    echo "usearch -unoise3 ">> ${RESULTSFOLDER}/NifMAP.log
    usearch9.2.64_i86linux32 -sortbylength ${RESULTSFOLDER}/unique.fa -fastaout ${WORKFOLDER}/unique_sorted.fa
    usearch9.2.64_i86linux32 -unoise3 ${WORKFOLDER}/unique_sorted.fa -zotus ${WORKFOLDER}/otus.fa -tabbedout ${RESULTSFOLDER}/unoise3.txt -threads $cores --log ${RESULTSFOLDER}/tmp.log
    cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
    echo "--">> ${RESULTSFOLDER}/NifMAP.log
else
    echo "usearch -cluster_fast">> ${RESULTSFOLDER}/NifMAP.log
    usearch9.2.64_i86linux32 -sortbylength ${RESULTSFOLDER}/unique.fa -fastaout ${WORKFOLDER}/unique_sorted.fa
    usearch9.2.64_i86linux32 -cluster_fast ${WORKFOLDER}/unique_sorted.fa -id $id -centroids ${WORKFOLDER}/otus.fa -uc ${WORKFOLDER}/clusters.uc -threads $cores --log ${RESULTSFOLDER}/tmp.log
    cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
    echo "--">> ${RESULTSFOLDER}/NifMAP.log
fi

# 5.Translate and correct OTU representatives using Framebot
#eVal_chL=1e-50
#score_chL=150
#eVal_bchX=1e-50
#score_bchX=150

java -jar ${FRAMEBOTPATH}/FrameBot.jar framebot -N -l 30 -i 0.4 -o ${WORKFOLDER}/nifH ${FRAMEBOTPATH}/Framebot/refset/nifh_prot_ref.fasta ${WORKFOLDER}/otus.fa
cp ${WORKFOLDER}/nifH_corr_nucl.fasta ${RESULTSFOLDER}

