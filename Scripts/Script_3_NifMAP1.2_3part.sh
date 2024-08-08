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
module load ncbiblast/2.11.0
module load mafft/7.453
module load RAxML/8.2.9

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

# 6. Make an OTU table based on corrected seqs
echo "usearch -usearch_global ">> ${RESULTSFOLDER}/NifMAP.log
usearch9.2.64_i86linux32 -usearch_global ${RESULTSFOLDER}/merged_reads.fa -db ${RESULTSFOLDER}/nifH_corr_nucl.fasta  -strand plus -id 0.97 -threads $cores -otutabout ${RESULTSFOLDER}/otu_table.txt --log ${RESULTSFOLDER}/tmp.log
cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
echo "--">> ${RESULTSFOLDER}/NifMAP.log

# 7. Filter out homologous genes (bchL; chlL; bchX; parA) using HMM:
# Corrected AA sequences are in nifH_corr_prot.fasta
# Screen with hmm to identify all hits
hmmscan --cpu $cores  --domtblout ${WORKFOLDER}/hmmOut2.out -o ${WORKFOLDER}/full_output_NifH_ChlL_bchX ${RESOURCEFOLDER}/nifH_chlL_bchX.hmm ${WORKFOLDER}/nifH_corr_prot.fasta
#Rscript --vanilla < ${HOMEFOLDER}/nifH_bch_hmmEvaluation.R ${WORKFOLDER}/hmmOut2.out
#mv nifH_bch_hmmEvaluation.pdf ${RESULTSFOLDER}
cat ${WORKFOLDER}/hmmOut2.out | awk 'NR>3{if($8>bitarray[$4]){bitarray[$4]=$8;outArray[$4]=$1"\t"$4}}END{for(entry in outArray){print outArray[entry]}}' > ${WORKFOLDER}/assignments.txt
cp ${WORKFOLDER}/hmmOut2.out ${RESULTSFOLDER}/nifH_bch_hmmEvaluation.hmm.out
grep "nifH" ${WORKFOLDER}/assignments.txt | awk '{print $2}' | sort > ${WORKFOLDER}/acceptable_hits
grep ">" ${RESULTSFOLDER}/nifH_corr_nucl.fasta | awk '{print $1}' | grep -v -F -f ${WORKFOLDER}/acceptable_hits | sed 's/>//'>${WORKFOLDER}/shitHits
totalOTUs=`grep ">" ${RESULTSFOLDER}/nifH_corr_nucl.fasta | wc -l`
totalAccepted=`cat ${WORKFOLDER}/acceptable_hits | wc -l`
totalRemoved=`cat ${WORKFOLDER}/shitHits | wc -l`
echo "FRAMEBOT and hmmscreen of nifH removed ${totalRemoved} sequences out of ${totalOTUs} ASVs. ${totalAccepted} ASV reps retained">>${RESULTSFOLDER}/NifMAP.log

#mv ${RESULTSFOLDER}/${inputReads} ${WORKFOLDER}/${inputReads}
cat ${RESULTSFOLDER}/nifH_corr_nucl.fasta | sed 's/ //g' | awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/nifH_corr_nucl_only_nifH.fasta
cat ${WORKFOLDER}/nifH_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/nifH_corr_prot.fasta
cat ${WORKFOLDER}/nifH_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/shitHits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/nifH_rej_prot.fasta

# Make an OTU table based on filter corrected seqs
echo "usearch -usearch_global ">> ${RESULTSFOLDER}/NifMAP.log
#usearch9.2.64_i86linux32 -usearch_global ${RESULTSFOLDER}/merged_reads.fa -db ${RESULTSFOLDER}/nifH_corr_nucl_only_nifH.fasta -id 1.0 -strand plus -threads $cores -otutabout ${RESULTSFOLDER}/otu_table.txt
usearch9.2.64_i86linux32 -otutab ${RESULTSFOLDER}/merged_reads.fa -otus ${RESULTSFOLDER}/nifH_corr_nucl_only_nifH.fasta -otutabout ${RESULTSFOLDER}/otu_table.txt -mapout ${WORKFOLDER}/map.txt -threads $cores

# 7. Place OTUs on a reference tree using RAxML-EPA
#RAxML-EPA
mkdir ${RESULTSFOLDER}/EPA
cp ${RESULTSFOLDER}/nifH_corr_prot.fasta ${RESULTSFOLDER}/EPA/
cp ${RESOURCEFOLDER}/${INPUTREFERENCEALIGNMENT} ${RESULTSFOLDER}/EPA/
cp ${RESOURCEFOLDER}/${INPUTREFERENCETREE} ${RESULTSFOLDER}/EPA/

cd ${RESULTSFOLDER}/EPA
mafft --add nifH_corr_prot.fasta --thread 16 ${INPUTREFERENCEALIGNMENT} > completeAlignment.fa
cat completeAlignment.fa | awk '{print $1}' | sed 's/|$//'> RAxML_compatible.fa
raxmlHPC-PTHREADS-SSE3 -f v -T $cores -s RAxML_compatible.fa -m PROTCATJTT -t ${INPUTREFERENCETREE} -n EPAplaced -p 123


