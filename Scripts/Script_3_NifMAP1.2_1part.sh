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

if [ -d ${WORKFOLDER} ]
then
 rm -r ${WORKFOLDER}
fi
if [ -d ${RESULTSFOLDER} ]
then
 rm -r ${RESULTSFOLDER}
fi
mkdir ${WORKFOLDER}
mkdir ${RESULTSFOLDER}

echo "NifMAP1.2"> ${RESULTSFOLDER}/NifMAP.log



# 1. Merge MiSeq aplicon reads
echo "usearch -fastq_mergepairs">> ${RESULTSFOLDER}/NifMAP.log
cd $inputDir

#for gzfile in *.fastq.gz 
#do
#	gunzip $gzfile
#done

usearch9.2.64_i86linux32 -fastq_mergepairs *_R1_001.fastq -fastqout merged_reads.fq  -relabel @  -threads $cores --log ${RESULTSFOLDER}/tmp.log
cd ../
cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
echo "--">> ${RESULTSFOLDER}/NifMAP.log

# 2. Quality filter merged reads
echo "usearch -fastq_filter">> ${RESULTSFOLDER}/NifMAP.log
#for inputReads in ${inputDir}/*.fastq; do
#    usearch -fastq_filter ${inputReads} -fastq_maxee 1.0 -fastq_minlen 200 -relabel @ -fastaout ${WORKFOLDER}/${inputReads%.*}.fasta --log ${RESULTSFOLDER}/tmp.log;
#    cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
#done
usearch9.2.64_i86linux32 -fastq_filter ${inputDir}/merged_reads.fq -fastq_maxee 1.0 -fastq_minlen $MINLEN -fastaout ${WORKFOLDER}/merged_reads.fa -threads $cores --log ${RESULTSFOLDER}/tmp.log
cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
echo "--" >> ${RESULTSFOLDER}/NifMAP.log

# 3. Filter sequences using nifH HMM
# Screen merged reads using HMMer
module load hmmer/3.3
echo "Screen merged reads using HMMer">> ${RESULTSFOLDER}/NifMAP.log
hmmsearch --cpu $cores --domtblout ${WORKFOLDER}/hmmOut1.out -o ${WORKFOLDER}/junk ${RESOURCEFOLDER}/hmm_nuc_1160_nifH.hmm ${WORKFOLDER}/merged_reads.fa
awk '{print $1}' ${WORKFOLDER}/hmmOut1.out | grep -v "#" > ${WORKFOLDER}/acceptable_hits
grep ">" ${WORKFOLDER}/merged_reads.fa | grep -v -F -f ${WORKFOLDER}/acceptable_hits >${WORKFOLDER}/shitHits
totalUnique=`grep ">" ${WORKFOLDER}/merged_reads.fa | wc -l`
totalAccepted=`cat ${WORKFOLDER}/acceptable_hits | wc -l`
totalRemoved=`cat ${WORKFOLDER}/shitHits | wc -l`
echo "hmmscreen of nifH removed ${totalRemoved} sequences out of ${totalUnique} sequences. ${totalAccepted} unique sequences retained" >> ${RESULTSFOLDER}/NifMAP.log
mv ${WORKFOLDER}/merged_reads.fa ${WORKFOLDER}/merged_reads_prehmm.fa
awk 'BEGIN{FS="\n";RS=">"};NR>1{print(">"$1);for(i=2;i<=NF;i++){printf($i)};print("")}' ${WORKFOLDER}//merged_reads_prehmm.fa | grep -A 1 -F -f ${WORKFOLDER}/acceptable_hits | grep -v "^\-\-$" >${RESULTSFOLDER}/merged_reads.fa
echo "--">> ${RESULTSFOLDER}/NifMAP.log
