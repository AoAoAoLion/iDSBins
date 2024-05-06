#!/bin/bash
##################################################################
# Author: Xin Wang                                                   
# Email: xin.wang@childrens.harvard.edu                            
# Copyright (c) 2021 Dr. Kaifu lab                                   
# PI: Kaifu Chen                                                   
# Description: 
#  DSBins contains four modules: 
#    iDSBquality, iDSBdetection, iDSBdeduplication and iDSBannotation. 
################################################################



### Usage function tells users how to run the software
helpFunction()
{
	echo "*********************************** how to use DSBins ***********************************"
	echo "Usage: sh $0 -a Sample ID -b Work Directory -c Forward Index -d Reverse Index -f forward reads -r reverse reads -p Software installed Directory"
	echo ""
	echo "DSBin consists of four parts: iDSBquality, iDSBdetection, iDSBdeduplication and iDSBannotation."
	echo -e "\t -- iDSBquality is used to eliminate the error index reads, phix reads and customed low quality reads."
	echo -e "\t -- iDSBdetection is used to detect the reads that contained large insertion events."
	echo -e "\t -- iDSBdeduplication is used to eliminate the duplicated large insertion events that caused by sequence errors, identify the representive read sequence and quality, and measure the read counts of final insertion events."
	echo -e "\t -- iDSBannotation is used to annoate the number of donors, the genetic feature of donors."
	echo ""
	echo ""
	echo -e "Request Parameters:"
	echo -e "\t-a Sample Id (Example: yYY398-B_S10)"
	echo -e "\t-b The working directory, where the raw read stored and we perform the analyses of the large insertion"
	echo -e "\t-c Index of forward reads(Example: CAG)"
	echo -e "\t-d Index of reverse reads (Example: ATC)"
	echo -e "\t-f forward reads (Example: SampleID_L001_R1_001.fastq)"
	echo -e "\t-r reverse reads (Example: SampleID_L001_R2_001.fastq)"
	echo -e "\t-p Software installed Path, Attation: This required to install DSBins, BLAST, PEAR, Bedtools in the same folder (Default:"")" 
	echo ""
	echo ""
	echo -e "Optional Parameters:"
	echo "" 
	echo -e "Optional Parameters -- overall requirements:"
	echo -e "\t-n Number of threads (Default: 15)"
	echo -e "\t-gs Genome sequence (Must corrected with chromosome ID)"
	echo -e "\t-ga Genome annotation (Default BED file format)"
	echo "" 
	echo -e "Optional Parameters of each step:"
	echo -e "iDSBquality: Optional paramter for the quality control of reads:"
	echo -e "\t-mq The minimum quality requirment of MAT region (Default: 25)"
	echo -e "\t-iq The minimum quality requirment of Inserted region (Default: 15)"
	echo ""

	echo -e "iDSBdetection: Define the MATA information:"
	echo -e "\t-il minimum length of large insertion (Default 10bp)"  
	echo -e "\t-ms Total size of whole MATA region (Default 90)"
	echo -e "\t-mc Mapped chromosme of MATA reference position (Default chrIII)"
	echo -e "\t-ms Mapped start site of MATA reference position (Default 294300)"
	echo -e "\t-me Mapped end site of MATA reference position (Default 294500)"
	echo "" 
	echo -e "iDSBdeduplication: parameters for the duplicated reads:"
	echo "" 
	echo -e "\tParameters for first cuting-edge deduplication:"
	echo -e "\t-cs cut-off of upstream and downstream, the cut-off we collected mata 11bp and inserted 19bp sequences for forward and reverse reads(default 30bp, 11+19bp)"
	echo -e "\t-cf The cutting start site of Forward reads (default 33)"
	echo -e "\t-cr The cutting start site of reverse reads (default 39)"
	echo "" 
	echo -e "\tParameters for reference-based deuplication:"
	echo -e "\t-di The minimum identity for two deduplicated reads(Default: 95)"
	echo -e "\t-do The minimum overlapping ratio for two deduplicated reads(Default: 0.95)"
	echo -e "\t-dd The maximum number of indels for two deduplicated reads, we set up the same values for different processes of deduplcations(Default: 6)"
	echo -e "\t-dm The maximum number of mismatches for two deduplicated reads, we set up the same values for different processes of deduplcations(Default: 6)"
	echo -e "\t-dl The maximum aligned gapsize of forward and reverse read that (Default: 30000)"
	echo -e "\t-dc The minimum coverage supports of each insertion events (Default:1)"
	echo -e "\t-dq The minimum quality supports of each insertion events (Default:1)"
	
	echo ""
	echo "iDSBannotation: parameters to detemine the number of donor and annotation of donors"
	echo -e "\t-hm The maxmium microhomologous length between donors (default 33)"
	echo -e "\t-ho Overlapping region between donor, less than half of min donor size (default 0.5)"
	echo ""
	
	echo -e "\t-h help"
	
	echo "For more detail information, please feel free to contact: xin.wang@childrens.harvard.edu"
	echo "**************************************************"
   exit 1 # Exit script after printing help
}


# Get absolute path for scripts and check if required scripts exist
in=$PWD

### Set default options
# some general paramaters including mata locus thread, and software path
# default threads
nproc=15
# default software path 
softwarepath=''

# Mata information
# default whole MAT length
Matasize=90
# default Mata chromosome 
Matachr="chrIII"
# default Mata start site
Matastart=294300
# default Mata end site
Mataend=294500



# default parameter for quality control
# default minimum MAT region read quality
Mqmin=25
# default minimum Inserted region read quality
Iqmin=15
# default minimum insertion length
Ilength=10

# default parameter for first cuting-edge deduplication
### the cut-off of upstream and downstream (default 30bp, 11+19bp)
Cutsize=30
CutstartF=33
CutstartR=39

# default parameter parameter for selfblast
### 
DepIden=95
DepGapsize=6
DepMismatch=6
DepCov=0.95

# default parmater for different donor identification
# The maxmium microhomologous length between donors(default 20)
MaxMicroHom=20
# the overlapping region between donor, less than half of min donor size (0.5)
OverProDonor=0.5



#default parameter for second reference-based deduplication
# The difference allowed between two reads, start site, end start site on the inserted reads (default: 6)
IDonorGap=6
# The difference allowed between two reads, start site, end start site on the chromosome locus of donor (default: 6)
CDonorGap=6


# default paramter for the estimated size of gap locus as single donor 3000
EstGap=3000


while getopts "a:b:c:d:p:f:r:gs;ga:n:ms:mc:mb:me:mq:iq:il:cs:cf:cr:di:do:dd:dm:dl:dc:dq:hm:ho:cg:ig" opt
do
   case "$opt" in
      a ) SampleID="$OPTARG" ;;
      b ) in="$OPTARG" ;;
	  c ) Findex="$OPTARG" ;;
      d ) Rindex="$OPTARG" ;;
	  p ) softpath="$OPTARG" ;;
  	  f ) Fread="$OPTARG" ;;
      r ) Rread="$OPTARG" ;;
	  
      gs) genomeseq="$OPTARG" ;;
      ga) genomeann="$OPTARG" ;;
	  
      n ) nproc="$OPTARG" ;;
      ms ) Matasize="$OPTARG" ;;
	  mc ) Matachr="$OPTARG" ;;
      mb ) Matastart="$OPTARG" ;;
	  me ) Mataend="$OPTARG" ;;
	  dl ) EstGap="$OPTARG" ;;
	  
  	  mq ) Mqmin="$OPTARG" ;;
      iq ) Iqmin="$OPTARG" ;; 	  
      il ) Ilength="$OPTARG" ;;
	  
      cs ) Cutsize="$OPTARG" ;;
	  cf ) CutstartF="$OPTARG" ;;
      cr ) CutstartR="$OPTARG" ;;
	  
	  di ) DepIden="$OPTARG" ;;
  	  do ) DepCov="$OPTARG" ;;
      dd ) CutstartR="$OPTARG" ;;
	  
      dm ) DepMismatch="$OPTARG" ;;

	  

	  dc ) FinsCount="$OPTARG" ;;
      dq ) FinsQual="$OPTARG" ;;
	  hm ) MaxMicroHom="$OPTARG" ;;
  	  ho ) OverProDonor="$OPTARG" ;;
	  
	  cg ) CDonorGap="$OPTARG" ;;
	  ig ) IDonorGap="$OPTARG" ;;

      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done


# defualt script path

srcDir=${softpath}/iDSBins/src
# default genome sequence
genomeseq=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_modified.fsa

# default genome annotation
genomeann=${softpath}/iDSBins/Database/saccharomyces_cerevisiae_R64-2-1_20150113_version6.bed

# Print helpFunction in case parameters are empty
if [ -z "${SampleID}" ] 
then
   echo "*** error: input Sample ID must be provided ***";
   helpFunction
fi

if [ -z "${in}" ] 
then
   echo "*** error: input work path must be provided ***";
   helpFunction
fi


if [ -z "${Findex}" ] 
then
   echo "*** error: input forward index must be defined ***";
   helpFunction
fi

if [ -z "${Rindex}" ] 
then
   echo "*** error: input reverse index must be defined ***";
   helpFunction
fi

if [ -z "${Fread}" ] 
then
   echo "*** error: input forward reads must be defined ***";
   helpFunction
fi

if [ -z "${Rread}" ] 
then
   echo "*** error: input reverse reads must be defined ***";
   helpFunction
fi

if [ -z "${softpath}" ] 
then
   echo "*** error: input software path must be provided ***";
   helpFunction
fi


# Begin script in case all parameters are correct
echo "${SampleID}"
echo "${in}"
echo "${Findex}"
echo "${Rindex}"
echo "${Fread}"
echo "${Rread}"
echo "${softpath}"

# 

echo "All the paramter are sucessfully provided, now let's detect the large insertion event"


#### Here is the last version of large insertion events


echo "The job array is started ..."
date


### Set up path file:

echo "Change to the Working Path, where you store your raw fastq reads"

cd ${in}

echo "Create a Sample ID folder that store all the results"
mkdir ${SampleID}
cd ${SampleID}

###############################################################################################################
#### 1. Quality Control
###############################################################################################################
echo ""
echo "iDSBquality: Quality Control is beginning ..."
date

#### 1.1 we removed trimmed low quailty sequences (less than 10)
### also removed all reads that have a 31-mer match to PhiX, allowing two mismatches. 
####detemrined the indext information

### First round the quality control
mkdir QualityControl
cd ${in}/${SampleID}/QualityControl

#fastqc ${in}/../${SampleID}_L001_R1_001.fastq   ${in}/../${SampleID}_L001_R2_001.fastq

### measure the raw reads number :
wc -l ${in}/${Fread} |awk '{ave=$1/4;print "RawReadnumber\t"ave}' >${in}/${SampleID}/${SampleID}_FinalReadStat.txt

perl  ${srcDir}/Extract_fastq_withindex.pl -f ${in}/${Fread} -r ${in}/${Rread} -i ${Findex} -g ${Rindex} -o ${SampleID}_filter

wc -l ${in}/${SampleID}/QualityControl/${SampleID}_filter.R1.fastq |awk '{ave=$1/4;print "RawRIndexReadnumber\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt


###measure the corrected index reads number:

#### 1.2 filter phix sequences

${softpath}/bbmap/bbduk.sh -in1=${SampleID}_filter.R1.fastq in2=${SampleID}_filter.R2.fastq out1=${SampleID}_filter.unmaped.R1.fastq out2=${SampleID}_filter.unmaped.R2.fastq outm1=${SampleID}_filter.mappedphix.R1.fastq outm2=${SampleID}_filter.mappedphix.R2.fastq ref=${softpath}/bbmap/resources/phix_adapters.fa.gz k=31 hdist=2 stats=stats.txt

wc -l ${in}/${SampleID}/QualityControl/${SampleID}_filter.mappedphix.R1.fastq |awk '{ave=$1/4;print "PhixReads\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt
wc -l ${in}/${SampleID}/QualityControl/${SampleID}_filter.unmaped.R1.fastq |awk '{ave=$1/4;print "RawFilterdReadnumber\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt

#### 1.3 set the stringint quality of beginning of 25 and ending of 15, here we also generate the fasta file of each type of reads 

perl ${srcDir}/Divide_high_lowQuality_v2.pl -f ${SampleID}_filter.unmaped.R1.fastq -r ${SampleID}_filter.unmaped.R2.fastq -u ${Mqmin} -g ${Iqmin} -o ${SampleID}_QC -t ${Matasize}

wc -l ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R1.fastq |awk '{ave=$1/4;print "FinalhighQualityReads\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt

### meausre the high quality reads number:
echo ""
echo "iDSBquality: Quality Control is finished ..."
date

###############################################################################################################
# 2 .Detection
###############################################################################################################
echo ""
echo "iDSBdetection: Detaction is beginning ..."
date

cd ${in}/${SampleID}/

mkdir Detection

cd ${in}/${SampleID}/Detection

#### 2.1  First round with the unmerged reads (In the aging software, we eliminated the unmerged strategy. We are going to update the method in the second verstion)

## ## We perform the blast for each forward and reverse read
# This step we firstly developed to compare the assembly and unassembly strategies.
# ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ../QualityControl/${SampleID}_QC.highquality.R1.fasta -out ${SampleID}.forward.blast.tbl  -db ${genomeseq} -num_threads 15  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1
#
# ${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ../QualityControl/${SampleID}_QC.highquality.R2.fasta -out ${SampleID}.reverse.blast.tbl -db ${genomeseq} -num_threads 15  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1
#
#
# ### Then detemine the insertion events based on the forward and reverse reads, we could set up different parameter accordingly.
# ### We generated a folder that contains the reads with large insertion events and a statistical file with the information of read quality, identity, mapped size.
# ### This statistical file will be further used to detemine the unique large insertion event.
# ### Here we used all the reads: the reads information put to the fastq file
#
# perl ${srcDir}/Determination_Insertion_nonInsertion-july1_v3.pl -i ${SampleID}.forward.blast.tbl -g ${SampleID}.reverse.blast.tbl  -f ../QualityControl/${SampleID}_QC.highquality.R1.fastq -r ../QualityControl/${SampleID}_QC.highquality.R2.fastq -q ../QualityControl/${SampleID}_QC.Allquality.stat -o ${SampleID}_detected_first -j First

#date

### 2.1 To further confirm the insertion event of each read, second round with the merge methods were applied.
## Assemble All the reads, in this case we would eliminate lots of sequence errors generated at the end of reads

# This will generate some assembled files and unassembled files
${softpath}/pear-0.9.11-linux-x86_64/bin/pear -f ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R1.fastq -r ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R2.fastq -o ${SampleID}_merged

grep ">" ${SampleID}_merged.assembled.fasta -c |awk '{print "MergedReads\t"$0}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt
wc -l ${SampleID}_merged.unassembled.forward.fastq |awk '{ave=$1/4;print "UnMergedReads\t"ave}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt


### tranfer into fastq
perl ${srcDir}/fastq2fasta.pl -i ${SampleID}_merged.assembled.fastq -o ${SampleID}_merged.assembled.fasta
perl ${srcDir}/fastq2fasta.pl -i ${SampleID}_merged.unassembled.forward.fastq -o ${SampleID}_merged.unassembled.forward.fasta
perl ${srcDir}/fastq2fasta.pl -i ${SampleID}_merged.unassembled.reverse.fastq -o ${SampleID}_merged.unassembled.reverse.fasta

 

### perform the blast analyes
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_merged.assembled.fasta -out ${SampleID}_merged.assembled.blast.tbl  -db  ${genomeseq} -num_threads ${nproc}  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_merged.unassembled.forward.fasta -out ${SampleID}_merged.unassembled.forward.blast.tbl  -db ${genomeseq} -num_threads ${nproc}  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_merged.unassembled.reverse.fasta -out ${SampleID}_merged.unassembled.reverse.blast.tbl -db ${genomeseq} -num_threads ${nproc} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' -max_target_seqs 1


cut -f 1 ${SampleID}_merged.assembled.blast.tbl |sort|uniq|wc |awk '{print "MergedMappedReads\t"$1}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt


### Measure the aligned reads:



### Detect the insertion for the second methods

# for assmbled one
perl ${srcDir}/Determination_Insertion_nonInsertion_single_blast_after_merge_v3.pl -i ${SampleID}_merged.assembled.blast.tbl -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -f ${SampleID}_merged.assembled.fasta -o ${SampleID}_detected_second_assembled -j SecondAssembly -c ${Matachr} -t ${Matastart} -e ${Mataend} -m ${Ilength} -n ${Matasize}

# for unassmbled ones
### Here we used unassembled the reads: the reads information put to the fastq file
perl ${srcDir}/Determination_Insertion_nonInsertion-july1_v3.pl -i ${SampleID}_merged.unassembled.forward.blast.tbl -g ${SampleID}_merged.unassembled.reverse.blast.tbl  -f ${SampleID}_merged.unassembled.forward.fastq -r ${SampleID}_merged.unassembled.forward.fastq -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -o  ${SampleID}_detected_second_unassembled -j SecondUnAssembly -c ${Matachr} -t ${Matastart} -e ${Mataend} -m ${Ilength} -n ${Matasize}

## 2.3 We further combined two round of large insertion events and confirm whether each read contain insertion event consistently. 

## If reads identified in both method, we confirm these reads as inserted reads. If reads only confirmed only in one method, we will check their read read and determined based on the identity (90) and quality ()

mkdir ${SampleID}_detected_final
cd ${SampleID}_detected_final

cp ../${SampleID}_detected_second_assembled/${SampleID}_detected_second_assembled.assmebled.Ainsertion.quality.txt ../${SampleID}_detected_second_unassembled/${SampleID}_detected_second_unassembled.Ainsertion.quality.txt  ./
cat ${SampleID}_detected_second_assembled.assmebled.Ainsertion.quality.txt ${SampleID}_detected_second_unassembled.Ainsertion.quality.txt >${SampleID}_detected_combined.highqual.txt
## check consistency, we will not check in the first method
#perl ${srcDir}/CheckConsistentInsertion.pl -i ${SampleID}_detected_first.Ainsertion.quality.txt -g ${SampleID}_detected_second.Ainsertion.quality.txt -o ${SampleID}_detected_combined
grep "ID" -v ${SampleID}_detected_combined.highqual.txt|wc |awk '{print "LargeInsertedReads\t"$1}' >>${in}/${SampleID}/${SampleID}_FinalReadStat.txt

echo ""
echo "iDSBdetection: First detection job array is finished ..."


#### 
echo ""
echo "iDSBdeduplication: Deduplication job array is started ..."
date
cd ${in}/${SampleID}/
mkdir DeDuplication
cd DeDuplication
mkdir DeduplicationFirst
cd ${in}/${SampleID}/DeDuplication/DeduplicationFirst

echo "First round of deduplication is started: junction based ..."

###########################################################################################################
#####  2. First deduplication:
#####.  cut the edge of two sides of each inserted reads and perform the first round of deduplication
#####.  criterion: the each junction sequence contains 19bp inserted sequence and 11bp MATA sequences
###########################################################################################################



cp ${in}/${SampleID}/Detection/${SampleID}_detected_final/${SampleID}_detected_combined.highqual.txt ./


### To reduce the speed of the following analyses, we firstly identify the complete identical fasta id:
perl ${srcDir}/Identical_ReadCaculation_v3.pl -g ${in}/${SampleID}/Detection/${SampleID}_merged.assembled.fastq -f ${in}/${SampleID}/Detection/${SampleID}_merged.unassembled.forward.fastq -r ${in}/${SampleID}/Detection/${SampleID}_merged.unassembled.reverse.fastq -i ${SampleID}_detected_combined.highqual.txt -o ${SampleID}_detected_combined.highqual2.txt


### 1.1  we combined the two junction sequences of all detected reads
perl ${srcDir}/Cut_twoedge_sequences_forunassembled_v2.pl  -f ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R1.fastq -r ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R2.fastq -i ${SampleID}_detected_combined.highqual2.txt -o ${SampleID}_detected_cut60.fasta -u ${Cutsize} -t ${CutstartF} -e ${CutstartR}


#### 1.2 Blast against itself
${softpath}/ncbi-blast-2.8.1+/bin/makeblastdb -dbtype nucl -in ${SampleID}_detected_cut60.fasta
 
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_detected_cut60.fasta -db ${SampleID}_detected_cut60.fasta -out ${SampleID}_detected_cut60.blast -num_threads ${nproc} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

### 1.3 We perform the first round of deduplications using the junction sequences, attention we also set up the assembly reads as the first priority if the unassmebled reads and assembeld reads at the same cluster.
perl ${srcDir}/SelfBlast_RemoveDuplication_eachsampleIdentitySorted_v2.pl  -b ${SampleID}_detected_cut60.blast -i ${SampleID}_detected_combined.highqual2.txt -f ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R1.fastq -r ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.R2.fastq -o ${SampleID}_detected_FirstDep -t ${DepIden} -g ${DepGapsize} -m ${DepMismatch} -c ${DepCov}
 

echo "First round of deduplication is ended: junction based ..."

###########################################################################################################
#####  2. Second deduplication:
#####.  we used the referenced based method to perform the deduplication, allowing the frameshift 
#####.  criterion: 
###########################################################################################################

echo "First round of deduplication is start: overall sequences based ..."

cd ${in}/${SampleID}/DeDuplication
mkdir DeduplicationSecond
cd DeduplicationSecond

cp ${in}/${SampleID}/DeDuplication/DeduplicationFirst/${SampleID}_detected_FirstDep.cls ./

### 2.1 Prepair the required information of second round deduplication: the unassembled reads, assembled reads, and their blast alignments. 


###    Extract the representive reads detected from the first round, then seperate them into assembled and unassembled reads.

perl ${srcDir}/SperatedInsertionbasedAssembly.pl -g ${in}/${SampleID}/Detection/${SampleID}_merged.assembled.fastq -f ${in}/${SampleID}/Detection/${SampleID}_merged.unassembled.forward.fastq -r ${in}/${SampleID}/Detection/${SampleID}_merged.unassembled.reverse.fastq -i ${SampleID}_detected_FirstDep.cls -o ${SampleID}_detected_assemblysep



#### 2.1.1 for the assmbled reads, we also perform the deduplicaton using the selfblast. The reason why we used here is because some MAT region due have lots of mismatches, sequence error that result in the frameshift of the junction based. And another reason we could choose the representive read for each cluster mainly based on the read counts ranking. Here we sorted the quality and identity firstly and final the read count for each clustering inserted reads.

### We used selfblast to perform the consequent deduplication
${softpath}/ncbi-blast-2.8.1+/bin/makeblastdb  -in ${SampleID}_detected_assemblysep.Insassembled.fasta -dbtype nucl
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${SampleID}_detected_assemblysep.Insassembled.fasta  -db ${SampleID}_detected_assemblysep.Insassembled.fasta  -out ${SampleID}_detected_assemblysep.Insassembled.self.blast -num_threads ${nproc}  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

### Here we could update the number of identity, gap size, overlapping rate of two reads (we set the default identity >95, overlapping rate >95, mismatches <=6, gapsize <=6)
perl ${srcDir}/SelfBlast_RemoveDuplication_eachsampleIdentitySorted_totalinsertion_v2.pl -b ${SampleID}_detected_assemblysep.Insassembled.self.blast -i ${SampleID}_detected_assemblysep.CoverageWithAssembly.txt -f ${SampleID}_detected_assemblysep.Insassembled.fasta -o ${SampleID}_detected_assemblysep_firstdedup -t ${DepIden} -g ${DepGapsize} -m ${DepMismatch} -c ${DepCov}
 
#### 2.1.3 For the unassembled reads, we cancatenated the forward and reverse reads and then perform the selfblast. Similar strategy were used for these reads
### We firstly concantenated the forward and reverse reads
perl ${srcDir}/ConcatenateForwardReverseReads.pl -f ${SampleID}_detected_assemblysep.InsUnassembled.R1.fastq -r ${SampleID}_detected_assemblysep.InsUnassembled.R2.fastq -o ${SampleID}_detected_assemblysep.InsUnassembled
${softpath}/ncbi-blast-2.8.1+/bin/makeblastdb  -in ${SampleID}_detected_assemblysep.InsUnassembled.Concatenated.fasta -dbtype nucl
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_detected_assemblysep.InsUnassembled.Concatenated.fasta -db ${SampleID}_detected_assemblysep.InsUnassembled.Concatenated.fasta -out ${SampleID}_detected_assemblysep.InsUnassembled.Concatenated.self.blast -num_threads ${nproc}  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'
perl  ${srcDir}/SelfBlast_RemoveDuplication_eachsampleIdentitySorted_Concanatedinsertion_v2.pl -b ${SampleID}_detected_assemblysep.InsUnassembled.Concatenated.self.blast -f ${SampleID}_detected_assemblysep.InsUnassembled.R1.fastq -r ${SampleID}_detected_assemblysep.InsUnassembled.R2.fastq -i ${SampleID}_detected_assemblysep.CoverageWithAssembly.txt -o ${SampleID}_detected_assemblysep.InsUnassembled.First


echo "First round of deduplication is ended: overall sequences based ..."

### 2.2 Divide the assembled file into different categories based on the number of insertions

### For assmebled we determine the number of large insertions

## double check the blastn results

echo "Second round of deduplication is ended: reference based ..."

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_detected_assemblysep_firstdedup.FAssuniq.fasta -out ${SampleID}_detected_assemblysep.Insassembled.blast_all_detected.tbl -db ${genomeseq} -num_threads ${nproc}  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

perl ${srcDir}/Determine_insertion_number_formutipleinserts_forallassembled_universal_updated_v3.pl -g ${SampleID}_detected_assemblysep.Insassembled.blast_all_detected.tbl -i ${SampleID}_detected_assemblysep_firstdedup.FAssuniq.fasta  -c ${SampleID}_detected_assemblysep_firstdedup.cls -o ${SampleID}_detected_number -t ${MaxMicroHom} -e ${OverProDonor} -m ${Ilength}

### Due to many insertions cannot identify blast results, we set a BLAST method using blastn-short 
## Reference : https://www.metagenomics.wiki/tools/blast/default-word-size
## Decreasing the word-size will increase the number of detected homologous sequences, we set up 11 as word size

${softpath}/ncbi-blast-2.8.1+/bin/blastn  -task blastn-short  -word_size 11 -evalue 0.001 -dust no -soft_masking false -gapopen 5 -penalty -1 -perc_identity 80 -query  ${SampleID}_detected_number.noaligned.fasta -out ${SampleID}_detected_number.noaligned.blast.tbl -db ${genomeseq} -num_threads ${nproc}  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

perl ${srcDir}/Determine_insertion_number_formutipleinserts_forallassembled_universal_updated_v3.pl  -g ${SampleID}_detected_number.noaligned.blast.tbl -i ${SampleID}_detected_number.noaligned.fasta -c ${SampleID}_detected_number.noaligned.txt -o ${SampleID}_detected_number.noaligned  -t ${MaxMicroHom} -e ${OverProDonor} -m ${Ilength} 


### We integrated all these insertion with one donor from assembly and unassembly


### The following should add the Description of each scripts


#### 1 cat single donor insertion, here included some mutiple inseritons but only detected one donor, single donor and single estimated donor### we allow for 6bp frameshife of 

cat ${SampleID}_detected_number.single.txt ${SampleID}_detected_number.noaligned.single.txt >${SampleID}_detected_second_singleassembly.txt
perl ${srcDir}/Single_confident_donor_coverage_deduplication_v4.pl -i ${SampleID}_detected_second_singleassembly.txt -a ${IDonorGap} -b ${CDonorGap}  -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.stat -o ${SampleID}_detected_second_singleassembly 


### 2. two donor insertion 

cat ${SampleID}_detected_number.twoinf.txt ${SampleID}_detected_number.noaligned.twoinf.txt >${SampleID}_detected_second_twoassembly.txt

perl ${srcDir}/Two_confident_donor_coverage_deduplication_universal_v4.pl -i ${SampleID}_detected_second_twoassembly.txt -a ${IDonorGap} -b ${CDonorGap} -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -o ${SampleID}_detected_second_twoassembly


### 3. Three donor insertion
cat ${SampleID}_detected_number.threeinf.txt ${SampleID}_detected_number.noaligned.threeinf.txt >${SampleID}_detected_second_threeassembly.txt
perl ${srcDir}/Three_confident_donor_coverage_deduplication_universal_v4.pl -i ${SampleID}_detected_second_threeassembly.txt -a ${IDonorGap} -b ${CDonorGap} -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -o ${SampleID}_detected_second_threeassembly


### 4. Four donor insertion

cat ${SampleID}_detected_number.fourinf.txt ${SampleID}_detected_number.noaligned.fourinf.txt >${SampleID}_detected_second_fourassembly.txt

perl ${srcDir}/Four_confident_donor_coverage_deduplication_universal_v4.pl -i ${SampleID}_detected_second_fourassembly.txt -a ${IDonorGap} -b ${CDonorGap} -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -o ${SampleID}_detected_second_fourassembly


## 5. Final Undetected insertions, we will double check whether they came from MAT

## here is also including the single insertion events from MAT

mkdir MAT
cd MAT
perl ${srcDir}/Seperate_unmappable_mat_uniq_v3.pl -f  ../${SampleID}_detected_number.noaligned.noaligned.fasta -b ../../../Detection/${SampleID}_merged.assembled.blast.tbl -o ${SampleID}_assembled_final_undetermined
### for the MAT deduplication:  ### we have a high standard for the insertion events
perl ${srcDir}/Single_confident_donor_coverage_deduplication_v4.pl -i ${SampleID}_assembled_final_undetermined.mat.txt -a 1 -b 1  -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.highquality.stat -o ${SampleID}_final_undetermined.mat



### For the Unmapped insertion , we require a high standard (Cov >=2 & Quality >=25)
perl -ne '{chomp; if (/>(\S+)\s+(\S+)\s+(\S+)/){$id =$1; $cov{$id}=$2; $qual{$id}=$3}else{$str{$id}.=$_;}}END {foreach my $i (keys %cov){if ($cov{$i}>1 && $qual{$i} >=25){print ">$i\t$cov{$i}\t$qual{$i}\n$str{$i}\n"}}}' ${SampleID}_assembled_final_undetermined.unmappable.fasta >${SampleID}_final_undetermined.fasta



##### 6 single estimated donor insertion

#### for the unassmbled sequences #### 
### for the unassembled sequence, we mainly estimated the large insertion if two reads blast against two side of simialr region within 3kb

### Double check the unassembled sequence and generated the confident blast file of large insertions that can identity the donor information in both reads.
### Here we used unassembled the reads: the reads information put to the fastq file
cd ../
perl ${srcDir}/Determination_Insertion_nonInsertion-july1_v3.pl -i ${in}/${SampleID}/Detection/${SampleID}_merged.unassembled.forward.blast.tbl -g ${in}/${SampleID}/Detection/${SampleID}_merged.unassembled.reverse.blast.tbl  -f ${SampleID}_detected_assemblysep.InsUnassembled.First.FUassuniq.forward.fastq -r ${SampleID}_detected_assemblysep.InsUnassembled.First.FUassuniq.reverse.fastq -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -o ${SampleID}_unassembledEstimate -j UnassEst  -c ${Matachr} -t ${Matastart} -e ${Mataend} -m ${Ilength} -n ${Matasize}

cd ${SampleID}_unassembledEstimate

  ##### estimate the single insertion  #### here we can only determine the single insertion
perl ${srcDir}/Yeast_overalapping_insertion_from_R1R2_v3.pl -i ${SampleID}_unassembledEstimate.insert.blast.txt -c ../${SampleID}_detected_assemblysep.CoverageWithAssembly.txt -g ${EstGap} -o ${SampleID}_unassembledEstimate

# Here we will not consider the insertions that cannot well aligned to the end of R1 and the begining of R2, we allow 5bp shift at each end 
perl ${srcDir}/Single_confident_donor_coverage_deduplication_estimated_v4.pl -i ${SampleID}_unassembledEstimate.oneinsert.txt -a 6 -b 6 -s 5 -q ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -o ${SampleID}_detected_second_estimated

echo "Referenced deduplication job array is finshed ..."

echo ""
echo "iDSBdeduplication: deduplication job array is finshed ..."
date

##############################################################################
### Part3.  extract the sequence and export the number of donors and annotation####
##################################################################################
echo ""
echo "iDSBannotation: annotation job array is starting ..."
echo ""
echo "Perform annotation of each categories of large insertions"

## 3.1 Prepair all the poteintial feature information of each assembled reads
# Identify all blast of each assembled reads with default paramter:


### For assembled: Blast aginst all assembled (First step using the default blast)
cd ../

#${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_detected_assemblysep.Insassembled.fasta -out ${SampleID}_detected_assemblysep.Insassembled.blast_all_detected.tbl -db ${genomeseq} -num_threads 15  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue' 

### integrate the default blast and the unaligned blast (shorter parameter):
cat ${SampleID}_detected_number.noaligned.blast.tbl ${SampleID}_detected_assemblysep.Insassembled.blast_all_detected.tbl >${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.tbl

### Generate the bedfile index:
perl ${srcDir}/Create_bedfile_forblast_determinemutiplelocus_v3.pl -i ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.tbl -o ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.bed

## sort and annotation
sort -k 1V,1 -k 2n,2 ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.bed >${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.sorted.bed
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ${genomeann} -a ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.sorted.bed >${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed

### Annotate the insertion with single donor, we need the blast result to confirm whether they have the alternative feature
cat ${SampleID}_detected_second_singleassembly.uniq.txt MAT/${SampleID}_final_undetermined.mat.uniq.txt >${SampleID}_detected_secondcombined_singleassembly.uniq.txt
perl ${srcDir}/Annotation_withupdownstream_insertion_For1donor_v3.pl -a ${SampleID}_detected_assemblysep.Insassembled.fasta -b ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -c ${SampleID}_detected_secondcombined_singleassembly.uniq.txt -d ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.tbl -e ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed -t ${SampleID} -o ${SampleID}_final_single


### Annotate the insertion for two donor
perl ${srcDir}/Annotation_withupdownstream_insertion_For2donors_v3.pl -a ${SampleID}_detected_assemblysep.Insassembled.fasta -b ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -c ${SampleID}_detected_second_twoassembly.uniq.txt -d ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.tbl -e ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed -t ${SampleID} -o ${SampleID}_final_two

### Annotate the insertion for three donor
perl ${srcDir}/Annotation_withupdownstream_insertion_For3donors_v3.pl -a ${SampleID}_detected_assemblysep.Insassembled.fasta -b ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -c ${SampleID}_detected_second_threeassembly.uniq.txt -d ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.tbl -e ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed -t ${SampleID} -o ${SampleID}_final_three

### Annotate the insertion for four donor

perl ${srcDir}/Annotation_withupdownstream_insertion_For4donors_v3.pl -a ${SampleID}_detected_assemblysep.Insassembled.fasta -b ${in}/${SampleID}/QualityControl/${SampleID}_QC.Allquality.stat -c ${SampleID}_detected_second_fourassembly.uniq.txt -e ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.annotation.bed -d ${SampleID}_detected_assemblysep.Insassembled.FinallyBlast.tbl -t ${SampleID} -o ${SampleID}_final_four

### Annotate the estimated insertion for esimated donor

cd ${SampleID}_unassembledEstimate

perl ${srcDir}/Change_ID_fordeduplicatedestmateone_insertion_tobedsequence_versionOct_v3.pl -g ${genomeseq} -m ${SampleID}_unassembledEstimate.Ainsertion.R1.fastq -n ${SampleID}_unassembledEstimate.Ainsertion.R2.fastq -f ${SampleID}_detected_second_estimated.uniq.txt -t ${SampleID} -o ${SampleID}_estimated

### double check whether they do have alternative feature

${softpath}/ncbi-blast-2.8.1+/bin/blastn -query ${SampleID}_estimated.wholeseq.fasta -out ${SampleID}_estimated.wholeseq.blast -db ${genomeseq} -num_threads 15  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'


### Annotate the estimated insertion.
perl -alne '{my ($start,$end)=($F[10]>$F[9])?($F[9],$F[10]):($F[10],$F[9]); my $strand= ($F[10]>$F[9])?"+":"-"; my $string= join "\t", $F[1],$start,$end; next if (exists $hash{$string}) ; $hash{$string}++; print "$string\t$F[0]\t0\t$strand"}' ${SampleID}_estimated.wholeseq.blast |sort -k 1,1V -k 2,2n >${SampleID}_estimated.wholeseq.bed
${softpath}/bedtools2/bin/bedtools closest -D b -t first -b ${genomeann} -a ${SampleID}_estimated.wholeseq.bed >${SampleID}_estimated.wholeseq.annotation.bed

### Add the alterative feature for estimated insertion
perl ${srcDir}/FurtherEstimatedFeature_annotations_v3.pl -i ${SampleID}_estimated.txt -b ${SampleID}_estimated.wholeseq.blast -c ${SampleID}_estimated.wholeseq.annotation.bed -o ${SampleID}_estimated.final.txt


echo "iDSBannotation: annotation job array is finished ..."
date
############################################
### Generate final insertion events
############################################


cd ../../
mkdir DeduplicationFinal
cd DeduplicationFinal



################################################################################################################
### Step 5 Further Blast to check the low coverage existence (Self Blast) it would be better to put here before 
### concatenate them together
################################################################################################################

echo "Generating final results with high standards ..."
cp ../DeduplicationSecond/${SampleID}_final_single.Ftable ../DeduplicationSecond/${SampleID}_final_two.Ftable ../DeduplicationSecond/${SampleID}_final_three.Ftable ../DeduplicationSecond/${SampleID}_final_four.Ftable ../DeduplicationSecond/MAT/${SampleID}_assembled_final_undetermined.unmappable.fasta  ../DeduplicationSecond/${SampleID}_unassembledEstimate/${SampleID}_estimated.final.txt ./

############################################
# Final round of self blast 
############################################
cat ../DeduplicationSecond/${SampleID}_final_single.Ftable ../DeduplicationSecond/${SampleID}_final_two.Ftable ../DeduplicationSecond/${SampleID}_final_three.Ftable ../DeduplicationSecond/${SampleID}_final_four.Ftable ../DeduplicationSecond/${SampleID}_unassembledEstimate/${SampleID}_estimated.final.txt >${SampleID}_mappable.txt
perl -ne '{chomp; my ($id,$string,$identity,$Rcount,$quality)=(split/\t/,$_)[1,3,7,8,9]; next if ($id eq "NID" || exists $hash{$id}); $hash{$id}++; print ">$id\t$Rcount\t$quality\t$identity\n$string\n" }' ${SampleID}_mappable.txt >${SampleID}_mappable.fasta

cp ../DeduplicationSecond/MAT/${SampleID}_final_undetermined.fasta ${SampleID}_unmappable_raw.fasta
perl -ne '{chomp; if ($_=~/>/){print "$_\t0\n"}else{print "$_\n"}}' ${SampleID}_unmappable_raw.fasta >${SampleID}_unmappable.fasta

#
cat ${SampleID}_mappable.fasta ${SampleID}_unmappable.fasta >${SampleID}_all.fasta

${softpath}/ncbi-blast-2.8.1+/bin/makeblastdb  -in ${SampleID}_all.fasta -dbtype nucl
${softpath}/ncbi-blast-2.8.1+/bin/blastn -query  ${SampleID}_all.fasta -db ${SampleID}_all.fasta -out ${SampleID}_all.blast -num_threads 15  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen bitscore evalue'

#### For the final round deduplication, we briefly adapt self-blast strategy to remove these duplicates. Considering the fact that the longer the reads are, the lower quality the reads will be. MiSeq read length is 300bp, error rate can be as high as 20%. Therefore, we specifically required a loose standard for long insertion deduplicates, the other short still used the same parameters. Regarding the long insertion events, we removed long insertions low-coverage duplicates that are with similar length (difference <=5), high identity (>80%), less mismatches/indels (<= 15% insertion length). 

perl ${srcDir}/SelfBlast_RemoveDuplication_eachsampleIdentitySorted_totalinsertion_finalround_v2.pl -b ${SampleID}_all.blast  -i ${SampleID}_mappable.txt -f ${SampleID}_all.fasta -o ${SampleID}_final


# perl ${srcDir}/SelfBlast_RemoveDuplication_eachsampleIdentitySorted_totalinsertion_finalround_v2.pl -b ${SampleID}_all.blast  -i ${SampleID}_mappable.txt -f ${SampleID}_all.fasta -o ${SampleID}_final
#
# #### Here we ignored the mismatches as long as the sequence identity
#
# #### Due to the mismatches inside the sequence, here we ignored the mismathes for low covereage insertions () but required sequence identity more than 85, shift less than bp
# perl -ne '{chomp; my @array=split; next if ($array[0] eq $array[1]); my $av=$array[3]/$array[8];my $dif=abs($array[8]-$array[11]); if ($av>=0.85 && $array[2]>=85 && $dif<=2) {print "$_\n"}}'  ${SampleID}_all.blast >${SampleID}_all.potential.txt
#
# perl /temp_work/ch220812/Project/script/DSBinsertion/Version_Jan27/SelfBlast_RemoveDuplication_eachsample.pl -i ${SampleID}_all.potential.txt -g ${SampleID}_all.fasta -o ${SampleID}_all.afterdup.txt
#

##### Updated version #####  srcDir=/lab-share/Cardio-Chen-e2/Public/xwang/tmhxxw9/tmhxxw9/project/Aging/Updated_V0903/UsedMutants/TestSoftWare/Software/iDSBins/src

############################################
# Generate final round of table
############################################
perl ${srcDir}/FinalInsertionQualityControl.pl  -g ${SampleID}_final.cls -i ${SampleID}_final.finalinsertion.txt -o ${SampleID}_final.finalinsertion.highquality
 
 awk 'NR==1; NR > 1 {print $0 | "sort  -k 11,11V -k 12,12n -k 13,13n"}' ${SampleID}_final.finalinsertion.highquality.One.txt|perl -ne '{chomp; my ($chr,$start,$end)=(split/\t/,$_)[10,11,12]; my $up; my $down;if (exists $hash{$chr}){$up=$start-$start0; $down=$end-$end0}else{$up=$start; $down=$end} print "$_\t$up\t$down\n"; $start0=$start; $end0=$end; $hash{$chr}++;}' >${SampleID}_final_One_updated.txt
 awk 'NR==1; NR > 1 {print $0 | "sort  -k 11,11V -k 12,12n -k 13,13n"}' ${SampleID}_final.finalinsertion.highquality.Multiple.txt  |perl -ne '{chomp; my ($chr,$start,$end)=(split/\t/,$_)[10,11,12];if ($start eq "Unknown" || $start  eq "NO"){print "$_\tUnknown\tUnknown\n"; next}; my $up; my $down;if (exists $hash{$chr}){$up=$start-$start0; $down=$end-$end0}else{$up=$start; $down=$end} print "$_\t$up\t$down\n"; $start0=$start; $end0=$end; $hash{$chr}++;}'  >${SampleID}_final_Multiple_updated.txt
 
 #### Then genetate the final result by removin the duplicates that have 0-10 bp shift and low coverage, and add the new read count number
 mkdir -p ${in}/${SampleID}/FinalInsertion
 ## for sinlge donor: 
 perl ${srcDir}/CorrectedSingleBasedOnUpdistDowndist.pl -i ${SampleID}_final_One_updated.txt -o ${SampleID}_final_One_updated
 awk 'NR==1; NR > 1 {print $0 | "sort  -k 11,11V -k 12,12n"}' ${SampleID}_final_One_updated.OneCorrected.txt | perl -ne '{chomp; $n++; my @array=split/\t/,$_; if($n==1){print "$_\n"}else{$m=$n-1;$array[0]="S$m";my $string=join "\t",@array;print "$string\n"}}' >${in}/${SampleID}/FinalInsertion/${SampleID}_final.sinlge.txt
 
 ## for multiple donor
 
 perl ${srcDir}/FinalMutipleInsertionSortUpdate.pl -i ${SampleID}_final_Multiple_updated.txt -o ${SampleID}_final_Multiple_updated.MuCorrected.txt
 
 ## Change ID
  perl -ne '{chomp; my @array=split/\t/,$_; if ($array[0] eq "CaseID"){print "$_\n"; next;}; my $string; if ($array[7] eq "Unknown"){$array[6] = "1orMore";} if (!exists $hash{$array[1]}){$n++;$array[0]="M$n"; $string=join "\t",@array; }else{my ($i,$m)=split/\./,$array[0]; $array[0]="M$n.$m"; $string=join "\t",@array}; print "$string\n"; $hash{$array[1]}++; }' ${SampleID}_final_Multiple_updated.MuCorrected.txt >${SampleID}_final.multiple.txt
 
 
 cp ${SampleID}_final.multiple.txt ${in}/${SampleID}/FinalInsertion/${SampleID}_final.multiple.txt

 ############################################
 # Double check the most abundance reads 
 #  and our final insertion events. 
 #      Here we can track reads ID.
 ############################################
 
cd ${in}/${SampleID}/FinalInsertion/

### you can use head 20 to track the ID print the top 20 
sort -k 7nr ${in}/${SampleID}/DeDuplication/DeduplicationFirst/${SampleID}_detected_combined.highqual2.txt |head -n 20 
 
 
date

echo "Congratulation! Insertion detection and deduplication job array is finished !"



