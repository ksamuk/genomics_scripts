#!/bin/bash

#USAGE: ./SB_gatk_pipe.sh <plate identifier> 
#Sept 2013 GLO
#Pipeline for using GATK and BWA to call snps
plate="$1"
Gpath='/home/samuk/bin'
Hpath='/home/samuk/gbs2015'
rawdata='/home/samuk/gbs2015/raw'
bwa='/home/samuk/bin/bwa/bwa-0.7.10'
stampy='/home/samuk/bin/stampy-1.0.23/stampy.py'
#demultiplex='/home/owens/bin/GBS_fastq_Demultiplexer_v8.GO.pl'
#barcodes="/home/owens/SB/Barcodes.$plate.txt"
ref='/home/samuk/review/ref/GA.broad1_75_dna_rm.fa'
stampyref='/home/samuk/review/ref/GA.broad1_75_dna_rm.fa.stampyref'
cleandata="/home/samuk/gbs2015/raw/combined"
trimmeddata="/home/samuk/gbs2015/raw/$plate.trimmed_data_paired"
unpaired="/home/ssamuk/gbs2015/raw/$plate.trimmed_data_unpaired"
sam="/home/samuk/gbs2015/raw/$plate.sam"
bam="/home/samuk/gbs2015/raw/$plate.bam"
log="/home/samuk/gbs2015/raw/$plate.log"
gvcf="/home/samuk/gbs2015/raw/gvcf"
tabix='/home/samuk/bin/tabix-0.2.6'
vcftools='/home/samuk/bin/vcftools_0.1.11/bin'
picardtools='/home/samuk/bin/picard-tools-1.96'
trim="/home/samuk/bin/Trimmomatic-0.32"
project='whtstbk_gbs_2015_L1'
#Prerequisites:
#Index the reference for GATK and BWA

#/home/samuk/bin/stampy-1.0.23/stampy.py -G GA.broad1_75_dna_rm.fa.stampyref ./GA.broad1_75_dna_rm.fa
#/home/samuk/bin/stampy-1.0.23/stampy.py -g GA.broad1_75_dna_rm.fa.stampyref -H GA.broad1_75_dna_rm.fa.stampyref

#Change the variables to fit your file structure and ensure the folders exist.

#set java opts
_JAVA_OPTIONS="-Xmx18g"

echo "Starting GBS processing pipeline on $plate."
#Make directories if they don't exist
if [ ! -d "$sam" ]; then
	mkdir $sam
fi
if [ ! -d "$bam" ]; then
	mkdir $bam
fi
if [ ! -d "$log" ]; then
	mkdir $log
fi
if [ ! -d "$gvcf" ]; then
        mkdir $gvcf
fi
if [ ! -d "$cleandata" ]; then
        mkdir $cleandata
fi
if [ ! -d "$trimmeddata" ]; then
        mkdir $trimmeddata
fi
if [ ! -d "$unpaired" ]; then
        mkdir $unpaired
fi

#Demultiplex
perl $demultiplex $barcodes $rawdata/"$plate"_R1.fastq $rawdata/"$plate"_R2.fastq $cleandata/$project


ls $cleandata | sed s/.fq// | uniq  > $Hpath/Samplelist.$plate.txt
# Trim the data using Trimmomatic. Removes bad reads and illumina adapter contamination.
while read prefix
do
java -jar $trim/trimmomatic-0.32.jar SE -phred33 $cleandata/"$prefix".fq $trimmeddata/"$prefix".fastq ILLUMINACLIP:$trim/adapters/TruSeq3-SE.fa:2:30:10:8:T TRAILING:3 SLIDINGWINDOW:4:15
done < Samplelist.$plate.txt


##Align using BWA. Turn from sam to bam. Sort by coordinate and add read group data.
while read prefix
do
	echo "Aligning $prefix"
        $bwa/bwa aln -t 8 $ref $trimmeddata/"$prefix"_R1.fastq 1> $trimmeddata/"$prefix"_R1.sai
        $bwa/bwa aln -t 8 $ref $trimmeddata/"$prefix"_R2.fastq 1> $trimmeddata/"$prefix"_R2.sai
        $bwa/bwa sampe $ref $trimmeddata/"$prefix"_R1.sai $trimmeddata/"$prefix"_R2.sai $trimmeddata/"$prefix"_R1.fastq $trimmeddata/"$prefix"_R2.fastq 1> $sam/$prefix.sam 2> $log/$prefix.bwasampe.log
        samtools view -Sb $sam/$prefix.sam > $bam/$prefix.bam
        $stampy -g $stampyref -h $stampyref -t8 --bamkeepgoodreads -M $bam/$prefix.bam -o $bam/$prefix.stampy.bam 2> $log/$prefix.stampy.log
        java -jar $picardtools/CleanSam.jar INPUT=$bam/$prefix.stampy.bam OUTPUT=$bam/$prefix.clean.bam 2> $log/$prefix.cleansam.log
        java -jar $picardtools/SortSam.jar INPUT=$bam/$prefix.clean.bam OUTPUT=$bam/$prefix.sort.bam SORT_ORDER=coordinate 2> $log/$prefix.sortsam.log
        java -jar $picardtools/AddOrReplaceReadGroups.jar I=$bam/$prefix.sort.bam O= $bam/$prefix.sortrg.bam SORT_ORDER=coordinate RGID=$prefix RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$prefix CREATE_INDEX=True 2> $log/$prefix.addRG.log
        rm $trimmeddata/"$prefix"_R1.sai
        rm $trimmeddata/"$prefix"_R2.sai
        rm $sam/$prefix.sam
        rm $bam/$prefix.bam
        rm $bam/$prefix.stampy.bam
        rm $bam/$prefix.clean.bam
        rm $bam/$prefix.sort.bam
done < Samplelist.$plate.txt

exit

#Make bam.list for GATK
ls -d $bam/*.* | grep sortrg.bam  > $Hpath/bamlist.$plate.list

#identify local indels

java -Xmx2g -jar $Gpath/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R $ref \
   -I $Hpath/bamlist.$plate.list \
   -nt 8 \
   -log $log/$plate.RealignerTargetCreator.log \
   -o $Hpath/$plate.realign.intervals

Realign around local indels
while read prefix
do
java -Xmx12g -d64 -jar $Gpath/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $ref \
	-I $bam/$prefix.sortrg.bam \
	-targetIntervals $Hpath/$plate.realign.intervals \
	-o $bam/$prefix.realign.bam \
	-log $log/$plate.$prefix.IndelRealigner.log 
done < $Hpath/Samplelist.$plate.txt

while read prefix
do

exit

#Call GATK HaplotypeCaller
java -Xmx18g -jar $Gpath/GenomeAnalysisTK.jar \
	-nct 4 \
	-l INFO \
	-R $ref \
	-log $log/$plate.$prefix.HaplotypeCaller.log \
	-T HaplotypeCaller \
	-I  $bam/$prefix.realign.bam \
	--emitRefConfidence GVCF \
	-variant_index_type LINEAR \
	-variant_index_parameter 128000 \
	-o $gvcf/$prefix.GATK.gvcf.vcf
done < $Hpath/Samplelist.$plate.txt

#Make input list for GATK GenotypeGVCFs
tmp=""
while read prefix
do
        tmp="$tmp --variant $gvcf/$prefix.GATK.gvcf.vcf"
done < $Hpath/Samplelist.$plate.txt

#Genotype all gvcf together into one vcf file
java -Xmx18g -jar $Gpath/GenomeAnalysisTK.jar \
	-nt 8 \
	-l INFO \
	-R $ref \
	-log $log/$plate.$prefix.GenotypeGVCFs.log \
	-T GenotypeGVCFs \
	$tmp \
	-o $Hpath/SB.GATK.total.vcf \
	-inv


exit

