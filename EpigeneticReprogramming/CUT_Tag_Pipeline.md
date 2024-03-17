## Trim Reads

```
adapter=adapter.fa
trimmomatic PE -threads 8 -phred33 *_R1.fastq.gz *_R2.fastq.gz ${PWD##*/}_clean_1.fastq.gz output_forward_unpaired.fq.gz ${PWD##*/}_clean_2.fastq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:$adapter:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36
```

## Map Reads

```
bowtie2 --local --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 2000 -x /root/vol_1T/Reference/bowtie2_Ms/mm10 -1 *clean_1.fastq.gz  -2 *clean_2.fastq.gz > ${PWD##*/}.sam

samtools view -@ 2 -Sb -q 10 ${PWD##*/}.sam | samtools sort -@ 2 - > ${PWD##*/}_sort.bam

samtools index ${PWD##*/}_sort.bam
```

## Remove Duplicates and Reads in Blacklist

```
java -jar ~/picard.jar MarkDuplicates -I ${PWD##*/}_sort.bam -O ${PWD##*/}_rmDup.bam -REMOVE_DUPLICATES true -M sample_RmDup_metrics.txt
samtools index ${PWD##*/}_rmDup.bam

bedtools intersect -v -abam ${PWD##*/}_rmDup.bam -b /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed > ${PWD##*/}_clean.bam
samtools index ${PWD##*/}_clean.bam
```

## Peak Calling Using MACS

```
macs3 callpeak -t ${PWD##*/}_clean.bam -f BAMPE --broad -g mm --keep-dup all --broad-cutoff 0.05 -n ${PWD##*/}
macs3 callpeak -t ${PWD##*/}_clean.bam -f BAMPE -g mm --keep-dup all -q 0.25 -n ${PWD##*/}
```

## BAM to Bigwig files

```
bamCoverage -b ${PWD##*/}_clean.bam --numberOfProcessors 2 --effectiveGenomeSize 2652783500 --normalizeUsing RPKM --extendReads 200 --binSize 50 --outFileFormat bigwig --outFileName ${PWD##*/}.bw
```

## Profile Plot

```
KC_GTF=KC_2019.gtf
SAMac_GTF=SAMac_2019.gtf

## PU.1
PU1_Met=/root/vol_1T/Fig3/PU1_Met/PU1_Met.bw
PU1_N=/root/vol_1T/Fig3/PU1_N/PU1_N.bw

computeMatrix scale-regions -S $PU1_N $PU1_Met -R $KC_GTF $SAMac_GTF --samplesLabel PU1_N PU1_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName PU1_KC_SAMac.mat.gz
plotProfile --matrixFile PU1_KC_SAMac.mat.gz --outFileName PU1_KC_SAMac.profile.pdf --perGroup --colors '#6a3d9a' '#ff7f00' --plotType se --plotFileFormat pdf

computeMatrix scale-regions -S $PU1_N $PU1_Met -R $KC_GTF --samplesLabel PU1_N PU1_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName PU1_KC.mat.gz
computeMatrix scale-regions -S $PU1_N $PU1_Met -R $SAMac_GTF --samplesLabel PU1_N PU1_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName PU1_SAMac.mat.gz
plotHeatmap --matrixFile PU1_KC.mat.gz --outFileName PU1_KC.heatmap.pdf --plotType se --plotFileFormat pdf
plotHeatmap --matrixFile PU1_SAMac.mat.gz --outFileName PU1_SAMac.heatmap.pdf --plotType se --plotFileFormat pdf

## H3K27me3
H3K27me3_Met=/root/vol_1T/Fig3/H3K27me3_Met/H3K27me3_Met.bw
H3K27me3_N=/root/vol_1T/Fig3/H3K27me3_N/H3K27me3_N.bw
computeMatrix scale-regions -S $H3K27me3_N $H3K27me3_Met -R $KC_GTF $SAMac_GTF --samplesLabel H3K27me3_N H3K27me3_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName H3K27me3_KC_SAMac.mat.gz
plotProfile --matrixFile H3K27me3_KC_SAMac.mat.gz --outFileName H3K27me3_KC_SAMac.profile.pdf --perGroup --colors '#6a3d9a' '#ff7f00' --plotType se --plotFileFormat pdf

computeMatrix scale-regions -S $H3K27me3_N $H3K27me3_Met -R $KC_GTF --samplesLabel H3K27me3_N H3K27me3_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName H3K27me3_KC.mat.gz
computeMatrix scale-regions -S $H3K27me3_N $H3K27me3_Met -R $SAMac_GTF --samplesLabel H3K27me3_N H3K27me3_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName H3K27me3_SAMac.mat.gz
plotHeatmap --matrixFile H3K27me3_KC.mat.gz --outFileName H3K27me3_KC.heatmap.pdf --plotType se --plotFileFormat pdf --sortUsing sum
plotHeatmap --matrixFile H3K27me3_SAMac.mat.gz --outFileName H3K27me3_SAMac.heatmap.pdf --plotType se --plotFileFormat pdf --sortUsing sum

## H3K27ac
H3K27ac_Met=/root/vol_1T/Fig3/H3K27ac_Met/H3K27ac_Met.bw
H3K27ac_N=/root/vol_1T/Fig3/H3K27ac_N/H3K27ac_N.bw
computeMatrix scale-regions -S $H3K27ac_N $H3K27ac_Met -R $KC_GTF $SAMac_GTF --samplesLabel H3K27ac_N H3K27ac_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName H3K27ac_KC_SAMac.mat.gz
plotProfile --matrixFile H3K27ac_KC_SAMac.mat.gz --outFileName H3K27ac_KC_SAMac.profile.pdf --perGroup --colors '#6a3d9a' '#ff7f00' --plotType se --plotFileFormat pdf

computeMatrix scale-regions -S $H3K27ac_N $H3K27ac_Met -R $KC_GTF --samplesLabel H3K27ac_N H3K27ac_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName H3K27ac_KC.mat.gz
computeMatrix scale-regions -S $H3K27ac_N $H3K27ac_Met -R $SAMac_GTF --samplesLabel H3K27ac_N H3K27ac_Met --upstream 3000 --regionBodyLength 5000 --downstream 3000 --skipZeros --blackListFileName /root/vol_1T/Reference/BlackList/mm10-blacklist.v2.bed --verbose --numberOfProcessors 2 --outFileName H3K27ac_SAMac.mat.gz
plotProfile --matrixFile H3K27ac_KC.mat.gz --outFileName H3K27ac_KC.profile.pdf --perGroup --colors '#6a3d9a' '#ff7f00' --plotType se --plotFileFormat pdf
plotProfile --matrixFile H3K27ac_SAMac.mat.gz --outFileName H3K27ac_SAMac.profile.pdf --perGroup --colors '#6a3d9a' '#ff7f00' --plotType se --plotFileFormat pdf
plotHeatmap --matrixFile H3K27ac_KC.mat.gz --outFileName H3K27ac_KC.heatmap.pdf --plotType se --plotFileFormat pdf --sortUsing sum
plotHeatmap --matrixFile H3K27ac_SAMac.mat.gz --outFileName H3K27ac_SAMac.heatmap.pdf --plotType se --plotFileFormat pdf --sortUsing sum
```