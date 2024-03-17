ngm_dir="/home/dell/miniconda3/bin/"
cutadapt_dir="/home/dell/miniconda3/bin/"
samblaster_dir="/home/dell/miniconda3/bin/"
bedtools_dir="/home/dell/miniconda3/bin/"
trimmomatic_dir="/home/dell/Documents/Trimmomatic-0.39/"
samtools_dir="/home/dell/miniconda3/bin/"
htslib_dir="/home/dell/miniconda3/bin/"
picard_dir="/home/dell/Documents/Picard.2.19.0/"
slamdunk_dir="/home/dell/miniconda3/bin/"
varscan_dir="/home/dell/Documents/VarScan.v2.3.9/"
seqtk_dir="/home/dell/miniconda3/bin/"
java_dir="/home/dell/miniconda3/bin/"
PATH=$samtools_dir:$PATH


refseq="/home/dell/Downloads/SLAM_seq_pipeline/GRCm38.primary_assembly.genome.fa"
#three_prime_utr_bed="ref.genome.three_prime_utr.bed"
min_map_qual=2
min_aln_identity=0.8
min_edit_distance=-1
min_cov_for_snp=10
min_var_freq=0.2
threads="24"

sample_id="MC38_KO1"
R1_reads="/home/dell/Downloads/novogene/P23020603/rawdata/MC38_KO1_Met_R1.fq.gz"
R2_reads="/home/dell/Downloads/novogene/P23020603/rawdata/MC38_KO1_Met_R2.fq.gz"

debug="no"

adapters="/home/dell/Documents/Reference/trimmomatic_adapters/*.fa"
cat $adapters |sed 's/>/\n>/' > adapter.fa

mkdir tmp

# index reference sequence
if [[ ! -e $refseq.fai ]]
then
    $samtools_dir/samtools faidx $refseq
fi

if [[ ! -e ref.genome.dict ]]
then
    $java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary  -REFERENCE $refseq -OUTPUT ref.genome.dict
fi

# trim the reads
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $R1_reads $R2_reads $sample_id.R1.trimmed.PE.fq.gz $sample_id.R1.trimmed.SE.fq.gz $sample_id.R2.trimmed.PE.fq.gz $sample_id.R2.trimmed.SE.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36

# reverse complement R2 reads
$seqtk_dir/seqtk seq -r $sample_id.R2.trimmed.PE.fq.gz |gzip -c > $sample_id.R2.trimmed.revcom.PE.fq.gz

$ngm_dir/ngm \
    -r $refseq \
    --qry1 $sample_id.R1.trimmed.PE.fq.gz \
    --qry2 $sample_id.R2.trimmed.revcom.PE.fq.gz \
    -t $threads \
    -n 1 \
    --strata \
    --bam \
    --slam-seq 2 \
    --no-progress \
    -o $sample_id.bam

# $samtools_dir/samtools index $sample_id.bam
$samtools_dir/samtools view -h -@ $threads $sample_id.bam | $samblaster_dir/samblaster | $samtools_dir/samtools sort -@ $threads -O bam -T $sample_id - >$sample_id.sort.bam
$samtools_dir/samtools index $sample_id.sort.bam

if [[ $debug == "no" ]]
then
    rm $sample_id.R1.trimmed.SE.fq.gz
    rm $sample_id.R2.trimmed.SE.fq.gz
    rm $sample_id.bam
fi

$slamdunk_dir/slamdunk filter -t $threads -mq $min_map_qual -mi $min_aln_identity -nm $min_edit_distance $sample_id.sort.bam -o  $sample_id.slamdunk_out

$samtools_dir/samtools mpileup -B -A -f $refseq --output-QNAME ./$sample_id.slamdunk_out/$sample_id.sort_filtered.bam  > $sample_id.filtered.mpileup 

$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $varscan_dir/VarScan.jar mpileup2snp $sample_id.filtered.mpileup  --strand-filter 1 --output-vcf 1 --min-var-freq $min_var_freq --min-coverage $min_cov_for_snp -t $threads --variants 1 > $sample_id.vcf

$htslib_dir/bgzip $sample_id.vcf
$htslib_dir/tabix -p vcf $sample_id.vcf.gz

perl vcf2reads.pl -m $sample_id.filtered.mpileup -v $sample_id.vcf.gz -o $sample_id.read_list.txt 

$seqtk_dir/seqtk subseq $R1_reads $sample_id.read_list.txt |gzip -c  > $sample_id.base_converted.R1.fastq.gz
$seqtk_dir/seqtk subseq $R2_reads $sample_id.read_list.txt |gzip -c  > $sample_id.base_converted.R2.fastq.gz

# trim the reads
$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $R1_reads $R2_reads $sample_id.base_converted.R1.trimmed.PE.fq.gz $sample_id.base_converted.R1.trimmed.SE.fq.gz $sample_id.base_converted.R2.trimmed.PE.fq.gz $sample_id.base_converted.R2.trimmed.SE.fq.gz ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36

# reverse complement R2 reads
$seqtk_dir/seqtk seq -r $sample_id.base_converted.R2.trimmed.PE.fq.gz |gzip -c > $sample_id.base_converted.R2.trimmed.revcom.PE.fq.gz

$ngm_dir/ngm \
    -r $refseq \
    --qry1 $sample_id.base_converted.R1.trimmed.PE.fq.gz \
    --qry2 $sample_id.base_converted.R2.trimmed.revcom.PE.fq.gz \
    -t $threads \
    -n 1 \
    --strata \
    --bam \
    --slam-seq 2 \
    --no-progress \
    -o $sample_id.base_converted.bam

# $samtools_dir/samtools index $sample_id.base_converted.bam
$samtools_dir/samtools view -h -@ $threads $sample_id.base_converted.bam | $samblaster_dir/samblaster | $samtools_dir/samtools sort -@ $threads -O bam -T $sample_id.base_converted - >$sample_id.base_converted.sort.bam
$samtools_dir/samtools index $sample_id.base_converted.sort.bam


if [[ $debug == "no" ]]
then
    rm $sample_id.base_converted.R1.trimmed.SE.fq.gz
    rm $sample_id.base_converted.R2.trimmed.SE.fq.gz
    rm -r $sample_id.sort.bam
    rm -r $sample_id.sort.bam.bai
    rm -r $sample_id.slamdunk_out
    rm $sample_id.vcf.gz
    rm $sample_id.vcf.gz.tbi
    rm $sample_id.filtered.mpileup
    rm $sample_id.base_converted.bam
    rm $sample_id.read_list.txt
    rm $sample_id.R1.trimmed.PE.fq.gz
    rm $sample_id.R2.trimmed.PE.fq.gz
    rm $sample_id.R2.trimmed.revcom.PE.fq.gz
    rm $sample_id.base_converted.R1.trimmed.PE.fq.gz
    rm $sample_id.base_converted.R2.trimmed.PE.fq.gz
    rm $sample_id.base_converted.R2.trimmed.revcom.PE.fq.gz
fi

mkdir $sample_id
mv $sample_id.* ./$sample_id


rm -r tmp
