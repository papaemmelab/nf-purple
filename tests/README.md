# Test data

```bash
# Build test BAM file
TUMOR_ID=IID_H195923_T01_01_WG01
ORIGINAL_TUMOR_BAM=`isabl get-bams $TUMOR_ID`
TEST_TUMOR_BAM=tumor.bam

NORMAL_ID=IID_H195923_N01_01_WG01
ORIGINAL_NORMAL_BAM=`isabl get-bams $NORMAL_ID`
TEST_NORMAL_BAM=normal.bam

REGION1=1:100000-200000
REGION2=2:300000-400000

samtools view -b -o $TEST_TUMOR_BAM $ORIGINAL_TUMOR_BAM $REGION1 $REGION2
samtools view -b -o $TEST_NORMAL_BAM $ORIGINAL_NORMAL_BAM $REGION1 $REGION2

samtools index $TEST_TUMOR_BAM
samtools index $TEST_NORMAL_BAM

# Build test Reference files
ORIGINAL_REF=`isabl get-reference GRCh37`
LOCAL_REF=gr37.fasta
TEST_REF=reference.fasta

cp $ORIGINAL_REF $LOCAL_REF
cp $ORIGINAL_REF.fai $LOCAL_REF.fai

samtools faidx $LOCAL_REF $REGION1 >> $TEST_REF
samtools faidx $LOCAL_REF $REGION2 >> $TEST_REF
samtools faidx $TEST_REF
samtools dict -a GRCh37 -s HUMAN $LOCAL_REF > $TEST_REF.dict
bwa index $TEST_REF
```
