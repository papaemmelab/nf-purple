# Create Test Files

Using `reference.fasta` regions:

```bash
$ cat tests/data/reference.fasta | grep "^>"
>1:100000-200000
>2:300000-400000
```

## Sage and Purple Ensembl Data Files

```bash
DIR_IN=tests/data/ensembl_data_original
DIR_OUT=tests/data/ensembl_data

# For
awk -F, 'NR==1 || ($3 == "1" && $5 >= 100000 && $6 <= 200000) || ($3 == "2" && $5 >= 300000 && $6 <= 400000)' $DIR_IN/ensemble_gene_data.csv > $DIR_OUT/ensemble_gene_data.csv
```

# Sage VCF:

```bash
## Define the input and output files
input_vcf="/data1/papaemme/isabl/home/svc_papaemme_bot/tmp/purple_matched/sage/IID_H211025_T01_01_WG01_vs_IID_H211025_N01_01_WG01.vcf.gz"
output_vcf="tests/data/TEST_TUMOR_vs_TEST_NORMAL.vcf"

## Extract tumor and normal sample names from the filename
filename=$(basename "$input_vcf" .vcf.gz)
tumor_sample=$(echo "$filename" | cut -d'_' -f1-5)
normal_sample=$(echo "$filename" | cut -d'_' -f7-11)

## Filter the VCF and replace sample names
zcat "$input_vcf" | \
awk 'BEGIN {FS=OFS="\t"} /^#/ {print; next} ($1 == "1" && $2 >= 100000 && $2 <= 200000) || ($1 == "2" && $2 >= 300000 && $2 <= 400000) {print}' | \
sed "s/$tumor_sample/TEST_TUMOR/g; s/$normal_sample/TEST_NORMAL/g" > "$output_vcf"

bgzip -f $output_vcf
tabix -f -p vcf $output_vcf.gz
```