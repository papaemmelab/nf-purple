[![nf-purple CI](https://github.com/papaemmelab/nf-purple/actions/workflows/ci.yml/badge.svg)](https://github.com/papaemmelab/nf-purple/actions/workflows/ci.yml)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)

# 🟣 nf-purple

Nextflow Pipeline to run [Purple](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md#tumor-only-mode) in *Tumor-Only* mode, uses [Amber](https://github.com/hartwigmedical/hmftools/tree/master/amber#tumor-only-mode) and [Cobalt](https://github.com/hartwigmedical/hmftools/tree/master/cobalt#tumor-only-mode) from HMFTools suite, of the Hartwig Foundation.

## 🚀 Run Pipeline

You need Nextflow installed.

```bash
module load java/jdk-11.0.11

nextflow papaemmelab/nf-purple \
    --tumor $tumor \
    --tumor_bam $TUMOR_BAM \
    --outdir $OUTDIR \
    ...refargs
```

## 🧬  Get Reference Data

Downloaded from [Purple Ref Data](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/dna_pipeline) for genome version 37.
More information on their [docs](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md).

## 📒 Tools Info

### Purple

>[Purple](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md) is a **pur**ity **pl**oidy **e**stimator for whole genome sequenced (WGS) data. It combines B-allele frequency (BAF) from AMBER, read depth ratios from COBALT, somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample. PURPLE supports both grch 37 and 38 reference assemblies.

### Amber

>[Amber](https://github.com/hartwigmedical/hmftools/tree/master/amber) is designed to generate a tumor BAF file for use in PURPLE from a provided VCF of likely heterozygous SNP sites.

### Cobalt

>[Cobalt](https://github.com/hartwigmedical/hmftools/tree/master/cobalt). **C**ount **ba**m **l**ines determines the read depth ratios of the supplied **t**umor and reference genomes.

## 🕵🏻‍♂️ Tests

```bash
sh tests/run_test.sh
```

### Docker

[Purple Docker image](https://hub.docker.com/r/papaemmelab/purple) was built for several platforms using [docker buildx](https://docs.docker.com/buildx/working-with-buildx/).

```bash
# Create a new builder instance
docker buildx create --name papaemmelab-builder --use

# Build the image
docker buildx build --platform linux/amd64,linux/arm64 -t papaemmelab/purple:v0.1.0 . --push
```
