[![nf-purple CI](https://github.com/papaemmelab/nf-purple/actions/workflows/ci.yml/badge.svg)](https://github.com/papaemmelab/nf-purple/actions/workflows/ci.yml)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)
# nf-purple

Nextflow Pipeline to run Purple from [HmfTools](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md#tumor-only-mode).

## Run Pipeline

You need Nextflow installed.

```bash
module load java/jdk-11.0.11

nextflow papaemmelab/nf-purple \
    --tumor $tumor \
    --tumor_bam $TUMOR_BAM \
    --outdir $OUTDIR
```

### Tests

```bash
sh tests/run_test.sh
```
# nf-purple
