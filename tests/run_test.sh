#!/bin/bash

# nextflow -dockerize ../main.nf \
#   -profile cloud \
#   -resume \
#   --sample TEST_01 \
#   --counts data/TEST_01.ReadsPerGene.out.tab \
#   --outdir outdir

docker pull -q nextflow/nextflow:latest

docker run --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v `pwd`:`pwd` \
  -w `pwd` \
  -e "NXF_HOME=`pwd`" \
  nextflow/nextflow:latest \
  nextflow run `pwd`/main.nf \
    -profile cloud \
    -resume \
    --sample TEST_01 \
    --counts `pwd`/tests/data/TEST_01.ReadsPerGene.out.tab \
    --outdir `pwd`/outdir

cat `pwd`/outdir/TEST_01_predictions.tsv
