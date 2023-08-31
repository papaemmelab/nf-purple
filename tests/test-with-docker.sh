#!/bin/bash

# Run tox inside the container.
echo "testing container..."

# make sure current dir is repo dir
cd $( dirname "${BASH_SOURCE[0]}" )

TEST_IMAGE="nf-heme-classifiers_test_image"

if [ "$1" = "--skip-build" ]; then
    echo "skipping build..."
else
    echo "building image, to skip run with --skip-build..."
    docker build -q -t $TEST_IMAGE .
fi

# see https://explainshell.com/explain?cmd=set+-euxo%20pipefail
set -euxo pipefail

# remove pytest cache
echo "testing docker image..."

# create alias
alias nextflow="docker run --rm $TEST_IMAGE"
alias nf-test="docker run --entrypoint '' test_nextflow /opt/bin/nf-test"

nextflow \
  -dockerize ../main.nf \
    -profile cloud \
    -resume \
    --sample TEST_01 \
    --counts data/TEST_01.ReadsPerGene.out.tab \
    --outdir outdir

echo "tests finished..."



# Ideas:
# https://github.com/nf-core/demultiplex/commit/709d21dd8aa9b2eb3aba376c8f951e97d52233e6
#
