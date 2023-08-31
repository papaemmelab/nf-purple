#!/bin/sh

TUMOR=$1
TUMOR_BAM=$2
OUTDIR=$3

REFGENOME=/work/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta
REFDIR=/work/isabl/ref/homo_sapiens/37/hmftools/reference_data/
LOCI=${REF_DIR}/copy_number/GermlineHetPon.37.vcf.gz
GCPROFILE=${REF_DIR}/copy_number/GC_profile.1000bp.37.cnp
ENSEMBLDATADIR=${REF_DIR}/common/ensembl_data
GENOMEVERSION=V37
CIRCOS=/opt/circos-0.69-2/bin/circos

nextflow main.nf \
    --tumor ${TUMOR} \
    --tumorBam ${TUMOR_BAM} \
    --outdir ${OUTDIR} \
    --loci ${LOCI} \
    --gcProfile ${GCPROFILE} \
    --ensemblDataDir ${ENSEMBLDATADIR} \
    --genomeVersion ${GENOMEVERSION} \
    --refGenome ${REFGENOME} \
    --circos ${CIRCOS}
