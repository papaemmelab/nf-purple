#!/bin/sh

ROOT=/data1/papaemme

TUMOR=IID_H211025_T01_01_WG01
NORMAL=IID_H211025_N01_01_WG01
TUMOR_BAM=`isabl get-bams ${TUMOR}`
NORMAL_BAM=`isabl get-bams ${NORMAL}`
SOMATIC_VCF=/data1/papaemme/isabl/home/liosisk/benchmarking/callers/sage/.sage.vcf.gz

NF_PURPLE=/data1/papaemme/isabl/home/svc_papaemme_bot/dev/nf-purple/main.nf
OUTDIR=/data1/papaemme/isabl/home/svc_papaemme_bot/tmp/purple_matched
REFGENOME=/data1/papaemme/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta
GENOMEVERSION=37

REFDIR=/data1/papaemme/isabl/ref/homo_sapiens/37/hmftools/v6_0/ref/37
LOCI=${REFDIR}/copy_number/GermlineHetPon.37.vcf.gz
GCPROFILE=${REFDIR}/copy_number/GC_profile.1000bp.37.cnp
DIPLOIDREGIONS=${REFDIR}/copy_number/DiploidRegions.37.bed.gz
ENSEMBLDATADIR=${REFDIR}/common/ensembl_data
CIRCOS=/opt/circos-0.69-2/bin/circos

nextflow run \
    -profile slurm \
    ${NF_PURPLE} \
    --tumor ${TUMOR} \
    --tumorBam ${TUMOR_BAM} \
    --normal ${NORMAL} \
    --normalBam ${NORMAL_BAM} \
    --somaticVcf $SOMATIC_VCF \
    --outdir ${OUTDIR} \
    --loci ${LOCI} \
    --gcProfile ${GCPROFILE} \
    --diploidRegions ${DIPLOIDREGIONS} \
    --ensemblDataDir ${ENSEMBLDATADIR} \
    --genomeVersion ${GENOMEVERSION} \
    --refGenome ${REFGENOME} \
    --circos ${CIRCOS} \
    --cores 16 \
    --memory '64G' \
    --binProbes 100 \
    --binLogR 0.5 \
    -resume
