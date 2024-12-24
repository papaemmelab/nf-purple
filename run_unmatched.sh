#!/bin/sh

ROOT=/data1/papaemme

TUMOR=IID_H211025_T01_01_WG01
TUMOR_BAM=`isabl get-bams ${TUMOR}`

NF_PURPLE=/data1/papaemme/isabl/home/svc_papaemme_bot/dev/nf-purple/main.nf
OUTDIR=/data1/papaemme/isabl/home/svc_papaemme_bot/tmp/purple_unmatched
REFGENOME=/data1/papaemme/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta
GENOMEVERSION=37

REFDIR=/data1/papaemme/isabl/ref/homo_sapiens/37/hmftools/v5_33
LOCI=${REFDIR}/copy_number/GermlineHetPon.37.vcf.gz
GCPROFILE=${REFDIR}/copy_number/GC_profile.1000bp.37.cnp
DIPLOIDREGIONS=${REFDIR}/copy_number/DiploidRegions.37.bed.gz
ENSEMBLDATADIR=${REFDIR}/common/ensembl_data
CIRCOS=/opt/circos-0.69-2/bin/circos

mkdir -p ${OUTDIR}

cd ${OUTDIR}

nextflow run \
    -profile slurm \
    ${NF_PURPLE} \
    --tumor ${TUMOR} \
    --tumorBam ${TUMOR_BAM} \
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
