#!/bin/sh

TUMOR=IID_H208153_T01_01_WG01
TUMOR_BAM=`isabl get-bams ${TUMOR}`
OUTDIR=/work/isabl/home/arangooj/run/purple/${TUMOR}

REFGENOME=/work/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta
GENOMEVERSION=V37
REFDIR=/work/isabl/ref/homo_sapiens/37/hmftools
LOCI=${REFDIR}/copy_number/GermlineHetPon.37.vcf.gz
GCPROFILE=${REFDIR}/copy_number/GC_profile.1000bp.37.cnp
DIPLOIDREGIONS=${REFDIR}/copy_number/DiploidRegions.37.bed.gz
ENSEMBLDATADIR=${REFDIR}/common/ensembl_data
CIRCOS=/opt/circos-0.69-2/bin/circos

nextflow run \
    -profile hpc \
    /work/isabl/home/arangooj/dev/nf-purple/main.nf \
    --tumor ${TUMOR} \
    --tumorBam ${TUMOR_BAM} \
    --outdir ${OUTDIR} \
    --loci ${LOCI} \
    --gcProfile ${GCPROFILE} \
    --diploidRegions ${DIPLOIDREGIONS} \
    --ensemblDataDir ${ENSEMBLDATADIR} \
    --genomeVersion ${GENOMEVERSION} \
    --refGenome ${REFGENOME} \
    --circos ${CIRCOS}
