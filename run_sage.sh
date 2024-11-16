#!/bin/sh

TUMOR=IID_H211025_T01_01_WG01
NORMAL=IID_H211025_N01_01_WG01
TUMOR_BAM=`isabl get-bams ${TUMOR}`
NORMAL_BAM=`isabl get-bams ${NORMAL}`
REFGENOME=/data1/papaemme/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta
GENOMEVERSION=V37
REFDIR=/data1/papaemme/isabl/home/liosisk/benchmarking/callers/sage
OUTDIR=/data1/papaemme/isabl/home/svc_papaemme_bot/tmp/purple_matched/sage_run
CORES=16

mkdir -p ${OUTDIR}

# singularity run \
#     --bind /data1:/data1 \
#     --bind /scratch:/scratch \
#     --bind /usersoftware:/usersoftware \
#     /data1/papaemme/isabl/home/liosisk/images/sage.sif \
#     java -Xms4G -Xmx64G -cp /sage_v3.0_beta.jar com.hartwig.hmftools.sage.SageApplication \
#         -threads ${CORES} \
#         -tumor ${TUMOR} \
#         -tumor_bam ${TUMOR_BAM} \
#         -reference ${NORMAL} \
#         -reference_bam ${NORMAL_BAM} \
#         -ref_genome_version ${GENOMEVERSION} \
#         -ref_genome ${REFGENOME} \
#         -hotspots ${REFDIR}/KnownHotspots.somatic.37.vcf.gz \
#         -panel_bed ${REFDIR}/ActionableCodingPanel.somatic.37.bed.gz \
#         -high_confidence_bed ${REFDIR}/GIAB-High-Conf/37/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed \
#         -ensembl_data_dir ${REFDIR}/Ensembl-Data-Cache \
#         -out "${REFDIR}/${TNAME}.sage.vcf.gz"


singularity run \
    --bind /data1:/data1 \
    --bind /scratch:/scratch \
    --bind /usersoftware:/usersoftware \
    /data1/papaemme/isabl/home/liosisk/images/sage.sif \
    java -Xms4G -Xmx64G -cp /sage_v3.0_beta.jar com.hartwig.hmftools.sage.SageApplication \
        -tumor ${TUMOR} \
        -tumor_bam ${TUMOR_BAM} \
        -reference ${NORMAL} \
        -reference_bam ${NORMAL_BAM} \
        -ref_genome_version ${GENOMEVERSION} \
        -ref_genome ${REFGENOME} \
        -output_vcf "${OUTDIR}/${TUMOR}.sage.vcf.gz"