params.cores = 1 
params.memory = '4 GB'

// Params Defaults in juno
params.refGenome = "/work/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta"
params.genomeVersion = "37"
params.circos = "/opt/circos-0.69-2/bin/circos"
params.loci = "/data/copy_number/GermlineHetPon.37.vcf.gz"
params.gcProfile = "/data/copy_number/GC_profile.1000bp.37.cnp"
params.ensemblDataDir = "/data/common/ensembl_data"
params.diploidRegions = "/data/copy_number/DiploidRegions.37.bed.gz"
params.normal = null
params.normalBam = null
params.binProbes = 0
params.binLogR = 0
params.minPurity = 0.08
params.maxPurity = 1.0


def logMessage = """\
    HMFTOOLS - PURPLE
    ========================================
    Running Mode: ${params.normal ? 'Matched' : 'Unmatched'}
    ----------------------------------------
    Params:
    ----------------------------------------
    tumor        : ${params.tumor}
    tumorBam     : ${params.tumorBam}
"""
logMessage += (params.normal && params.normalBam) ? """\
    normal       : ${params.normal}
    normalBam    : ${params.normalBam}
""" : ""
logMessage += """\
    somaticVcf   : ${params.somaticVcf}
    germlineVcf  : ${params.germlineVcf}
    outdir       : ${params.outdir}
    cores        : ${params.cores}
    memory       : ${params.memory}
    binProbes    : ${params.binProbes}
    binLogR      : ${params.binLogR}
    minPurity    : ${params.minPurity}
    maxPurity    : ${params.maxPurity}
    ========================================
    Workflow:
    ----------------------------------------
    Project    : ${workflow.projectDir}
    Cmd line   : ${workflow.commandLine}
"""

log.info(logMessage.stripIndent())

// See https://github.com/hartwigmedical/hmftools/tree/master/amber
process runAmber {
    tag "AMBER on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/amber", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    val tumor
    val normal
    path tumorBam
    path normalBam

    output:
    path "${tumor}.amber.baf.tsv.gz", emit: amber_baf_tsv
    path "${tumor}.amber.baf.pcf", emit: amber_baf_pcf
    path "${tumor}.amber.qc", emit: amber_qc
    path "${tumor}.amber.contamination.vcf.gz", emit: amber_contamination_vcf

    script:
    def reference_args = normal ? """-reference ${normal} \\\n    -reference_bam ${normalBam}"""  : ""
    """
    if [ -f "${params.outdir}/amber/${tumor}.amber.baf.tsv.gz" ] && \\
       [ -f "${params.outdir}/amber/${tumor}.amber.baf.pcf" ] && \\
       [ -f "${params.outdir}/amber/${tumor}.amber.qc" ]  && \\
       [ -f "${params.outdir}/amber/${tumor}.amber.contamination.vcf.gz" ]; then
        echo "Output files already exist. Skipping amber execution."
        ln -fs ${params.outdir}/amber/${tumor}.amber.baf.tsv.gz ${tumor}.amber.baf.tsv.gz
        ln -fs ${params.outdir}/amber/${tumor}.amber.baf.pcf ${tumor}.amber.baf.pcf
        ln -fs ${params.outdir}/amber/${tumor}.amber.qc ${tumor}.amber.qc
        ln -fs ${params.outdir}/amber/${tumor}.amber.contamination.vcf.gz ${tumor}.amber.contamination.vcf.gz
    else
        amber \\
            -tumor ${tumor} \\
            -tumor_bam ${tumorBam} \\
            ${reference_args} \\
            -output_dir \$PWD \\
            -threads ${params.cores} \\
            -loci ${params.loci} \\
            -ref_genome_version V${params.genomeVersion}
    fi
    """.stripIndent()
}

// See https://github.com/hartwigmedical/hmftools/tree/master/cobalt#mandatory-arguments
process runCobalt {
    tag "COBALT on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/cobalt", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    val tumor
    val normal
    path tumorBam
    path normalBam

    output:
    path "${tumor}.cobalt.ratio.tsv.gz", emit: cobalt_tumor_ratio_tsv
    path "${tumor}.cobalt.ratio.pcf", emit: cobalt_tumor_ratio_pcf
    path "${normal}.cobalt.ratio.pcf", emit: cobalt_normal_ratio_pcf, optional: true

    script:
    def reference_args = normal ? "\t\t\t-reference ${normal} \\\n    -reference_bam ${normalBam}"  : "\t\t\t-tumor_only_diploid_bed ${params.diploidRegions}"
    """
    if [ -f "${params.outdir}/cobalt/${tumor}.cobalt.ratio.tsv.gz" ] && \
       [ -f "${params.outdir}/cobalt/${tumor}.cobalt.ratio.pcf" ]; then
        echo "Output files already exist. Skipping cobalt execution."
        ln -s ${params.outdir}/cobalt/${tumor}.cobalt.ratio.tsv.gz ${tumor}.cobalt.ratio.tsv.gz
        ln -s ${params.outdir}/cobalt/${tumor}.cobalt.ratio.pcf ${tumor}.cobalt.ratio.pcf

        if [ -f "${params.outdir}/cobalt/${normal}.cobalt.ratio.pcf" ]; then
            ln -s ${params.outdir}/cobalt/${normal}.cobalt.ratio.pcf ${normal}.cobalt.ratio.pcf
        fi
    else
        cobalt \\
            -tumor ${tumor} \\
            -tumor_bam ${tumorBam} \\
            ${reference_args} \\
            -output_dir \$PWD \\
            -threads ${params.cores} \\
            -gc_profile ${params.gcProfile}
    fi
    """.stripIndent()
}

process binCobalt {
    tag "COBALT BIN on ${params.tumor}" + (params.normal ? " and ${params.normal}" : "")
    publishDir "${params.outdir}/cobalt/binned_${params.binProbes}_probes_${params.binLogR}_LogR", mode: 'copy'
    cpus 1
    memory '4G'
    time '1h'

    input:
    val tumor
    val normal
    val binProbes
    val binLogR
    path cobalt_tumor_ratio_tsv
    path cobalt_tumor_ratio_pcf
    path cobalt_normal_ratio_pcf

    output:
    path "${tumor}.cobalt.ratio.tsv.gz", emit: cobalt_tumor_ratio_tsv
    path "${tumor}.cobalt.ratio.pcf", emit: cobalt_tumor_ratio_pcf
    path "${normal}.cobalt.ratio.pcf", emit: cobalt_normal_ratio_pcf

    script:
    """
    # Bin Cobalt Tumor Probes
    bin_cobalt.py \\
        --in_pcf ${cobalt_tumor_ratio_pcf} \\
        --bin_probes ${binProbes} \\
        --bin_log_r ${binLogR}

    # Bin Cobalt Normal probes
    if [ -f "${cobalt_normal_ratio_pcf}" ]; then
        bin_cobalt.py \\
            --in_pcf ${cobalt_normal_ratio_pcf} \\
            --bin_probes ${binProbes} \\
            --bin_log_r ${binLogR}
    fi
    """.stripIndent()
}

// See https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md#arguments
process runPurple {
    tag "PURPLE on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/purple/purple_${params.minPurity}_${params.maxPurity}", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    val tumor
    val normal
    path somatic_vcf
    path germline_vcf
    path amber_baf_tsv
    path amber_baf_pcf
    path amber_qc
    path amber_contamination_vcf
    path cobalt_tumor_ratio_tsv
    path cobalt_tumor_ratio_pcf
    path cobalt_normal_ratio_pcf
    path cobalt_path

    output:
    path "${tumor}.purple.purity.tsv", emit: purple_purity_tsv
    path "${tumor}.purple.qc", emit: purple_qc
    path "${tumor}.purple.purity.range.tsv", emit: purple_purity_range_tsv
    path "${tumor}.purple.cnv.somatic.tsv", emit: purple_cnv_somatic_tsv
    path "${tumor}.purple.cnv.gene.tsv", emit: purple_cnv_gene_tsv
    path "${tumor}.purple.segment.tsv", emit: purple_segment_tsv
    path "${tumor}.purple.somatic.clonality.tsv", emit: purple_somatic_clonality_tsv
    path "plot/${tumor}.segment.png", emit: purple_segment_png
    path "plot/${tumor}.copynumber.png", emit: purple_copynumber_png
    path "plot/${tumor}.circos.png", emit: purple_circos_png
    path "plot/${tumor}.map.png", emit: purple_map_png
    path "plot/${tumor}.input.png", emit: purple_input_png
    path "plot/${tumor}.purity.range.png", emit: purple_purity_range_png

    script:
    def reference_args = normal ? "-reference ${normal}" : ""
    def somatic_vcf_args = somatic_vcf ? "-somatic_vcf ${somatic_vcf}" : ""
    def germline_vcf_args = germline_vcf ? "-germline_vcf ${germline_vcf}" : ""
    """
    purple \\
        -tumor ${tumor} \\
        ${reference_args} \\
        ${somatic_vcf_args} \\
        ${germline_vcf_args} \\
        -amber ${params.outdir}/amber \\
        -cobalt ${cobalt_path} \\
        -output_dir \$PWD \\
        -gc_profile ${params.gcProfile} \\
        -ref_genome ${params.refGenome} \\
        -ref_genome_version ${params.genomeVersion} \\
        -ensembl_data_dir ${params.ensemblDataDir} \\
        -circos ${params.circos} \\
        -min_purity ${params.minPurity} \\
        -max_purity ${params.maxPurity}

    rsync -a --no-links \$PWD/ ${params.outdir}/purple/
    """.stripIndent()
}

workflow {
    // Arguments
    tumor = Channel.value(params.tumor)
    normal = Channel.value(params.normal)
    tumorBam = Channel.fromPath(params.tumorBam)
    normalBam = Channel.fromPath(params.normalBam)
    binProbes = Channel.value(params.binProbes)
    binLogR = Channel.value(params.binLogR)
    somaticVcf = Channel.value(params.somaticVcf)
    germlineVcf = Channel.value(params.germlineVcf)

    // Run Amber and Cobalt
    amberOutput = runAmber(tumor, normal, tumorBam, normalBam)
    cobaltOutput = runCobalt(tumor, normal, tumorBam, normalBam)

    // Bin Cobalt if expected
    postCobaltOutput = (binProbes != 0 || binLogR != 0)
        ? binCobalt(tumor, normal, binProbes, binLogR, cobaltOutput.cobalt_tumor_ratio_tsv, cobaltOutput.cobalt_tumor_ratio_pcf, cobaltOutput.cobalt_normal_ratio_pcf)
        : cobaltOutput

    cobaltOutdir = (binProbes != 0 || binLogR != 0)
        ? "${params.outdir}/cobalt/binned_${params.binProbes}_probes_${params.binLogR}_LogR"
        : "${params.outdir}/cobalt"

    // Run Purple
    runPurple(
        tumor,
        normal,
        somaticVcf,
        germlineVcf,
        amberOutput.amber_baf_tsv,
        amberOutput.amber_baf_pcf,
        amberOutput.amber_qc,
        amberOutput.amber_contamination_vcf,
        postCobaltOutput.cobalt_tumor_ratio_tsv,
        postCobaltOutput.cobalt_tumor_ratio_pcf,
        postCobaltOutput.cobalt_normal_ratio_pcf,
        cobaltOutdir,
    )
}

workflow.onComplete {
    log.info (
      workflow.success
        ? "\nDone! Purple ran successfully. See the results in: ${params.outdir}\n"
        : "\nOops .. something went wrong\n"
    )
}
