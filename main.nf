params.cores = 1 
params.memory = '4 GB'

// Params Defaults
params.genomeVersion = 37
params.circos = "/opt/circos-0.69-2/bin/circos"
params.loci = "/data/copy_number/GermlineHetPon.37.vcf.gz"
params.gcProfile = "/data/copy_number/GC_profile.1000bp.37.cnp"
params.ensemblDataDir = "/data/common/ensembl_data"
params.diploidRegions = "/data/copy_number/DiploidRegions.37.bed.gz"
params.normal = null
params.normalBam = "assets/NO_FILE"
params.binProbes = 0
params.binLogR = 0
params.minPurity = 0.08
params.maxPurity = 1.0

// https://nextflow-io.github.io/patterns/optional-input
def NO_FILE = file("${projectDir}/assets/NO_FILE")
params.normal = null
params.normalBam = NO_FILE


def logMessage = """\
    HMFTOOLS - PURPLE
    ========================================
    Running Mode   : ${params.normal ? 'Matched' : 'Unmatched'}
    ----------------------------------------
    Params:
    ----------------------------------------
    tumor          : ${params.tumor}
    tumorBam       : ${params.tumorBam}
"""
logMessage += (params.normal && params.normalBam) ? """\
    normal        : ${params.normal}
    normalBam     : ${params.normalBam}
""" : ""
logMessage += """\
    outdir         : ${params.outdir}
    cores          : ${params.cores}
    memory         : ${params.memory}
    binProbes      : ${params.binProbes}
    binLogR        : ${params.binLogR}
    minPurity      : ${params.minPurity}
    maxPurity      : ${params.maxPurity}
    ensemblDataDir : ${params.ensemblDataDir}
    ========================================
    Workflow:
    ----------------------------------------
    Project        : ${workflow.projectDir}
    Cmd line       : ${workflow.commandLine}
"""

log.info(logMessage.stripIndent())

// See https://github.com/hartwigmedical/hmftools/tree/master/amber
process runAmber {
    tag "AMBER on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/amber", mode: 'copy', overwrite: false, cache: true
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    path tumorBam
    path normalBam

    output:
    path "${params.tumor}.amber.baf.tsv.gz", emit: amber_baf_tsv
    path "${params.tumor}.amber.baf.pcf", emit: amber_baf_pcf
    path "${params.tumor}.amber.qc", emit: amber_qc
    path "${params.tumor}.amber.contamination.vcf.gz", emit: amber_contamination_vcf, optional: true

    script:
    def reference_args = params.normal ? """\\
            -reference ${params.normal} \\
            -reference_bam ${normalBam} """  : ""

    """
    amber \\
        -tumor ${params.tumor} \\
        -tumor_bam ${tumorBam} ${reference_args} \\
        -output_dir \$PWD \\
        -threads ${params.cores} \\
        -loci ${params.loci} \\
        -ref_genome_version V${params.genomeVersion}
    """.stripIndent()
}

// See https://github.com/hartwigmedical/hmftools/tree/master/cobalt#mandatory-arguments
process runCobalt {
    tag "COBALT on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/cobalt", mode: 'copy', overwrite: false,  cache: true
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    path tumorBam
    path normalBam

    output:
    path "${params.tumor}.cobalt.ratio.tsv.gz", emit: cobalt_tumor_ratio_tsv
    path "${params.tumor}.cobalt.ratio.pcf", emit: cobalt_tumor_ratio_pcf
    path "${params.normal}.cobalt.ratio.pcf", emit: cobalt_normal_ratio_pcf, optional: true

    script:
    def reference_args = params.normal ? """\\
            -reference ${params.normal} \\
            -reference_bam ${normalBam}""" : """\\
            -tumor_only_diploid_bed ${params.diploidRegions}"""

    """
    cobalt \\
        -tumor ${params.tumor} \\
        -tumor_bam ${tumorBam} ${reference_args} \\
        -output_dir \$PWD \\
        -threads ${params.cores} \\
        -gc_profile ${params.gcProfile}
    """.stripIndent()
}

process binCobalt {
    // Only when Unmatched
    tag "COBALT BIN on ${params.tumor}" + (params.normal ? " and ${params.normal}" : "")
    publishDir "${params.outdir}/cobalt/binned_${params.binProbes}_probes_${params.binLogR}_LogR", mode: 'copy'
    cpus 1
    memory '4G'
    time '1h'

    input:
    path cobalt_tumor_ratio_tsv
    path cobalt_tumor_ratio_pcf

    output:
    path "${params.tumor}.cobalt.ratio.tsv.gz", emit: cobalt_tumor_ratio_tsv
    path "${params.tumor}.cobalt.ratio.pcf", emit: cobalt_tumor_ratio_pcf

    script:
    """
    bin_cobalt.py \\
        --in_pcf ${cobalt_tumor_ratio_pcf} \\
        --bin_probes ${params.binProbes} \\
        --bin_log_r ${params.binLogR}
    """.stripIndent()
}

// See https://github.com/hartwigmedical/hmftools/blob/master/sage/README.md#usage
process runSage {
    // Only when matched
    tag "SAGE on ${params.tumor}" + (params.normal ? " vs ${params.normal}" : "")
    publishDir "${params.outdir}/sage", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '4h'

    input:
    path tumorBam
    path normalBam

    output:
    path "${params.tumor}_vs_${params.normal}.vcf.gz", emit: sage_vcf

    script:
    """
    sage \\
        -tumor ${params.tumor} \\
        -tumor_bam ${tumorBam} \\
        -reference ${params.normal} \\
        -reference_bam ${normalBam} \\
        -ref_genome ${params.refGenome} \\
        -ref_genome_version ${params.genomeVersion} \\
        -output_vcf \$PWD/${params.tumor}_vs_${params.normal}.vcf.gz \\
        -threads ${params.cores} \\
        -ensembl_data_dir ${params.ensemblDataDir}
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
    path amber_baf_tsv
    path amber_baf_pcf
    path amber_qc
    path amber_contamination_vcf
    path cobalt_tumor_ratio_tsv
    path cobalt_tumor_ratio_pcf
    path cobalt_normal_ratio_pcf
    path cobalt_path
    path somatic_vcf

    output:
    path "${params.tumor}.purple.purity.tsv", emit: purple_purity_tsv
    path "${params.tumor}.purple.qc", emit: purple_qc
    path "${params.tumor}.purple.purity.range.tsv", emit: purple_purity_range_tsv
    path "${params.tumor}.purple.cnv.somatic.tsv", emit: purple_cnv_somatic_tsv
    path "${params.tumor}.purple.cnv.gene.tsv", emit: purple_cnv_gene_tsv
    path "${params.tumor}.purple.segment.tsv", emit: purple_segment_tsv
    path "${params.tumor}.purple.somatic.clonality.tsv", emit: purple_somatic_clonality_tsv
    path "plot/${params.tumor}.segment.png", emit: purple_segment_png
    path "plot/${params.tumor}.copynumber.png", emit: purple_copynumber_png
    path "plot/${params.tumor}.circos.png", emit: purple_circos_png
    path "plot/${params.tumor}.map.png", emit: purple_map_png
    path "plot/${params.tumor}.input.png", emit: purple_input_png
    path "plot/${params.tumor}.purity.range.png", emit: purple_purity_range_png

    script:
    def reference_args = params.normal ? """\\
        -reference ${params.normal}""" : ""
    def somatic_vcf_args = params.normal && somatic_vcf ? """\\
        -somatic_vcf ${somatic_vcf}""" : ""

    """
    purple \\
        -tumor ${params.tumor} ${reference_args} \\
        -amber ${params.outdir}/amber \\
        -cobalt ${cobalt_path} \\
        -output_dir \$PWD ${somatic_vcf_args} \\
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
    // Input Bams
    tumorBam = params.tumorBam ? Channel.fromPath(params.tumorBam) : error("tumorBam is required")
    normalBam = params.normalBam ? Channel.fromPath(params.normalBam) : Channel.empty()

    // Run Amber and Cobalt
    amberOutput = runAmber(tumorBam, normalBam)
    cobaltOutput = runCobalt(tumorBam, normalBam)
    
    cobaltOutdir = "${params.outdir}/cobalt"

    // // Bin Cobalt, if unmatched and if any bin param is provided
    // if ((!params.normal) && (params.binProbes != 0 || params.binLogR != 0)) {
    //     cobaltOutdir = "${params.outdir}/cobalt"
    // } else {
    //     cobaltOutput = binCobalt(cobaltOutput.cobalt_tumor_ratio_tsv, cobaltOutput.cobalt_tumor_ratio_pcf)
    //     cobaltOutput.cobalt_normal_ratio_pcf = NO_FILE
    //     cobaltOutdir = "${params.outdir}/cobalt/binned_${params.binProbes}_probes_${params.binLogR}_LogR"
    // }

    // Run Sage, if matched
    if (params.normal) {
        sageOutput = runSage(tumorBam, normalBam)
        somatic_vcf = sageOutput.sage_vcf
        // amber_contamination_vcf = amberOutput.amber_contamination_vcf
    } else {
        somatic_vcf = null
        // amber_contamination_vcf = null
    }

    // Run Purple
    runPurple(
        amberOutput.amber_baf_tsv,
        amberOutput.amber_baf_pcf,
        amberOutput.amber_qc,
        amberOutput.amber_contamination_vcf ?: NO_FILE,
        cobaltOutput.cobalt_tumor_ratio_tsv,
        cobaltOutput.cobalt_tumor_ratio_pcf,
        cobaltOutput.cobalt_normal_ratio_pcf ?: NO_FILE,
        cobaltOutdir,
        somatic_vcf ?: NO_FILE,
    )
}

workflow.onComplete {
    log.info (
      workflow.success
        ? "\nDone! Purple ran successfully. See the results in: ${params.outdir}\n"
        : "\nOops .. something went wrong\n"
    )
}
