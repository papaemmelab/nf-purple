params.cores = 4
params.memory = '32 GB'

// Params Defaults in juno
params.refGenome = "/work/isabl/ref/homo_sapiens/GRCh37d5/gr37.fasta"
params.genomeVersion = "V37"
params.circos = "/opt/circos-0.69-2/bin/circos"
params.loci = "/data/copy_number/GermlineHetPon.37.vcf.gz"
params.gcProfile = "/data/copy_number/GC_profile.1000bp.37.cnp"
params.ensemblDataDir = "/data/common/ensembl_data"
params.diploidRegions = "/data/copy_number/DiploidRegions.37.bed.gz"
params.binProbes = 0
params.binLogR = 0
params.minPurity = 0.08
params.maxPurity = 1.0


log.info """\
    HMFTOOLS - PURPLE
    ========================================
    Params:
    ----------------------------------------
    tumor        : ${params.tumor}
    tumorBam     : ${params.tumorBam}
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
    .stripIndent()


process runAmber {
    tag "AMBER on ${params.tumor}"
    publishDir "${params.outdir}/amber", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    val tumor
    path tumorBam

    output:
    path "${tumor}.amber.baf.tsv.gz", emit: amber_baf_tsv
    path "${tumor}.amber.baf.pcf", emit: amber_baf_pcf
    path "${tumor}.amber.qc", emit: amber_qc

    script:
    """
    if [ -f "${params.outdir}/amber/${tumor}.amber.baf.tsv.gz" ] && \
       [ -f "${params.outdir}/amber/${tumor}.amber.baf.pcf" ] && \
       [ -f "${params.outdir}/amber/${tumor}.amber.qc" ]; then
        echo "Output files already exist. Skipping amber execution."
        ln -s ${params.outdir}/amber/${tumor}.amber.baf.tsv.gz ${tumor}.amber.baf.tsv.gz
        ln -s ${params.outdir}/amber/${tumor}.amber.baf.pcf ${tumor}.amber.baf.pcf
        ln -s ${params.outdir}/amber/${tumor}.amber.qc ${tumor}.amber.qc
    else
        amber \
            -tumor ${tumor} \
            -tumor_bam ${tumorBam} \
            -output_dir \$PWD \
            -threads ${params.cores} \
            -loci ${params.loci} \
            -ref_genome_version ${params.genomeVersion}
    fi
    """.stripIndent()
}

process runCobalt {
    tag "COBALT on ${params.tumor}"
    publishDir "${params.outdir}/cobalt", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    val tumor
    path tumorBam

    output:
    path "${tumor}.cobalt.ratio.tsv.gz", emit: cobalt_ratio_tsv
    path "${tumor}.cobalt.ratio.pcf", emit: cobalt_ratio_pcf

    script:
    """
    if [ -f "${params.outdir}/cobalt/${tumor}.cobalt.ratio.tsv.gz" ] && \
       [ -f "${params.outdir}/cobalt/${tumor}.cobalt.ratio.pcf" ]; then
        echo "Output files already exist. Skipping cobalt execution."
        ln -s ${params.outdir}/cobalt/${tumor}.cobalt.ratio.tsv.gz ${tumor}.cobalt.ratio.tsv.gz
        ln -s ${params.outdir}/cobalt/${tumor}.cobalt.ratio.pcf ${tumor}.cobalt.ratio.pcf
    else
        cobalt \
        -tumor ${tumor} \
        -tumor_bam ${tumorBam} \
        -output_dir \$PWD \
        -threads ${params.cores} \
        -gc_profile ${params.gcProfile} \
        -tumor_only_diploid_bed ${params.diploidRegions}
    fi
    """.stripIndent()
}

process binCobalt {
    tag "COBALT BIN on ${params.tumor}"
    publishDir "${params.outdir}/cobalt/binned_${params.binProbes}_probes_${params.binLogR}_LogR", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    val tumor
    val binProbes
    val binLogR
    path cobalt_ratio_pcf
    path cobalt_ratio_tsv

    output:
    path "${tumor}.cobalt.ratio.tsv.gz", emit: cobalt_ratio_tsv
    path "${tumor}.cobalt.ratio.pcf", emit: cobalt_ratio_pcf

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import numpy as np
 
    cobalt_ratio_pcf = pd.read_csv('${cobalt_ratio_pcf}', sep='\\t')
    cobalt_ratio_pcf_probes = pd.DataFrame(columns=cobalt_ratio_pcf.columns)

    # First bin by probes
    chrom_arm = None
    last_idx = None
    for idx, seg in cobalt_ratio_pcf.iterrows():
        if chrom_arm != '_'.join(seg[['chrom','arm']]):
            chrom_arm = '_'.join(seg[['chrom','arm']])
            cobalt_ratio_pcf_probes = cobalt_ratio_pcf_probes.append(seg, ignore_index=True)
            last_idx = cobalt_ratio_pcf_probes.index[-1]
            continue
        if (
            cobalt_ratio_pcf_probes.loc[last_idx, 'n.probes'] <= ${binProbes}
        ) or (
            seg['n.probes'] <= ${binProbes}
        ):
            means = [cobalt_ratio_pcf_probes.loc[last_idx, 'mean']] * cobalt_ratio_pcf_probes.loc[last_idx, 'n.probes']
            means.extend([seg['mean']] * seg['n.probes'])
            cobalt_ratio_pcf_probes.loc[last_idx, 'mean'] = np.mean(means)
            cobalt_ratio_pcf_probes.loc[last_idx, 'n.probes'] += seg['n.probes']
            cobalt_ratio_pcf_probes.loc[last_idx, 'end.pos'] = seg['end.pos']
        else:
            cobalt_ratio_pcf_probes = cobalt_ratio_pcf_probes.append(seg, ignore_index=True)
            last_idx = cobalt_ratio_pcf_probes.index[-1]

    # Then bin by logR mean
    cobalt_ratio_pcf_probes = cobalt_ratio_pcf_probes.reset_index().drop(columns="index")
    cobalt_ratio_pcf_probes_logR = pd.DataFrame(columns=cobalt_ratio_pcf_probes.columns)
    chrom_arm = None
    for idx, seg in cobalt_ratio_pcf_probes.iterrows():
        if chrom_arm != '_'.join(seg[['chrom','arm']]):
            chrom_arm = '_'.join(seg[['chrom','arm']])
            cobalt_ratio_pcf_probes_logR = cobalt_ratio_pcf_probes_logR.append(seg, ignore_index=True)
            last_idx = cobalt_ratio_pcf_probes_logR.index[-1]
            continue
        if abs(cobalt_ratio_pcf_probes.loc[last_idx, 'mean'] - seg['mean']) <= ${binLogR}:
            means = [cobalt_ratio_pcf_probes_logR.loc[last_idx, 'mean']] * cobalt_ratio_pcf_probes_logR.loc[last_idx, 'n.probes']
            means.extend([seg['mean']] * seg['n.probes'])
            cobalt_ratio_pcf_probes_logR.loc[last_idx, 'mean'] = np.mean(means)
            cobalt_ratio_pcf_probes_logR.loc[last_idx, 'n.probes'] += seg['n.probes']
            cobalt_ratio_pcf_probes_logR.loc[last_idx, 'end.pos'] = seg['end.pos']
        else:
            cobalt_ratio_pcf_probes_logR = cobalt_ratio_pcf_probes_logR.append(seg, ignore_index=True)
            last_idx = cobalt_ratio_pcf_probes_logR.index[-1]
 
    cobalt_ratio_pcf_probes_logR.to_csv("${tumor}.cobalt.ratio.pcf", sep='\\t', index=False)
    """
}

process runPurple {
    tag "PURPLE on ${params.tumor}"
    publishDir "${params.outdir}/purple/purple_${params.minPurity}_${params.maxPurity}", mode: 'copy'
    cpus params.cores
    memory params.memory
    time '1h'

    input:
    val tumor
    path amber_baf_tsv
    path amber_baf_pcf
    path amber_qc
    path cobalt_ratio_tsv
    path cobalt_ratio_pcf
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
    """
    purple \
        -tumor ${tumor} \
        -amber ${params.outdir}/amber \
        -cobalt ${cobalt_path} \
        -output_dir \$PWD \
        -gc_profile ${params.gcProfile} \
        -ref_genome ${params.refGenome} \
        -ref_genome_version ${params.genomeVersion} \
        -ensembl_data_dir ${params.ensemblDataDir} \
        -circos ${params.circos} \
        -min_purity ${params.minPurity} \
        -max_purity ${params.maxPurity}
    """.stripIndent()
}

workflow {
    tumor = Channel.value(params.tumor)
    tumorBam = Channel.fromPath(params.tumorBam)
    binProbes = Channel.value(params.binProbes)
    binLogR = Channel.value(params.binLogR)

    amberOutput = runAmber(tumor, tumorBam)
    cobaltOutput = runCobalt(tumor, tumorBam)

    if (binProbes != 0 || binLogR != 0) {
        binCobaltOutput = binCobalt(tumor, binProbes, binLogR, cobaltOutput.cobalt_ratio_pcf, cobaltOutput.cobalt_ratio_tsv)
        runPurple(tumor, amberOutput, binCobaltOutput, "${params.outdir}/cobalt/binned_${params.binProbes}_probes_${params.binLogR}_LogR")
    } else {
        runPurple(tumor, amberOutput, cobaltOutput, "${params.outdir}/cobalt")
    }
}

