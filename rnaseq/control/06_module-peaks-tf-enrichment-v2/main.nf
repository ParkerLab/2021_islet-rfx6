#! /usr/bin/env nextflow

if (params.results == null) {
    println "fatal: key parameter is unset."
    System.exit(1)
}
nextflow.enable.dsl = 2

Channel.fromPath("data/**/hg19.tss.1kb_ext/*.bed")
    .map { it -> [it.getParent().getParent().getName(), it.getParent().getName(), it.getBaseName(), it] }
    .filter { it -> it[2] =~ /Module/ }
    .set { beds }

// Temporary; only because I ran MACS2 pipeline already
Channel.fromPath("work/sciatac/peaks/*_summits.noblacklist.bed")
    .map { it -> [it.getBaseName().split("_")[0], it] }
    .set { summits }

Channel.fromPath(params.motif_db)
    .map { it -> [it.getBaseName(), it] }
    .set { motif_dbs }

Channel.fromPath("data/sciatac_bams/*.bed")
    .map{ it -> [it.getBaseName(), it] }
    .set{ sciatac_bams }

Channel.from(params.extend_bp).set { extend_bp }

Channel.fromPath("data/tss_defs/*.bed")
    .map { it -> [it.getBaseName(), it] }
    .set { tss_defs }

def get_module_res(id) {
    module_res = [
        "Alpha":"/lab/work/vivekrai/2020-01_vanderbilt_rna/work/wgcna-explore/2021-04-07/Alpha.cell.RNA/power-80_cut-0.2/module-assignment.rds",
        "Beta":"/lab/work/vivekrai/2020-01_vanderbilt_rna/work/wgcna-explore/2021-04-07/Beta.cell.RNA/power-80_cut-0.2/module-assignment.rds"
    ]
    return module_res[id]
}

workflow {
    // call_narrow_peaks(sciatac_bams)
    summits.combine(tss_defs).combine(extend_bp)
        | extend_and_annotate_summits
        | module_linked_summits

    shuff_in = module_linked_summits.out.map {
        it -> it[3].collect { nit -> [it[0], it[1], it[2], nit.getBaseName(), nit] }
    }.flatMap { n -> n }

    ame_run_shuffle(shuff_in, motif_dbs)

    // module_linked_summits.out
    // .map {
    //     it -> it[3].collect { nit -> [it[0], it[1], it[2], nit] }
    // }
    // .flatMap { n -> n }
    // .filter { it[3] =~ /NO-MODULE/ }
    // .set { control_seqs }

    // module_linked_summits.out
    // .map {
    //     it -> it[3].collect { nit -> [it[0], it[1], it[2], nit.getBaseName(), nit] }
    // }
    // .flatMap { n -> n }
    // .filter{ !(it[3] =~ /NO-MODULE/) }
    // .combine(control_seqs, by: [0,1,2])
    //     | ame_run_control_seq

}


process call_narrow_peaks_noblacklist {
    storeDir "${params.results}/sciatac/peaks"
    memory '8G'
    time '2h'

    input:
        tuple val(id), path(bed)
    output:
        tuple val(id), path("macs2_out/${id}_summits.bed")
        path("macs2_out/*.{xls,gappedPeak,bdg,narrowPeak}")

    """
    macs2 callpeak -t $bed --outdir macs2_out -n $id -f BED -g hs --nomodel \
        --shift -100 --extsize 200 --seed 2020 --keep-dup all \
        --call-summits

    bedtools intersect -a ${id}_peaks.narrowPeak -b $params.blacklist -v > macs2_out/${id}_peaks.noblacklist.narrowPeak
    bedtools intersect -a ${id}_summits.bed -b $params.blacklist -v > macs2_out/${id}_summits.noblacklist.bed
    """
}

process extend_and_annotate_summits {
    publishDir "${params.results}/annotated_summits/${tss_id}", mode: 'copy'
    input:
        tuple val(id), path(summit_bed), val(tss_id), path(tss_file), val(extend_bp)
    output:
        tuple val(id), val(tss_id), val(extend_bp), path("${id}_summits-ext-${extend_bp}.bed")

    """
    bedtools slop -b $extend_bp -i $summit_bed \
        -g ${params.hg19_sizes} > ${id}_summits_ext.bed

    bedtools intersect -a ${id}_summits_ext.bed \
        -b $tss_file -wa -wb > ${id}_summits-ext-${extend_bp}.bed
    """
}

process module_linked_summits {
    publishDir "${params.results}/tf_enrichment/${id}/${tss_id}/${extend_bp}_bp_ext", mode: 'copy'
    input:
        tuple val(id), val(tss_id), val(extend_bp), path(summit_tss_bed)
    output:
        tuple val(id), val(tss_id), val(extend_bp), path("summit-linked-peaks/*.bed")

    """
    get-module-linked-summits.R --bed $summit_tss_bed \
        --module-res ${get_module_res(id)} --out-dir summit-linked-peaks
    """
}

process ame_run_shuffle {
    publishDir "${params.results}/tf_enrichment/${id}/$tss_id/${extend_bp}_bp_ext/ame_shuff/$mod_id", mode: 'copy'

    input:
        tuple val(id), val(tss_id), val(extend_bp), val(mod_id), path(mod_bed)
        each motif_db
    output:
        tuple path("${motif_db[0]}/*")

    """
    bedtools getfasta -fi ${params.hg19_fasta} -bed $mod_bed > ${mod_id}.fa

    /home/vivekrai/bin/meme/bin/ame --o ${motif_db[0]} --control --shuffle-- ${mod_id}.fa ${motif_db[1]}
    """
}

process ame_run_control_seq {
    publishDir "${params.results}/tf_enrichment/$id/$tss_id/${extend_bp}_bp_ext/ame_control/{$mod_id", mode: "copy"

    input:
        tuple val(id), val(tss_id), val(extend_bp), val(mod_id), path(mod_bed), path(control_bed)
        each motif_db
    output:
        tuple path("${motif_db[0]}/*")

    """
    bedtools getfasta -fi ${params.hg19_fasta} -bed $mod_bed > ${mod_id}.fa
    bedtools getfasta -fi ${params.hg19_fasta} -bed ${control_bed} > control.fa

    /home/vivekrai/bin/meme/bin/ame --o ${motif_db[0]} --control control.fa ${mod_id}.fa ${motif_db[1]}
    """
}
