libraries = params.libraries

bams = []
qc = []
for (lib in libraries) {
    library = lib.key

    bam = file("${params.pruned_bam_dir}/${library}-hg19-mCherry-mKate2.before-dedup.bam", checkIfExists: true)
        index = file("${params.pruned_bam_dir}/${library}-hg19-mCherry-mKate2.before-dedup.bam.bai", checkIfExists: true )
        qcfile = file("${params.qc_dir}/${library}-hg19-mCherry-mKate2.qc.txt", checkIfExists: true)

        bams << [library, bam, index]
        qc << [library, qcfile]
}

autosomes = (1..22).collect({"chr" + it}) + ["chrX"]

// From RNA QC, select nuclei somewhat leniently to reduce the bam file size
// Then Prep bam to run demuxlet and further select singlets.


process select_nuclei {
        /* Removing completely junk nuclei here with min total reads or min UMIs helps reduce the bam
         file size for further steps because demuxlet is very slow with large bam files
         Also split selected nuclei for each library further into chunks of 2500/5000 nuclei with roughly same number of reads */
        publishDir "${params.results}/selected_nuclei", mode: 'copy'
        container params.general_container

        input:
        set val(library), path(qcfile) from Channel.fromList(qc)

        output:
        set val(library), path("${library}.selected_barcodes_*") into select_nuclei_out

        """
        split_droplets.py --qc $qcfile --umi-threshold ${params.umi_threshold} \
        --mitochondrial-fraction ${params.mitochondrial_threshold} \
        --nbarcodes 1000 --out-prefix ${library}.selected_barcodes
        """
}

bamc = Channel.fromList(bams)
filter_bam_in = select_nuclei_out.transpose().map{
    it -> [it[0], it[1].name.replaceFirst(~/$/, ''), it[1] ]
}.combine(bamc, by: 0)


process filter_bam {
        /* Removing completely junk nuclei here with min total reads or min UMIs helps reduce the bam
         file size for further steps
         This process though seems to not take much memory while running, looks like it writes bam at once.
         So that's when it needs more memory */
        publishDir "${params.results}/filter-bam", mode: 'copy'
        errorStrategy 'retry'
        maxRetries 2
        cpus 10
        time { 10.hour * task.attempt }
        memory { 8.GB * task.attempt }
        container params.general_container


        input:
        set val(library), val(library_nuclei), path(selected_nuclei), path(bam), path(bamindex) from filter_bam_in

        output:
        set val(library_nuclei), path("${library_nuclei}.filtered.bam"), path("${library_nuclei}.filtered.bam.bai") into filter_bam_out

        """
        subset-bam --bam $bam \
            --bam-tag CB --cell-barcodes ${selected_nuclei} \
            --cores 10 --log-level info --out-bam ${library_nuclei}.filtered.bam

        samtools index ${library_nuclei}.filtered.bam
        """
}


process prep_vcf {
        /*For RNA, select vcf to contain gencode basic transcripts (intron + exon),
         blacklist regions removed, MAF 5% in the samples of the respective batch. */
        publishDir "${params.results}/vcf", mode: 'copy'
        container params.general_container

        input:
        path(full_vcf) from Channel.fromPath(params.full_vcf)
        path(selected_regions) from Channel.fromPath(params.selected_regions)

        output:
        path("genotypes.maf05.selected_regions.recode.vcf") into prep_vcf_out

        """
        vcftools --vcf ${full_vcf} --maf 0.05 --bed ${selected_regions} --recode --out genotypes.maf05.selected_regions
        """

}


process demuxlet {
        /*Copied over the vcf file to localscratch and reading it from params instead from input
         so it doesn't need to be copied to work directory with the bam file for each run.*/
        errorStrategy 'retry'
        maxRetries 1
        time { 36.hour * task.attempt }
        memory { 16.GB * task.attempt }
        publishDir "${params.results}/demuxlet", mode: 'copy'
        container params.demuxlet_container

        input:
        set val(library_nuclei), path(bam), path(bam_index), path(vcf) from filter_bam_out.combine(prep_vcf_out)

        output:
        path("${library_nuclei}.best") into demuxlet_out
        set val(library_nuclei), path("${library_nuclei}.single"), path("${library_nuclei}.sing2") into demuxlet_out_others

        """
        demuxlet --sam $bam --tag-group CB --tag-UMI UB --vcf ${vcf} --field GT \
            --out ${library_nuclei}
        """

}

process plot_qc {
        errorStrategy 'ignore'
        maxRetries 2
        time { 2.hour * task.attempt }
        memory { 4.GB * task.attempt }
        storeDir "${params.results}/figures"
        container params.general_container

        input:
        path(demux) from demuxlet_out.collect()

        output:
        set path("demuxlet.tsv"), path("fig.demuxlet_qc.png"), path("fig.quality_nuclei_per_sample_by_weight.pdf"), path("fig.demuxlet_qc_faceted_pointdensity.png")

        """
        plot-demuxlet-qc.py --input ${demux.join(" ")} --weights ${params.weights} --data demuxlet.tsv --qc-dir ${params.qc_dir} ;
        plot-demuxlet-qc.R demuxlet.tsv fig.demuxlet_qc_faceted_pointdensity.png
        """
}

workflow.onComplete {
    if (workflow.success){
        subject = "Demuxlet execution complete"
    }
    else {
        subject = "Demuxlet execution error"
    }

    recipient = params.email

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

