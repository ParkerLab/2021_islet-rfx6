#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'
libraries = params.libraries

bams = []
qc = []
for (library in libraries) {
    bam = new File(params.md_bam_dir, "${library}-hg19.md.bam")
	index = new File(params.md_bam_dir, "${library}-hg19.md.bam.bai")
	qcfile = new File(params.qc_dir, "${library}-hg19.ataqv.txt")
	bams << [library, file(bam), file(index)]
	qc << [library, file(qcfile)]
}


// From RNA QC, select nuclei somewhat leniently to reduce the bam file size
// Then Prep bam to run demuxlet and further select singlets.

process prune {
    /* No need to sort here because chrom sorted fasta was used for bwa, so the bam file is already sorted */
	storeDir "${params.results}/prune-store"
	errorStrategy 'retry'
	maxRetries 2
	cpus 24
	memory { 50.GB * task.attempt }
	time { 36.h * task.attempt }
	container params.general_container
	
	input:
	set val(library), path(md_bam), path(bam_index) from Channel.fromList(bams)

	output:
	set val(library), path("${library}.pruned.bam"), path("${library}.pruned.bam.bai") into prune_out

	"""
	${IONICE} samtools view -h -b -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $md_bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 \
    chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 > ${library}.pruned.bam
    samtools index ${library}.pruned.bam
	"""

}

process select_nuclei {
	/* Removing completely junk nuclei here with min total reads or min UMIs helps reduce the bam
	 file size for further steps because demuxlet is very slow with large bam files
	 Also split selected nuclei for each library further into chunks of 1000/2500/5000 nuclei with roughly same number of reads */
	publishDir "${params.results}/selected_nuclei"
	container params.general_container
	
	input:
	set val(library), path(qcfile) from Channel.fromList(qc)

	output:
	set val(library), path("${library}.selected_barcodes_*") into select_nuclei_out

	"""
    split_droplets.py --qc $qcfile --hqaa-threshold ${params.hqaa_threshold}  --fraction-mitochondrial-threshold ${params.mitochondrial_threshold} \
    --tss-threshold ${params.tss_threshold}  --nbarcodes 1000 --out-prefix ${library}.selected_barcodes
    """
}


filter_bam_in = select_nuclei_out.transpose().map{it -> [it[0], it[1].name.replaceFirst(~/$/, ''), it[1] ]}.combine(prune_out, by: 0)


process filter_bam {
	/* Removing completely junk nuclei here with min total reads or min UMIs helps reduce the bam
	 file size for further steps 
	 This process though seems to not take much memory while running, looks like it writes bam at once. 
	 So that's when it needs more memory */
	publishDir "${params.results}/filter-bam"
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
    subset-bam --bam $bam --bam-tag CB --cell-barcodes ${selected_nuclei} \
     --cores 10 --log-level info --out-bam ${library_nuclei}.filtered.bam;
    samtools index ${library_nuclei}.filtered.bam
    """
}

filter_bam_out.into { demuxlet_batch_in; demuxlet_full_in }

process prep_vcf {
	/*For RNA, select vcf to contain gencode basic transcripts (intron + exon),
	 blacklist regions removed, MAF 5% in the samples of the respective batch. */
	publishDir "${params.results}/vcf"
	container params.general_container
	
	input:
	path(full_vcf) from Channel.fromPath(params.full_vcf)
	path(selected_regions) from Channel.fromPath(params.selected_regions)
	path(sample_list) from Channel.fromPath(params.sample_list)

	output:
	path("genotypes.maf05.selected_regions.batch.recode.vcf") into prep_vcf_out

	"""
    vcftools --gzvcf ${full_vcf} --maf 0.05 --bed ${selected_regions} --keep ${sample_list} --recode --out genotypes.maf05.selected_regions.batch  
    """

}
// demuxlet
process demuxlet_batch {
	/*Copied over the vcf file to localscratch and reading it from params instead from input
	 so it doesn't need to be copied to work directory with the bam file for each run.*/
	errorStrategy 'retry'
	maxRetries 2
	time { 36.hour * task.attempt }
	memory { 8.GB * task.attempt }
	storeDir "${params.results}/demuxlet"
	container params.demuxlet_container
	
	input:
	set val(library_nuclei), path(bam), path(bam_index), path(vcf) from demuxlet_batch_in.combine(prep_vcf_out)
	
	output:
	path("${library_nuclei}.batch.best") into demuxlet_out
	set val(library_nuclei), path("${library_nuclei}.batch.single"), path("${library_nuclei}.batch.sing2") into demuxlet_out_others

	"""
    demuxlet --sam $bam --tag-group CB  --vcf ${vcf} --field GT --out ${library_nuclei}.batch 
    """

}

process plot_qc {
	time { 2.hour * task.attempt }
	memory { 4.GB * task.attempt }
	storeDir "${params.results}/figures"
	container params.general_container
	
	input:
	path(demux) from demuxlet_out.collect()
	
	output:
	set path("demuxlet.tsv"), path("fig.demuxlet_qc_mito.png"), path("fig.demuxlet_qc_tss.png"), path("fig.quality_nuclei_per_sample_by_weight.pdf"), path("fig.demuxlet_qc_faceted_pointdensity_mito.png"), path("fig.demuxlet_qc_faceted_pointdensity_tss.png")

	"""
    plot-demuxlet-qc.py --input ${demux.join(" ")} --weights ${params.weights} --batch ${params.batch} --data demuxlet.tsv --qc-dir ${params.qc_dir} ;
    plot-demuxlet-qc.R demuxlet.tsv fig.demuxlet_qc_faceted_pointdensity_mito.png fig.demuxlet_qc_faceted_pointdensity_tss.png 
    """
	
}

workflow.onComplete {
    if (workflow.success){
        subject = "Demuxlet ATAC execution complete"
    }
    else {
        subject = "Demuxlet ATAC execution error"
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

