#! /usr/bin/env nextflow

nextflow.enable.dsl = 2

params.atac_barcodes = "/lab/work/vivekrai/2021-06_pilot-rfx6/work/nuclei-qc/atac/atac_barcodes_post-qc.txt"
params.rna_barcodes = "/lab/work/vivekrai/2021-06_pilot-rfx6/work/nuclei-qc/rna/rna_barcodes_post-qc.txt"

Channel.fromPath("/lab/work/vivekrai/2021-06_pilot-rfx6/work/atacseq/gene-counts/*.txt")
    .map { it -> [it.getBaseName().split("-")[0].replaceAll(/_ATAC/, ""), it] }
    .set { atac_libs }

Channel.fromPath("/lab/work/vivekrai/2021-06_pilot-rfx6/work/rnaseq/starsolo/*", type: 'dir')
    .map { it -> [it.getBaseName().split("-")[0], file(it + "/Solo.out")] }
    .set { rna_libs }


workflow {
    get_atac_counts(atac_libs)
    get_rna_counts(rna_libs)
}

process get_atac_counts {
    storeDir "${params.results}/atac"
    memory '10 GB'
    time '1h'

    input:
    tuple val(library), path(counts)

    output:
    path("*.hdf5")

    """
    grep -w $library ${params.atac_barcodes} | cut -d ' ' -f1 > keep-barcodes.txt
    atac-peak-counts-to-hdf5.py --counts $counts --barcodes keep-barcodes.txt --out ${library}.hdf5
    """
}

process get_rna_counts {
    storeDir "${params.results}/rna"
    memory '10 GB'
    time '1h'

    input:
    tuple val(library), path(solo_dir)

    output:
    path("*.hdf5")

    """
    grep -w $library ${params.rna_barcodes} | cut -d ' ' -f1 > keep-barcodes.txt
    starsolo-counts-to-hdf5-gene-name.py --solo-dir $solo_dir \
        --barcodes keep-barcodes.txt --library $library --out ${library}.hdf5
    """
}
