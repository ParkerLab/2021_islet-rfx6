nextflow.preview.dsl = 2
params.dev = false

Channel.fromPath("${params.root}/data/bams/*.bam")
    .take(params.dev ? 3 : -1)
    /*
    .map{it ->
        [
            (it.getSimpleName() - /.Aligned.sortedByCoord.out.bam/),
            file(it, checkIfExists: true)
        ]
    }
    */
    .map{
        name = (it.getSimpleName() - /.Aligned.sortedByCoord.out.bam/)
        [
            name,
            "${params.root}/work/2020-01-26_remap/fastq/${name}.1.fq.gz",
            "${params.root}/work/2020-01-26_remap/fastq/${name}.2.fq.gz",
        ]
    }
    .set{ fastqs }


workflow {
    bam_to_fastq(bams) | star | map { it -> 
        [it[0], it[1], it[2]]
    } | (feature_count & count_biotype & tin & qorts)

    merge_counts("exons", feature_count.out.toList().collect())
    merge_counts("biotype", count_biotype.out.toList().collect())
}

def get_star_index(key) {
    return params.star_index[key]
}

def get_gtf(key) {
    return params.gtf[key]
}


process bam_to_fastq {
    storeDir "${params.outdir}/fastq",
    saveAs: {filename -> "${id}.${filename}"}
    input:
        tuple id, path(bam)
    output:
        tuple id, '1.fq.gz', '2.fq.gz'

    """
    samtools sort -n -@ $task.cpus $bam | \
        samtools fastq -1 1.fq -2 2.fq -0 /dev/null -s /dev/null -

    pigz *.fq
    """
}

process star {
    label 'map'
    storeDir "${params.outdir}/star/$id"
    input:
        tuple id, path(fq1), path(fq2)

    output:
        tuple id, "Aligned.sortedByCoord.out.bam", "Aligned.sortedByCoord.out.bam.bai", "Log.*", "SJ.out.tab"

    """
    STAR --genomeDir ${get_star_index('hg19')} \
         --runMode alignReads \
         --sjdbGTFfile ${get_gtf('hg19')} \
         --readFilesIn $fq1 $fq2 \
         --runThreadN $task.cpus \
         --genomeLoad $params.star_genome_load \
         --runRNGseed 2020 \
         --readFilesCommand zcat \
         --outSAMunmapped Within KeepPairs \
         --outSAMtype BAM SortedByCoordinate

    samtools index Aligned.sortedByCoord.out.bam
    """
}

process qorts {
    label 'avg_memory'
    storeDir "${params.outdir}/qorts/"
    input:
        tuple id, path(bam), path(bai)
    output:
        path "$id/*"

    """
    java $params.qorts_java_opts -jar \$(which QoRTs.jar) QC --quiet \
        --generatePlots --randomSeed $params.seed --title $id \
        $bam ${get_gtf('hg19')} $id
    """
}

process tin {
    time 24.h
    storeDir "${params.outdir}/tin/$id"
    input:
        tuple id, path(bam), path(bai)
    output:
        path("*.{txt,xls}")

    """tin.py -i $bam -r ${params.gene_model['hg19']}"""
}

process feature_count {
    memory 4.GB
    time 1.h
    storeDir "${params.outdir}/counts/exon"
    input:
        tuple id, path(bam), path(bai)
    output:
        tuple id, "${id}_featureCounts.txt", "${id}_featureCounts.txt.summary"

    """
    featureCounts -T $task.cpus -p -t $params.fc_count_type \
        -g $params.fc_group_features -A $params.fc_extra_attrib \
        -a ${params.gtf['hg19']} -s $params.fc_strandedness -o ${id}_featureCounts.txt $bam
    """
}


process count_biotype {
    memory 4.GB
    time 1.h
    storeDir "${params.outdir}/counts/biotype"
    input:
        tuple id, path(bam), path(bai)
    output:
        tuple id, "${id}_featureCounts_biotype.txt", "${id}_featureCounts_biotype.txt.summary"

    """
    featureCounts -T $task.cpus -p -g $params.fc_group_biotype \
        -a ${params.gtf['hg19']} -s $params.fc_strandedness -o ${id}_featureCounts_biotype.txt $bam
    """
}


process merge_counts {
    time 15.m
    storeDir "${params.outdir}/counts/collated"
    input:
        val prefix
        val counts
    output:
        file "gene_counts_${prefix}.txt"

    script:
    // 1: gene id
    // ...
    // 7: counts
    gene_ids = "<(tail -n +2 ${counts[0][1]} | cut -f1)"
    counts_list = counts.collect{item ->
        "<(tail -n +2 ${item[1]} | sed 's:Aligned.sortedByCoord.out.bam:${item[0]}:' | cut -f7)"
    }.join(" ")

    """paste $gene_ids $counts_list > gene_counts_${prefix}.txt"""
}
