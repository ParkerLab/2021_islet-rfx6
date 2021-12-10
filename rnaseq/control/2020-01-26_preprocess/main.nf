nextflow.preview.dsl = 2
params.dev = false

Channel.fromPath("${params.root}/data/bams/*.bam")
    .take(params.dev ? 3 : -1)
    .map{it ->
        [
            (it.getSimpleName() - /.Aligned.sortedByCoord.out.bam/),
            file(it, checkIfExists: true)
        ]
    }.set{ bams }


workflow {
    // Run QoRTs pre and post-prune
    prune(bams)
    bams | mix(prune.out) | map {it ->
        // if filename has 'prune' in it; it is already pruned
        if(it[1] =~ /prune/) [it[0], 'post-prune', it[1]]
        // else not pruned
        else [it[0], 'pre-prune', it[1]]
    } | qorts
}


process prune {
    publishDir "${params.outdir}/prune", mode: 'copy',
    saveAs: {filename -> "${id}.$filename"}
    input:
        tuple id, path(bam)
    output:
        tuple id, 'pruned.bam'

    """
    samtools view -h -b -q 255   \
        -@ $task.cpus            \
        -f $params.include_flags \
        -F $params.exclude_flags \
        -o pruned.bam            \
        $bam
    """
}

process qorts {
    memory "14 G"

    publishDir "${params.outdir}/qorts/$stage/$id", mode: 'copy'
    input:
        tuple id, stage, path(bam)
    output:
        path "$id/*"

    """
    java $params.qorts_java_opts -jar \$(which QoRTs.jar) QC --quiet \
        --generatePlots --randomSeed $params.seed --title $id \
        $bam $params.gtf $id
    """
}
