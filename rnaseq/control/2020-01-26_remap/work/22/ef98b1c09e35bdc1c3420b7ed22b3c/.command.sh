#!/bin/bash -ue
samtools sort -n -@ 1 MG-0057_b-ND-53.Aligned.sortedByCoord.out.bam |         samtools fastq -1 1.fq -2 2.fq -0 /dev/null -s /dev/null -

pigz *.fq