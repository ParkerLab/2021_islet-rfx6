rna_libs = {
    "3528_CV1-hg19-mCherry-mKate2": "/lab/work/vivekrai/2021-06_pilot-rfx6/work/rnaseq/starsolo/3528_CV1-hg19-mCherry-mKate2/Aligned.sortedByCoord.out.bam",
    "3528_CV2-hg19-mCherry-mKate2": "/lab/work/vivekrai/2021-06_pilot-rfx6/work/rnaseq/starsolo/3528_CV2-hg19-mCherry-mKate2/Aligned.sortedByCoord.out.bam",
}

for k,v in rna_libs.items():
    print(f"""samtools view  "chr11:2181008:2185439" {rna_libs[k]} | cut -f5""")
