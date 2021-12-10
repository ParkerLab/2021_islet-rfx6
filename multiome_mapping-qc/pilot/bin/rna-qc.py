work_dir="/lab/work/vivekrai/2021-06_pilot-rfx6/work/rnaseq"
libraries=["3528_CV1", "3528_CV2"]
libraries_atac = ["3528_ATAC_CV1", "3528_ATAC_CV2"]
genome="hg19-mCherry-mKate2"

# print("#drmr:job processor_memory=25000 processors=1 node_properties=wolverine")
# for i in libraries:
#     star_path = f"{work_dir}/starsolo/{i}-{genome}"
#     print(f"""/home/vivekrai/analyses/2021-06_pilot-rfx6/control/rnaseq/bin/qc-from-starsolo.py {star_path}/Aligned.sortedByCoord.out.bam \
#                 {star_path}/Solo.out/GeneFull/raw/matrix.mtx \
#                 {star_path}/Solo.out/GeneFull/raw/barcodes.tsv \
#                 > {work_dir}/qc/{i}-{genome}.qc.txt
#     """)


print("#drmr:job processor_memory=100000 processors=1 node_properties=wolverine")
for i in libraries_atac:
    atac_path = "/lab/work/vivekrai/2021-06_pilot-rfx6/work/atacseq/ataqv"

    print(f"""/home/vivekrai/analyses/2021-06_pilot-rfx6/control/atacseq/bin/extractAtaqvMetric.py --files {atac_path}/{i}-{genome}.ataqv.json.gz --metrics \
            tss_enrichment percent_hqaa hqaa total_reads total_autosomal_reads percent_mitochondrial \
            percent_autosomal_duplicate percent_duplicate max_fraction_reads_from_single_autosome \
            | perl -pe 's@.*.ataqv.json.gz\t@{i}-@' > {atac_path}/{i}.metrics.txt
    """)
