work_dir="/lab/work/vivekrai/2021-06_pilot-rfx6/work/atacseq"
libraries_atac = ["3528_ATAC_CV1", "3528_ATAC_CV2"]
genome="hg19-mCherry-mKate2"


print("#drmr:job processor_memory=100000 processors=1 node_properties=wolverine")
for i in libraries_atac:
    print(f"""/home/porchard/github/pltools/bin/extractAtaqvMetric.py \
            --files {work_dir}/ataqv/{i}-{genome}.ataqv.custom.json.gz --metrics tss_enrichment percent_hqaa \
            hqaa total_reads total_autosomal_reads percent_mitochondrial percent_autosomal_duplicate percent_duplicate \
            max_fraction_reads_from_single_autosome > {work_dir}/ataqv/{i}.custom-ref.metrics.txt
    """)
    print(f"""/home/porchard/github/pltools/bin/extractAtaqvMetric.py \
            --files {work_dir}/ataqv/{i}-{genome}.ataqv.json.gz --metrics tss_enrichment percent_hqaa \
            hqaa total_reads total_autosomal_reads percent_mitochondrial percent_autosomal_duplicate percent_duplicate \
            max_fraction_reads_from_single_autosome > {work_dir}/ataqv/{i}.metrics.txt
    """)
