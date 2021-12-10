#! /usr/bin/env python3

libraries_atac = ["3528_ATAC_CV1", "3528_ATAC_CV2"]
genome="hg19-mCherry-mKate2"


print("#drmr:job processor_memory=100GB processors=1 node_properties=wolverine")
for i in libraries_atac:
    atac_path = "/lab/work/vivekrai/2021-06_pilot-rfx6/work/atacseq/ataqv"

    print(f"""/home/vivekrai/analyses/2021-06_pilot-rfx6/control/atacseq/bin/extractChromosomeCounts.py --files {atac_path}/{i}-{genome}.ataqv.custom.json.gz > {atac_path}/{i}.chrom.metrics.txt
    """)
