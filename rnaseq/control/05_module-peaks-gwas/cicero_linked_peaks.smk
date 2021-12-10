from functools import partial
from datetime import date

run_date = "2021-06-09"
work = partial(os.path.join, "work", run_date)

# CUTOFFS = ["0.15", "0.2", "0.25", "0.3"]
CUTOFFS = ["0.2"]

# CONNS = {
#     "beta": "/lab/work/vivekrai/2019_sciATAC_new/work/2019-07-26_Cicero-2/beta/beta_gt-0.05.filtered.txt",
#     "alpha": "/lab/work/vivekrai/2019_sciATAC_new/work/2019-07-26_Cicero-2/alpha/alpha_gt-0.05.filtered.txt",
# }

# MODULE_RES = {
#   "beta": "/lab/work/vivekrai/2020-01_vanderbilt_rna/work/wgcna/2021-02-15/Beta.cell.RNA/power-80-limma-TRUE/module-assignment.rds",
#   "alpha": "/lab/work/vivekrai/2020-01_vanderbilt_rna/work/wgcna/2021-02-15/Alpha.cell.RNA/power-80-limma-TRUE/module-assignment.rds"
# }

CONNS = dict()
for k in ["beta", "alpha"]:
    for cutoff in CUTOFFS:
        CONNS[f"{k}-{cutoff}"] = f"/lab/work/vivekrai/2019_sciATAC_new/work/2019-07-26_Cicero-2/{k}/{k}_gt-0.05.filtered.txt"

MODULE_RES = dict()
for k in ["Beta.cell.RNA", "Alpha.cell.RNA"]:
    for cutoff in CUTOFFS:
        dir_label = k.replace(".cell.RNA", "").lower()
        MODULE_RES[f"{dir_label}-{cutoff}"] = f"/lab/work/vivekrai/2020-01_vanderbilt_rna/work/wgcna-explore/2021-04-07/{k}/power-80_cut-{cutoff}/module-assignment.rds"


annotations = glob_wildcards("data/{id}.bed")

CICERO_SCRIPT = "bin/annotate_cicero_by_bed.py"
MODULE_PEAK_SCRIPT = "bin/get_module_linked_peaks.R"
CREATE_ANNOT_LINK_SCRIPT = "bin/create_annot_link_file.R"


rule all:
    input:
        expand(
            work("{type}/{annot}/peaks/{annot}.txt"),
            type=CONNS.keys(),
            annot=annotations.id,
        ),
        expand(
            work("{type}/{annot}/peaks/annotation_link_file.txt"),
            type=CONNS.keys(),
            annot=annotations.id,
        ),


rule annotate_peaks:
    input:
        cicero_conn = lambda wildcards: CONNS[wildcards.type],
        annotation = "data/{annot}.bed"
    output:
        work("{type}", "{annot}", "peaks", "{annot}.txt")
    shell:
        """
        python3 {CICERO_SCRIPT} {input.cicero_conn} {input.annotation} $HG19_SIZES > {output}
        """

rule get_module_linked_peaks:
    input:
        cicero_annot = work("{type}/{annot}/peaks/{annot}.txt"),
        module_res = lambda wildcards: MODULE_RES[wildcards.type]
    output:
        annot_link_file = work("{type}/{annot}/peaks/annotation_link_file.txt")
    params:
        module_peaks_dir = directory(os.path.abspath(work("{type}/{annot}/peaks")))
    shell:
        """
          Rscript {MODULE_PEAK_SCRIPT} --cicero {input.cicero_annot} --celltype-label {wildcards.type} \
          --module-res {input.module_res} --out-dir-prefix {params.module_peaks_dir}

          ln -s /home/vivekrai/analyses/2019_sciATAC_new/control/2019-10-07_forReviewers/2019-10-15_compareBulkPeaks/compare_peaks/Alpha.bed {params.module_peaks_dir}/Alpha_sciATAC.bed
          ln -s /home/vivekrai/analyses/2019_sciATAC_new/control/2019-10-07_forReviewers/2019-10-15_compareBulkPeaks/compare_peaks/Alpha_Arda.bed {params.module_peaks_dir}/Alpha_Arda.bed
          ln -s /home/vivekrai/analyses/2019_sciATAC_new/control/2019-10-07_forReviewers/2019-10-15_compareBulkPeaks/compare_peaks/Alpha_Ackermann.bed {params.module_peaks_dir}/Alpha_Ackermann.bed

          ln -s /home/vivekrai/analyses/2019_sciATAC_new/control/2019-10-07_forReviewers/2019-10-15_compareBulkPeaks/compare_peaks/Beta.bed {params.module_peaks_dir}/Beta_sciATAC.bed
          ln -s /home/vivekrai/analyses/2019_sciATAC_new/control/2019-10-07_forReviewers/2019-10-15_compareBulkPeaks/compare_peaks/Beta_Arda.bed {params.module_peaks_dir}/Beta_Arda.bed
          ln -s /home/vivekrai/analyses/2019_sciATAC_new/control/2019-10-07_forReviewers/2019-10-15_compareBulkPeaks/compare_peaks/Beta_Ackermann.bed {params.module_peaks_dir}/Beta_Ackermann.bed

          ln -s /home/vivekrai/analyses/2019_sciATAC_new/control/2019-10-07_forReviewers/2019-10-15_compareBulkPeaks/compare_peaks/sci-ATAC-seq.broadPeak.noblacklist {params.module_peaks_dir}/Bulk_sci-ATAC-seq.bed

          Rscript {CREATE_ANNOT_LINK_SCRIPT} --input {params.module_peaks_dir} \
            --outfile {output.annot_link_file} --celltype {wildcards.type}
        """
