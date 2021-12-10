#! /usr/bin/env python3

import os
import pandas
import glob

from functools import partial
from datetime import date
from collections import namedtuple

# run_date = str(date.today())
run_date = "2021-06-09"
work = partial(os.path.join, "work", run_date)

# GARFIELD data files; Rarely changed.
# ------------------------------------
GARFIELD_DATA = config["GARFIELD"]["DATA"]
GARFIELD_SCRIPTS = config["GARFIELD"]["SCRIPTS"]
PRUNETAGSDIR = os.path.join(GARFIELD_DATA, "tags/r01")
CLUMPTAGSDIR = os.path.join(GARFIELD_DATA, "tags/r08")
MAFTSSDDIR = os.path.join(GARFIELD_DATA, "maftssd")
PVALDIR = os.path.join(GARFIELD_DATA, "pval")

# Change these values based on your analysis
# ------------------------------------------
PARAMETERS = config["PARAMETERS"]

PGWAS = PARAMETERS["pthresh"].split(",")
TRAITS = PARAMETERS["traits"]
CHROM = range(1, 23)

assert len(TRAITS) > 0

_annot_tuple = namedtuple("Annotation", ["celltype", "tss_def", "link_file"])

ANNOT_GROUPS_TUPLE = []
for x in glob.glob(f"work/{run_date}/*/hg19.tss.10kb_ext/peaks/annotation_link_file.txt"):
    _path = x.split("/")
    ANNOT_GROUPS_TUPLE.append(_annot_tuple(_path[2], _path[3], x))


def get_annot_link_file(wildcards):
    for val in ANNOT_GROUPS_TUPLE:
        if val.celltype == wildcards.celltype and val.tss_def == wildcards.tss_def:
            return val.link_file

def get_annotations_list(wildcards):
    df = pandas.read_csv(get_annot_link_file(wildcards), sep="\t")
    annots = df["Annotation"].tolist()
    assert len(annots) > 0
    return annots

def get_annot_bed_file(wildcards):
    df = pandas.read_csv(get_annot_link_file(wildcards), sep="\t")
    return df[df["Annotation"] == wildcards.annotation].iloc[0]["Path"]

def get_padj(filename):
    with open(filename) as m:
        padj = float(m.readlines()[1].split("\t")[1].rstrip())
    return padj


# BEGIN Pipeline
# --------------

wildcard_constraints:
    chrom="\d+"

rule all:
    """Single Annotation Workflow, includes function definitions from included Snakefile """
    input:
        # Prepare trait GWAS files
        gwas_files=expand(
            os.path.join(PVALDIR, "{trait}", "chr{chrom}"),
            trait=TRAITS,
            chrom=CHROM
        ),
        ## Calculations below are done "PER ANNOTATION"
        annotation_files=expand(
            work("{group.celltype}", "{group.tss_def}", "annotations", "chr{chrom}"),
            group=ANNOT_GROUPS_TUPLE,
            chrom=CHROM
        ),
        enrichment=expand(
            work("{group.celltype}", "{group.tss_def}", "results", "garfield.enrichment.{trait}.out"),
            group=ANNOT_GROUPS_TUPLE,
            trait=TRAITS
        ),
        meff=expand(
            work("{group.celltype}", "{group.tss_def}", "results", "garfield.Meff.{trait}.out"),
            group=ANNOT_GROUPS_TUPLE,
            trait=TRAITS
        )
        # significant_annots=expand(work("{celltype}", "{tss_def}",
        #     "results", "garfield.significant_annots.pgwas{pgwas}.padj.{trait}.out"
        #     ),

        # variants=work("{celltype}", "{tss_def}",
        #     "results", "garfield.variants.pgwas{pgwas}.padj.{trait}.out"
        # )

rule prepare_annotations:
    input:
        snp=os.path.join(GARFIELD_DATA, "annotation", "chr{chrom}"), # SNP reference
        annot=lambda wildcards: get_annot_bed_file(wildcards)
    output:
        work("{celltype}", "{tss_def}", "annotations", "intersects", "chr{chrom}", "{annotation}.bed")
    shell:
        """
        awk '{{print "chr{wildcards.chrom}\t"$1-1,$1}}' OFS='\t' {input.snp} | \
            intersectBed -a stdin -b {input.annot} -c | \
                awk '{{ if (($4>1)) {{$4=1}}; print $3,$4 }}' OFS='\t' > {output}
        """

rule collapseIntersectResults:
    """Assemble all input into one file. Don't include TSS distance column for now"""
    input:
        snp=os.path.join(GARFIELD_DATA, "annotation", "chr{chrom}"),
        annot_link_file=lambda wildcards: get_annot_link_file(wildcards),
        intersect_results=lambda wildcards: expand(
            work("{celltype}", "{tss_def}", "annotations", "intersects", "chr{chrom}", "{annotation}.bed"),
            annotation=get_annotations_list(wildcards),
            allow_missing=True
        ),
    output:
        work("{celltype}", "{tss_def}", "annotations", "chr{chrom}")
    params:
        script="bin/prep_annotation_overlap.py"
    shell:
        """
        python {params.script} --link {input.annot_link_file} \
            --annots {input.intersect_results} --snp {input.snp} \
            --output {output}
        """

rule prune_and_clump:
    """Prune SNPs r2 0.1 and clump. Runs ./garfield-prep-chr which prepares
    data for garfield-test.Randgarfield-Meff-Padj.R."""
    input:
        # GARFIELD related data; includes trait information as well
        garfield_prep_chr=os.path.join(GARFIELD_SCRIPTS, "garfield-prep-chr"),
        prune=os.path.join(PRUNETAGSDIR, "chr{chrom}"),
        clump=os.path.join(CLUMPTAGSDIR, "chr{chrom}"),
        maf=os.path.join(MAFTSSDDIR, "chr{chrom}"),
        pval=os.path.join(PVALDIR, "{trait}", "chr{chrom}"),
        # Output from previous steps
        annot=work("{celltype}", "{tss_def}", "annotations", "chr{chrom}"),
    output:
        work("{celltype}", "{tss_def}", "results", "{trait}/prep/chr{chrom}")
    params:
        exclude=(
            f" -excl {PARAMETERS['exclude_annots']} "
            if PARAMETERS["exclude_annots"] != ""
            else ""
        ),
    shell:
        # ' for CHR in `seq 1 22`; do '
        "{input.garfield_prep_chr} "
        " -ptags {input.prune} "
        " -ctags {input.clump} "
        " -maftss {input.maf} "
        " -pval {input.pval} "
        " -ann {input.annot} "
        " {params.exclude} "
        " -chr {wildcards.chrom} "
        ' -o {output} || {{ echo "Failure!"; }} ;'


rule concat_prep:
    """Concat prepped data for chroms """
    input:
        inputs=lambda wildcards: expand(
            work("{celltype}", "{tss_def}", "results", "{trait}/prep/chr{chrom}"),
            chrom=CHROM,
            allow_missing=True
        ),
    output:
        work("{celltype}", "{tss_def}", "results", "garfield.prep.{trait}.out")
    shell:
        """
        cat {input.inputs} | sort -t' ' -k1,1 -n > {output}
        """


rule calculate_effective_annots:
    """Calculate effective number of annotations"""
    input:
        garfield_meff=os.path.join(GARFIELD_SCRIPTS, "garfield-Meff-Padj.R"),
        prep=rules.concat_prep.output,
    output:
        work("{celltype}", "{tss_def}", "results", "garfield.Meff.{trait}.out")
    shell:
        " Rscript {input.garfield_meff} "
        " -i {input.prep} "
        " -o {output} "


rule test_enrichment:
    """Calculate effective number of annotations"""
    input:
        garfield_test=os.path.join(GARFIELD_SCRIPTS, "garfield-test.R"),
        prep=rules.concat_prep.output,
        annot_link_file=lambda wildcards: get_annot_link_file(wildcards)
    output:
        work("{celltype}", "{tss_def}", "results", "garfield.enrichment.{trait}.out")
    shell:
        " Rscript {input.garfield_test} "
        " -i {input.prep} "
        " -o {output} "
        " -l {input.annot_link_file} "
        " -pt {PARAMETERS[pthresh]} "
        " -b {PARAMETERS[binning]} "
        " -s {PARAMETERS[subset_annots]} "
        " -c 0 "

rule extract_variants:
    input:
        extract=os.path.join(
            GARFIELD_SCRIPTS,
            "garfield_extract_variants_overlapping_enriched_annotations.sh",
        ),
        prep=rules.concat_prep.output,
        meff=rules.calculate_effective_annots.output,
        enrichment=rules.test_enrichment.output,
    output:
        significant_annots=work("{celltype}", "{tss_def}",
            "results", "garfield.significant_annots.pgwas{pgwas}.padj.{trait}.out"
        ),
        variants=work("{celltype}", "{tss_def}",
            "results", "garfield.variants.pgwas{pgwas}.padj.{trait}.out"
        )
    params:
        pvaldir=os.path.join(PVALDIR, "{trait}"),
    shell:
        """ padj=`less {input.meff} | awk '{{if ((NR==2)) print $2}}'` ;"""
        " bash {input.extract} "
        " {PRUNETAGSDIR} "
        " {CLUMPTAGSDIR} "
        " {ANNOTDIR} "
        " {params.pvaldir} "
        " {input.prep} "
        " {input.enrichment} "
        " {output.significant_annots} "
        " {output.variants} "
        " {wildcards.pgwas} "
        " $padj "


# rule create_plots:
#     """MAke default and custom GARFIELD plots"""
#     input:
#         custom_plot=SCRIPTS["plot"],
#         meff=rules.calculate_effective_annots.output.main,
#         enrichment=rules.test_enrichment.output.main,
#         annot_link_file=PARAMETERS["annot_link_file"],
#     output:
#         custom_plot=expand(
#             os.path.join(d_figures, "{{trait}}.enrichment-ci95.{pgwas}.{format}"),
#             pgwas=PGWAS,
#             format=["png", "pdf"],
#         ),
#     params:
#         script="scripts/plot_garfield.R",
#         prefix=os.path.join(d_figures, "{trait}.enrichment-ci95"),
#     shell:
#         """ padj=`less {input.meff} | awk '{{if ((NR==2)) print $2}}'` ;
#         Rscript {params.script} {input.enrichment} {params.prefix} $padj
#         """


# rule conditional_analysis:
#     """Prioritize relevant annotations by conditional analysis"""
#     input:
#         garfield_test=os.path.join(GARFIELD_SCRIPTS, "garfield-test.R"),
#         prep=rules.concat_prep.output.main,
#         annot_link_file=PARAMETERS["annot_link_file"],
#         meff=rules.calculate_effective_annots.output.main,
#         enrichment=rules.test_enrichment.output.main,
#     params:
#         subset_annots=(
#             " -s {PARAMETERS['subset_annots']}"
#             if PARAMETERS["subset_annots"] != ""
#             else ""
#         ),
#     output:
#         main=os.path.join(
#             d_garfield, "garfield.enrichment.{trait}.out.0.05.cond.indep"
#         ),
#     shell:
#         """ padj=`less {input.meff} | awk '{{if ((NR==2)) print $2}}'` ;"""
#         " Rscript {input.garfield_test} "
#         " -i {input.prep} "
#         " -o {input.enrichment} "
#         " -l {input.annot_link_file} "
#         " -pt {PARAMETERS[pthresh]} "
#         " -b {PARAMETERS[binning]} "
#         " -s {PARAMETERS[subset_annots]} "
#         " -c 1 "
#         " -padj $padj "
#         " -ct 0.05 "
