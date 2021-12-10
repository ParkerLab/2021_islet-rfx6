#!/usr/bin/env python3
# coding: utf-8

import pandas
import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import glob

dpi = 150
matplotlib.rcParams["figure.dpi"] = dpi
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams.update({"font.size": 12})

import pandas_extra
import argparse

# In[41]:
def getOpts():
    parser = argparse.ArgumentParser(description="Make RNA qc plot for channels")
    parser.add_argument(
        "--input", required=True, nargs="+", help=""".best outputs from demuxlet"""
    )
    parser.add_argument("--weights", help="""data with sample names and weights""")
    parser.add_argument("--batch", type=int, help="""Batch number""")
    parser.add_argument("--data", required=True, help="""Output data file name.""")
    parser.add_argument("--qc-file", required=True, help="""Path to Ataqv .txt files""")
    parser.add_argument(
        "--fig-dir", default=".", help="""Directory to output figures."""
    )
    args = parser.parse_args()
    return args


def assemble(f):
    t = pandas.read_csv(
        f,
        sep="\t",
        usecols=["BARCODE", "RD.UNIQ", "BEST", "SNG.1ST", "SNG.LLK1", "PRB.SNG1"],
    )
    t.rename(columns={"BARCODE": "barcode"}, inplace=True)
    t["out"] = t["BEST"].map(lambda x: x.split("-")[0])

    library = os.path.basename(f).split(".")[0]
    split = os.path.basename(f).split(".")[1]
    q = pandas.read_csv(
        args.qc_file,
        sep="\t",
        header=None,
        names=["library", "barcode", "metric", "value"],
    )
    metrics = [
        "hqaa",
        "tss_enrichment",
        "percent_mitochondrial",
        "total_autosomal_reads",
    ]
    q = q[q["metric"].isin(metrics)]
    q = q[q["value"] != "None"]
    q["value"] = q["value"].astype(float)
    q["batch"], q["library"] = zip(
        *q["library"].map(lambda x: getlib(x.replace(".ataqv.json.gz", "")))
    )

    d = q.pivot_table(
        index=["batch", "library", "barcode"], columns="metric", values="value"
    ).reset_index()
    d["fraction_mitochondrial"] = d["percent_mitochondrial"].map(lambda x: x / 100)

    t = pandas.merge(
        t,
        d[["barcode", "hqaa", "tss_enrichment", "fraction_mitochondrial"]],
        how="inner",
        on="barcode",
    )
    t.sort_values("out", inplace=True)
    t["library"] = library
    t["split"] = split
    return t


def getlib(x):
    s = x.split("-")
    batch = f"{s[0]}-{s[1]}"
    library = s[2]
    return [batch, library]


if __name__ == "__main__":
    args = getOpts()

    pe = pandas_extra.ExtraFunctions(args.fig_dir)

    # Read different info
    # info = pandas.read_csv(args.weights, sep='\t', usecols=['batch', 'sample','Tissue amount'])
    # info['sample'] = info['sample'].map(lambda x: x.replace('M', "") if x.startswith('M') else x)
    # info = info[info['batch']==args.batch]
    # info['Tissue amount'] = info['Tissue amount'].map(lambda x: float(str(x).replace(" mg", "")))

    d = pandas.concat([assemble(f) for f in args.input])

    d["log_hqaa"] = d["hqaa"].map(lambda x: math.log(x, 10))
    d["log_mito"] = d["fraction_mitochondrial"].map(
        lambda x: math.log((x + 0.0001), 10)
    )
    d["log_tss"] = d["tss_enrichment"].map(lambda x: math.log((x + 0.0001), 10))
    d.to_csv(args.data, sep="\t", na_rep="NA", index=False)
    # d['out'] = d['out'].str.replace('AMB', 'zAMB')
    # d.sort_values('out', inplace=True)
    # scatterplot colored by demuxlet outcome
    plt.figure(figsize=(5, 4))
    sns.scatterplot(
        data=d,
        x="log_hqaa",
        y="log_mito",
        hue="out",
        **{"alpha": 0.1},
        hue_order=["DBL", "SNG", "AMB"],
    )
    plt.xlabel("log10(HQAA)")
    plt.ylabel("log10(fraction mitochondrial)")
    pe.save("fig.demuxlet_qc_mito.png")

    plt.figure(figsize=(5, 4))
    sns.scatterplot(
        data=d,
        x="log_hqaa",
        y="log_tss",
        hue="out",
        **{"alpha": 0.1},
        hue_order=["DBL", "SNG", "AMB"],
    )
    plt.xlabel("log10(HQAA)")
    plt.ylabel("log10(TSS enrichment)")
    pe.save("fig.demuxlet_qc_tss.png")

    sub = d[d["out"] == "SNG"]
    df = pandas.DataFrame(sub.groupby("SNG.1ST").size()).reset_index()
    df.rename(columns={0: "N"}, inplace=True)
    df["SNG.1ST"] = df["SNG.1ST"].astype(str)

    # df = pandas.merge(df, info, how="inner", left_on="SNG.1ST", right_on="sample")
    # df.sort_values("N", inplace=True)

    # Compare nuclei counts with input sample weights

    # sns.lmplot(data=df, x="N", y="Tissue amount")
    # plt.xlabel("# Nuclei assigned as singlet ")
    # plt.ylabel("Tissue amount (mg)")
    # pe.save("fig.quality_nuclei_per_sample_by_weight.pdf")
