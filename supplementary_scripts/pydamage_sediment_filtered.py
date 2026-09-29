#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
from pathlib import Path
import sys
import os


def run_cmd(cmd, dry_run=False):
    print("[CMD]", " ".join(cmd))
    if dry_run:
        return
    subprocess.run(cmd, check=True)


def find_pydamage_result(outdir):
    """
    自动寻找 pyDamage 输出结果表。
    不同版本 pyDamage 输出文件名可能略有差异，因此这里递归寻找 csv/tsv。
    """
    outdir = Path(outdir)
    candidates = list(outdir.rglob("*.csv")) + list(outdir.rglob("*.tsv"))

    if not candidates:
        raise FileNotFoundError(
            f"没有在 {outdir} 下面找到 pyDamage 输出的 csv/tsv 文件。"
        )

    # 优先选择文件名中包含 result 的文件
    result_like = [
        x for x in candidates
        if "result" in x.name.lower() or "pydamage" in x.name.lower()
    ]

    if result_like:
        return result_like[0]

    return candidates[0]


def read_table_auto(path):
    path = Path(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    else:
        return pd.read_csv(path, sep="\t")


def find_column(df, candidates, required=True):
    """
    在 pyDamage 结果表中自动寻找列名。
    """
    lower_map = {c.lower(): c for c in df.columns}

    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]

    # 模糊匹配
    for col in df.columns:
        col_l = col.lower()
        for cand in candidates:
            if cand.lower() in col_l:
                return col

    if required:
        raise ValueError(
            f"无法在结果表中找到列：{candidates}\n"
            f"当前结果表包含列：{list(df.columns)}"
        )
    return None


def parse_pydamage_results(
    result_file,
    qvalue_threshold=0.05,
    damage_threshold=0.01,
    accuracy_threshold=0.8,
    min_reads=None,
    contig_col=None,
    qvalue_col=None,
    damage_col=None,
    accuracy_col=None,
    reads_col=None
):
    df = read_table_auto(result_file)

    print("\n[INFO] pyDamage 结果表：", result_file)
    print("[INFO] 结果表列名：")
    for c in df.columns:
        print("  -", c)

    # 自动识别常见列名
    if contig_col is None:
        contig_col = find_column(
            df,
            ["contig", "reference", "reference_id", "seqid", "sequence", "name"]
        )

    if qvalue_col is None:
        qvalue_col = find_column(
            df,
            ["qvalue", "q_value", "q-value", "q", "qval", "padj", "adjusted_pvalue"],
            required=False
        )

    if damage_col is None:
        damage_col = find_column(
            df,
            [
                "damage",
                "damage_rate",
                "damage_level",
                "damage_estimate",
                "d_max",
                "D_max",
                "lambda",
                "damage_model_p"
            ],
            required=False
        )

    if accuracy_col is None:
        accuracy_col = find_column(
            df,
            ["predicted_accuracy", "accuracy", "probability", "posterior_probability"],
            required=False
        )

    if reads_col is None:
        reads_col = find_column(
            df,
            ["reads", "n_reads", "num_reads", "coverage", "aligned_reads"],
            required=False
        )

    print("\n[INFO] 自动识别到的关键列：")
    print("  contig_col  =", contig_col)
    print("  qvalue_col  =", qvalue_col)
    print("  damage_col  =", damage_col)
    print("  accuracy_col=", accuracy_col)
    print("  reads_col   =", reads_col)

    filtered = df.copy()

    # 过滤 q-value / FDR
    if qvalue_col is not None:
        filtered = filtered[pd.to_numeric(filtered[qvalue_col], errors="coerce") <= qvalue_threshold]
    else:
        print("[WARN] 没有找到 q-value 列，跳过 q-value 过滤。")

    # 过滤损伤强度
    if damage_col is not None:
        filtered = filtered[pd.to_numeric(filtered[damage_col], errors="coerce") >= damage_threshold]
    else:
        print("[WARN] 没有找到 damage 列，跳过 damage 过滤。")

    # 过滤 predicted accuracy / probability
    if accuracy_col is not None:
        filtered = filtered[pd.to_numeric(filtered[accuracy_col], errors="coerce") >= accuracy_threshold]
    else:
        print("[WARN] 没有找到 accuracy/probability 列，跳过 accuracy 过滤。")

    # 过滤 read 数或 coverage
    if min_reads is not None and reads_col is not None:
        filtered = filtered[pd.to_numeric(filtered[reads_col], errors="coerce") >= min_reads]
    elif min_reads is not None and reads_col is None:
        print("[WARN] 指定了 min_reads，但没有找到 reads/coverage 列，跳过该过滤。")

    # 加一列标记
    df["ancient_candidate"] = df[contig_col].isin(filtered[contig_col])

    return df, filtered, contig_col


def extract_reads_by_contigs(bam, selected_contigs, out_bam, tmp_bed):
    """
    根据高可信 damaged contig 提取 reads。
    用 samtools view -L BED 实现。
    """
    if len(selected_contigs) == 0:
        print("[WARN] 没有高可信 contig，跳过 reads 提取。")
        return

    # 从 BAM header 中读取 contig 长度
    header_cmd = ["samtools", "view", "-H", bam]
    header = subprocess.check_output(header_cmd, text=True)

    contig_lengths = {}
    for line in header.splitlines():
        if line.startswith("@SQ"):
            fields = line.split("\t")
            sn = None
            ln = None
            for f in fields:
                if f.startswith("SN:"):
                    sn = f.replace("SN:", "")
                elif f.startswith("LN:"):
                    ln = int(f.replace("LN:", ""))
            if sn is not None and ln is not None:
                contig_lengths[sn] = ln

    selected_contigs = set(selected_contigs)

    with open(tmp_bed, "w") as f:
        for ctg in selected_contigs:
            if ctg in contig_lengths:
                f.write(f"{ctg}\t0\t{contig_lengths[ctg]}\n")
            else:
                print(f"[WARN] contig 不在 BAM header 中，跳过：{ctg}")

    run_cmd(["samtools", "view", "-b", "-L", tmp_bed, "-o", out_bam, bam])
    run_cmd(["samtools", "index", out_bam])


def summarize_by_taxonomy(filtered_df, taxonomy_file, contig_col, out_file):
    tax = pd.read_csv(taxonomy_file, sep="\t", header=None)

    if tax.shape[1] < 2:
        raise ValueError("taxonomy 文件至少需要两列：contig_id 和 taxon")

    tax = tax.iloc[:, :2]
    tax.columns = [contig_col, "taxon"]

    merged = filtered_df.merge(tax, on=contig_col, how="left")
    merged["taxon"] = merged["taxon"].fillna("Unknown")

    summary = (
        merged.groupby("taxon")
        .size()
        .reset_index(name="ancient_candidate_contig_count")
        .sort_values("ancient_candidate_contig_count", ascending=False)
    )

    summary.to_csv(out_file, sep="\t", index=False)
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Run pyDamage on sediment DNA BAM and filter ancient DNA candidate contigs."
    )

    parser.add_argument(
        "--bam",
        required=True,
        help="输入排序并建立索引的 BAM 文件，例如 sample.sorted.bam"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="输出目录"
    )

    parser.add_argument(
        "--run-pydamage",
        action="store_true",
        help="是否运行 pyDamage。如果已经运行过，可以不加该参数，脚本会直接解析结果。"
    )

    parser.add_argument(
        "--pydamage-result",
        default=None,
        help="已有 pyDamage 结果表。如果不提供，则自动在 outdir 中寻找。"
    )

    parser.add_argument(
        "--pydamage-cmd-template",
        default="pydamage analyze {bam} -o {outdir}",
        help=(
            "pyDamage 命令模板。不同版本 pyDamage 参数可能不同，"
            "如有需要可自行修改。默认：'pydamage analyze {bam} -o {outdir}'"
        )
    )

    parser.add_argument(
        "--qvalue",
        type=float,
        default=0.05,
        help="q-value/FDR 阈值，默认 0.05"
    )

    parser.add_argument(
        "--damage",
        type=float,
        default=0.01,
        help="损伤水平阈值，默认 0.01"
    )

    parser.add_argument(
        "--accuracy",
        type=float,
        default=0.8,
        help="predicted accuracy/probability 阈值，默认 0.8"
    )

    parser.add_argument(
        "--min-reads",
        type=float,
        default=None,
        help="最小 reads 数或 coverage 阈值，可选"
    )

    parser.add_argument(
        "--taxonomy",
        default=None,
        help="可选，contig 分类文件，两列：contig_id 和 taxon"
    )

    parser.add_argument(
        "--extract-reads",
        action="store_true",
        help="是否提取比对到高可信 ancient contig 的 reads"
    )

    parser.add_argument(
        "--contig-col",
        default=None,
        help="手动指定 pyDamage 结果表中的 contig 列名"
    )

    parser.add_argument(
        "--qvalue-col",
        default=None,
        help="手动指定 q-value 列名"
    )

    parser.add_argument(
        "--damage-col",
        default=None,
        help="手动指定 damage 列名"
    )

    parser.add_argument(
        "--accuracy-col",
        default=None,
        help="手动指定 accuracy/probability 列名"
    )

    parser.add_argument(
        "--reads-col",
        default=None,
        help="手动指定 reads/coverage 列名"
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pydamage_outdir = outdir / "pydamage"

    if args.run_pydamage:
        pydamage_outdir.mkdir(parents=True, exist_ok=True)

        cmd_string = args.pydamage_cmd_template.format(
            bam=args.bam,
            outdir=str(pydamage_outdir)
        )
        cmd = cmd_string.split()
        run_cmd(cmd)

    if args.pydamage_result is not None:
        result_file = Path(args.pydamage_result)
    else:
        result_file = find_pydamage_result(pydamage_outdir)

    all_df, filtered_df, contig_col = parse_pydamage_results(
        result_file=result_file,
        qvalue_threshold=args.qvalue,
        damage_threshold=args.damage,
        accuracy_threshold=args.accuracy,
        min_reads=args.min_reads,
        contig_col=args.contig_col,
        qvalue_col=args.qvalue_col,
        damage_col=args.damage_col,
        accuracy_col=args.accuracy_col,
        reads_col=args.reads_col
    )

    all_out = outdir / "all_pydamage_results_with_flag.tsv"
    high_out = outdir / "high_confidence_ancient_contigs.tsv"
    contig_list_out = outdir / "high_confidence_ancient_contigs.list"

    all_df.to_csv(all_out, sep="\t", index=False)
    filtered_df.to_csv(high_out, sep="\t", index=False)

    selected_contigs = filtered_df[contig_col].astype(str).tolist()

    with open(contig_list_out, "w") as f:
        for c in selected_contigs:
            f.write(c + "\n")

    print("\n[RESULT] 总 contig 数：", all_df.shape[0])
    print("[RESULT] 高可信 ancient candidate contig 数：", len(selected_contigs))
    print("[OUTPUT]", all_out)
    print("[OUTPUT]", high_out)
    print("[OUTPUT]", contig_list_out)

    if args.taxonomy is not None:
        tax_out = outdir / "ancient_signal_by_taxon.tsv"
        summary = summarize_by_taxonomy(
            filtered_df=filtered_df,
            taxonomy_file=args.taxonomy,
            contig_col=contig_col,
            out_file=tax_out
        )
        print("[OUTPUT]", tax_out)
        print("\n[INFO] 分类群古 DNA 信号概览：")
        print(summary.head(20).to_string(index=False))

    if args.extract_reads:
        out_bam = outdir / "reads_mapped_to_high_confidence_ancient_contigs.bam"
        tmp_bed = outdir / "high_confidence_ancient_contigs.bed"
        extract_reads_by_contigs(
            bam=args.bam,
            selected_contigs=selected_contigs,
            out_bam=str(out_bam),
            tmp_bed=str(tmp_bed)
        )
        print("[OUTPUT]", out_bam)


if __name__ == "__main__":
    main()
