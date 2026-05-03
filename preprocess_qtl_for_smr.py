#!/usr/bin/env python3
import argparse
import re
import subprocess
from pathlib import Path
import polars as pl
from liftover import get_lifter

def load_gene_positions_from_gtf(gtf_path: str) -> pl.DataFrame:
    gtf = pl.read_csv(
        gtf_path,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        new_columns=[
            "chrom", "source", "feature", "start", "end",
            "score", "strand", "frame", "attributes"
        ],
        truncate_ragged_lines=True
    )

    gene_id_pattern = r'gene_id "([^"]+)"'
    gene_name_pattern = r'gene_name "([^"]+)"'

    genes = (
        gtf
        .filter(pl.col("feature") == "gene")
        .with_columns([
            pl.col("attributes").map_elements(
                lambda x: re.search(gene_id_pattern, x).group(1) if re.search(gene_id_pattern, x) else None,
                return_dtype=pl.Utf8
            ).alias("gtf_gene_id"),
            pl.col("attributes").map_elements(
                lambda x: re.search(gene_name_pattern, x).group(1) if re.search(gene_name_pattern, x) else None,
                return_dtype=pl.Utf8
            ).alias("gtf_gene_symbol")
        ])
        .select(["gtf_gene_symbol", "gtf_gene_id", "chrom", "start", "end", "strand"])
        .unique(subset=["gtf_gene_symbol"])
    )

    return genes


def liftover_dataframe(
    df: pl.DataFrame,
    chr_col: str,
    bp_col: str,
    src_build: str,
    dst_build: str
) -> pl.DataFrame:

    lifter = get_lifter(src_build, dst_build)

    def liftover_row(row):
        chrom = row[chr_col]
        pos = row[bp_col]

        if chrom is None or pos is None:
            return {"CHR_harmonised": None, "BP_harmonised": None}

        chrom = str(chrom).strip()
        pos = str(pos).strip()

        if chrom == "" or pos == "":
            return {"CHR_harmonised": None, "BP_harmonised": None}

        chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)

        try:
            pos = int(float(pos))
        except ValueError:
            return {"CHR_harmonised": None, "BP_harmonised": None}

        try:
            hits = lifter[chrom][pos]
        except KeyError:
            try:
                hits = lifter[f"chr{chrom}"][pos]
            except Exception:
                return {"CHR_harmonised": None, "BP_harmonised": None}

        if not hits:
            return {"CHR_harmonised": None, "BP_harmonised": None}

        hit_chr = hits[0][0]
        hit_pos = hits[0][1]
        hit_chr = re.sub(r"^chr", "", str(hit_chr), flags=re.IGNORECASE)

        return {"CHR_harmonised": hit_chr, "BP_harmonised": int(hit_pos)}

    lifted = (
        df
        .with_columns(
            pl.struct([chr_col, bp_col]).map_elements(
                liftover_row,
                return_dtype=pl.Struct([
                    pl.Field("CHR_harmonised", pl.Utf8),
                    pl.Field("BP_harmonised", pl.Int64)
                ])
            ).alias("lifted")
        )
        .unnest("lifted")
    )

    return lifted


def prepare_qtl_with_gene_positions(
    input_path: str,
    gtf_path: str,
    gene_col: str,
    snp_col: str,
    chr_col: str,
    bp_col: str,
    a1_col: str,
    a2_col: str,
    freq_col: str,
    beta_col: str,
    se_col: str,
    p_col: str,
    input_build: str,
    output_build: str,
    default_freq: float | None
) -> pl.DataFrame:

    print("Reading QTL input...", flush=True)

    schema_overrides = {
        chr_col: pl.Utf8,
        bp_col: pl.Utf8,
        gene_col: pl.Utf8,
        snp_col: pl.Utf8,
        a1_col: pl.Utf8,
        a2_col: pl.Utf8,
        beta_col: pl.Utf8,
        se_col: pl.Utf8,
        p_col: pl.Utf8,
    }

    if freq_col != "NA":
        schema_overrides[freq_col] = pl.Utf8

    qtl = pl.read_csv(
        input_path,
        separator="\t",
        has_header=True,
        infer_schema_length=10000,
        truncate_ragged_lines=True,
        schema_overrides=schema_overrides
    )

    required_cols = [
        gene_col, snp_col, chr_col, bp_col,
        a1_col, a2_col, beta_col, se_col, p_col
    ]

    if freq_col != "NA":
        required_cols.append(freq_col)

    missing = []

    for col in required_cols:
        if col not in qtl.columns:
            missing.append(col)

    if len(missing) > 0:
        raise ValueError(f"Missing required columns in input: {missing}")

    if input_build != output_build:
        print(f"Lifting over coordinates: {input_build} -> {output_build}...", flush=True)

        qtl = liftover_dataframe(
            df=qtl,
            chr_col=chr_col,
            bp_col=bp_col,
            src_build=input_build,
            dst_build=output_build
        )
    else:
        qtl = qtl.with_columns([
            pl.col(chr_col).cast(pl.Utf8).alias("CHR_harmonised"),
            pl.col(bp_col).cast(pl.Utf8).str.strip_chars().cast(pl.Float64, strict=False).cast(pl.Int64, strict=False).alias("BP_harmonised")
        ])

    if freq_col == "NA":
        if default_freq is None:
            raise ValueError("freq_col is NA, but --default-freq was not supplied.")

        qtl = qtl.with_columns([
            pl.lit(float(default_freq)).alias("Freq_input")
        ])

        freq_expr = pl.col("Freq_input").cast(pl.Float64, strict=False)

    else:
        freq_expr = pl.col(freq_col).cast(pl.Utf8).str.strip_chars().cast(pl.Float64, strict=False)

    print("Loading gene positions from GTF...", flush=True)
    gene_positions = load_gene_positions_from_gtf(gtf_path)

    print("Joining QTLs to gene positions...", flush=True)

    qtl = (
        qtl
        .with_columns([
            pl.col(gene_col).cast(pl.Utf8).str.strip_chars().alias("gene_symbol_join")
        ])
        .join(
            gene_positions,
            left_on="gene_symbol_join",
            right_on="gtf_gene_symbol",
            how="left"
        )
        .select([
            pl.col("CHR_harmonised").cast(pl.Utf8).str.replace("^chr", "").alias("Chr"),
            pl.col(snp_col).cast(pl.Utf8).str.strip_chars().alias("SNP"),
            pl.col("BP_harmonised").cast(pl.Int64).alias("Bp"),
            pl.col(a1_col).cast(pl.Utf8).str.strip_chars().alias("A1"),
            pl.col(a2_col).cast(pl.Utf8).str.strip_chars().alias("A2"),
            freq_expr.alias("Freq"),
            pl.col(beta_col).cast(pl.Utf8).str.strip_chars().cast(pl.Float64, strict=False).alias("Beta"),
            pl.col(se_col).cast(pl.Utf8).str.strip_chars().cast(pl.Float64, strict=False).alias("se"),
            pl.col(p_col).cast(pl.Utf8).str.strip_chars().cast(pl.Float64, strict=False).alias("p"),
            pl.col(gene_col).cast(pl.Utf8).str.strip_chars().alias("Probe"),
            pl.col("gtf_gene_id").cast(pl.Utf8).alias("ProbeID"),
            pl.col("chrom").cast(pl.Utf8).str.replace("^chr", "").alias("ProbeChr"),
            pl.col("start").cast(pl.Int64).alias("ProbeBp"),
            pl.col("strand").cast(pl.Utf8).alias("Orientation"),
        ])
        .drop_nulls([
            "Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p",
            "Probe", "ProbeChr", "ProbeBp"
        ])
        .filter(
            (pl.col("SNP") != "") &
            (pl.col("A1") != "") &
            (pl.col("A2") != "") &
            (pl.col("Probe") != "")
        )
        .filter(
            pl.col("Chr").is_in([str(i) for i in range(1, 23)] + ["X", "Y"])
        )
    )

    print(f"Harmonised rows kept: {qtl.height}", flush=True)

    return qtl


def write_smr_input_files(qtl: pl.DataFrame, out_dir: str, output_label: str) -> Path:
    out_dir = Path(out_dir)
    esd_dir = out_dir / f"{output_label}_esd"

    out_dir.mkdir(parents=True, exist_ok=True)
    esd_dir.mkdir(parents=True, exist_ok=True)

    print("Writing probe-level ESD files...", flush=True)

    flist_rows = []

    for probe_df in qtl.partition_by("Probe", as_dict=False):
        probe = probe_df["Probe"][0]
        probe_id = probe_df["ProbeID"][0]
        probe_chr = probe_df["ProbeChr"][0]
        probe_bp = probe_df["ProbeBp"][0]
        orientation = probe_df["Orientation"][0]

        safe_probe = re.sub(r"[^A-Za-z0-9._-]", "_", str(probe))
        esd_path = esd_dir / f"{safe_probe}.esd"

        esd_df = (
            probe_df
            .select(["Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p"])
            .drop_nulls()
            .filter(
                (pl.col("SNP") != "") &
                (pl.col("A1") != "") &
                (pl.col("A2") != "") &
                (pl.col("A1") != pl.col("A2")) &
                pl.col("Freq").is_finite() &
                pl.col("Beta").is_finite() &
                pl.col("se").is_finite() &
                pl.col("p").is_finite() &
                (pl.col("Freq") >= 0) &
                (pl.col("Freq") <= 1) &
                (pl.col("se") > 0) &
                (pl.col("p") > 0) &
                (pl.col("p") <= 1)
            )
            .unique(subset=["SNP", "Bp", "A1", "A2"])
            .sort(["Chr", "Bp", "SNP"])
        )

        if esd_df.height == 0:
            continue

        esd_pd = esd_df.to_pandas()

        esd_pd.to_csv(
            esd_path,
            sep="\t",
            index=False,
            lineterminator="\n",
            encoding="utf-8"
        )

        with open(esd_path, "rb") as f:
            raw = f.read()

        if raw.startswith(b"\xef\xbb\xbf"):
            raw = raw[3:]

        raw = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")

        with open(esd_path, "wb") as f:
            f.write(raw)

        flist_rows.append({
            "Chr": probe_chr,
            "ProbeID": probe_id if probe_id is not None else probe,
            "GeneticDistance": 0,
            "ProbeBp": probe_bp,
            "Gene": probe,
            "Orientation": orientation if orientation is not None else "NA",
            "PathOfEsd": str(esd_path.resolve())
        })

    if not flist_rows:
        raise ValueError("No valid probes remained after harmonisation and filtering.")

    flist = pl.DataFrame(flist_rows).sort(["Chr", "ProbeBp", "Gene"])
    flist_path = out_dir / f"{output_label}.flist"
    flist.write_csv(flist_path, separator="\t")

    print(f"Probe ESD files written: {len(flist_rows)}", flush=True)

    return flist_path


def build_output_label(qtl_type: str, dataset_id: str, cell_type: str) -> str:
    if qtl_type == "sc":
        if cell_type == "bulk":
            raise ValueError("--cell-type cannot be 'bulk' when --qtl-type sc")
        return f"{dataset_id}_{cell_type}"

    if qtl_type == "bulk":
        return dataset_id

    raise ValueError("--qtl-type must be 'sc' or 'bulk'")


def main():
    parser = argparse.ArgumentParser(description="Prepare bulk or single-cell QTL summary statistics for SMR BESD generation.")
    parser.add_argument("--input", required=True, help="Input QTL TSV file")
    parser.add_argument("--gtf", required=True, help="GENCODE GTF file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--qtl-type", required=True, choices=["sc", "bulk"], help="QTL type: sc or bulk")
    parser.add_argument("--dataset-id", required=True, help="Dataset ID, e.g. fujita_eqtl or brainmeta_sqtl")
    parser.add_argument("--cell-type", default="bulk", help="Cell type for sc QTLs; use bulk for bulk QTLs")
    parser.add_argument("--gene-col", required=True)
    parser.add_argument("--snp-col", required=True)
    parser.add_argument("--chr-col", required=True)
    parser.add_argument("--bp-col", required=True)
    parser.add_argument("--a1-col", required=True)
    parser.add_argument("--a2-col", required=True)
    parser.add_argument("--freq-col", required=True)
    parser.add_argument("--beta-col", required=True)
    parser.add_argument("--se-col", required=True)
    parser.add_argument("--p-col", required=True)
    parser.add_argument("--default-freq", type=float, default=None, help="Fallback frequency if --freq-col NA")
    parser.add_argument(
        "--input-build",
        required=True,
        choices=["hg19", "hg38"],
        help="Genome build of input coordinates"
    )

    parser.add_argument(
        "--output-build",
        default="hg19",
        choices=["hg19", "hg38"],
        help="Target genome build for output coordinates"
    )

    parser.add_argument(
        "--run-smr",
        action="store_true",
        help="Run SMR --make-besd after writing flist/esd files"
    )

    parser.add_argument(
        "--smr-bin",
        default="smr",
        help="Path to SMR binary"
    )

    args = parser.parse_args()

    output_label = build_output_label(
        qtl_type=args.qtl_type,
        dataset_id=args.dataset_id,
        cell_type=args.cell_type
    )

    qtl = prepare_qtl_with_gene_positions(
        input_path=args.input,
        gtf_path=args.gtf,
        gene_col=args.gene_col,
        snp_col=args.snp_col,
        chr_col=args.chr_col,
        bp_col=args.bp_col,
        a1_col=args.a1_col,
        a2_col=args.a2_col,
        freq_col=args.freq_col,
        beta_col=args.beta_col,
        se_col=args.se_col,
        p_col=args.p_col,
        input_build=args.input_build,
        output_build=args.output_build,
        default_freq=args.default_freq
    )

    flist_path = write_smr_input_files(
        qtl=qtl,
        out_dir=args.outdir,
        output_label=output_label
    )

    if args.run_smr:
        print("Running SMR --make-besd...", flush=True)

        subprocess.run([
            args.smr_bin,
            "--eqtl-flist", str(flist_path),
            "--make-besd",
            "--out", str(Path(args.outdir) / output_label)
        ], check=True)

    print(f"QTL type: {args.qtl_type}", flush=True)
    print(f"Dataset ID: {args.dataset_id}", flush=True)
    print(f"Cell type: {args.cell_type}", flush=True)
    print(f"Output label: {output_label}", flush=True)
    print(f"flist: {flist_path}", flush=True)
    print(qtl.head(), flush=True)
    print(qtl.columns, flush=True)


if __name__ == "__main__":
    main()
