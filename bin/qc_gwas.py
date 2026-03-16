#!/usr/bin/env python3
import argparse
from pathlib import Path
import polars as pl

bases = ["A", "C", "G", "T"]

def count(label: str, df: pl.DataFrame) -> None:
    print(f"{label}: {df.height:,}")

def read_table(path: Path, sep: str, has_header: bool = True) -> pl.DataFrame:
    return pl.read_csv(
        path,
        separator=sep,
        has_header=has_header,
        infer_schema_length=5000,
        ignore_errors=True,
        null_values=["NA", "NaN", "nan", "", "null", "NULL"],
    )

def as_upper(df: pl.DataFrame, col: str) -> pl.DataFrame:
    if col not in df.columns:
        return df
    return df.with_columns(
        pl.col(col).cast(pl.Utf8, strict=False).str.strip_chars().str.to_uppercase().alias(col)
    )

def to_int(df: pl.DataFrame, col: str) -> pl.DataFrame:
    if col not in df.columns:
        return df
    return df.with_columns(pl.col(col).cast(pl.Int64, strict=False).alias(col))

def to_float(df: pl.DataFrame, col: str) -> pl.DataFrame:
    if col not in df.columns:
        return df
    return df.with_columns(pl.col(col).cast(pl.Float64, strict=False).alias(col))

def exclude_region(df: pl.DataFrame, chr_col: str, pos_col: str, chrom: int, start: int, end: int) -> pl.DataFrame:
    if chr_col not in df.columns or pos_col not in df.columns:
        return df
    bad = (pl.col(chr_col) == chrom) & (pl.col(pos_col) >= start) & (pl.col(pos_col) <= end)
    return df.filter(~bad)

def keep_snps_only(df: pl.DataFrame, a1_col: str, a2_col: str) -> pl.DataFrame:
    if a1_col not in df.columns or a2_col not in df.columns:
        return df
    a1 = pl.col(a1_col)
    a2 = pl.col(a2_col)
    ok_len = (a1.str.len_chars() == 1) & (a2.str.len_chars() == 1)
    ok_bases = a1.is_in(bases) & a2.is_in(bases)
    no_gap = ~a1.str.contains("-") & ~a2.str.contains("-")
    return df.filter(ok_len & ok_bases & no_gap)

def drop_palindromes(df: pl.DataFrame, a1_col: str, a2_col: str) -> pl.DataFrame:
    if a1_col not in df.columns or a2_col not in df.columns:
        return df
    a1 = pl.col(a1_col)
    a2 = pl.col(a2_col)
    pal = (
        ((a1 == "A") & (a2 == "T")) | ((a1 == "T") & (a2 == "A")) |
        ((a1 == "C") & (a2 == "G")) | ((a1 == "G") & (a2 == "C"))
    )
    return df.filter(~pal)

def add_frq_maf(
    df: pl.DataFrame,
    frq_col: str,
    eaf_col: str | None,
    fcas_col: str | None,
    fcon_col: str | None,
    ncas_col: str | None,
    ncon_col: str | None,
) -> pl.DataFrame:
    if frq_col in df.columns:
        pass
    elif eaf_col and (eaf_col in df.columns):
        df = df.with_columns(pl.col(eaf_col).cast(pl.Float64, strict=False).alias(frq_col))
    elif fcas_col and fcon_col and (fcas_col in df.columns) and (fcon_col in df.columns):
        fcas = pl.col(fcas_col)
        fcon = pl.col(fcon_col)
        if ncas_col and ncon_col and (ncas_col in df.columns) and (ncon_col in df.columns):
            ncas = pl.col(ncas_col)
            ncon = pl.col(ncon_col)
            w = ncas + ncon
            frq = pl.when(w > 0).then((fcas * ncas + fcon * ncon) / w).otherwise((fcas + fcon) / 2.0)
        else:
            frq = (fcas + fcon) / 2.0
        df = df.with_columns(frq.alias(frq_col))
    else:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias(frq_col))
    df = df.with_columns(
        pl.min_horizontal(pl.col(frq_col), (1.0 - pl.col(frq_col))).alias("MAF")
    )
    return df

def add_n(
    df: pl.DataFrame,
    n_out_col: str,
    n_col: str | None,
    ncas_col: str | None,
    ncon_col: str | None,
) -> pl.DataFrame:
    if n_out_col in df.columns:
        return df
    if n_col and (n_col in df.columns):
        return df.with_columns(pl.col(n_col).cast(pl.Float64, strict=False).alias(n_out_col))
    if ncas_col and ncon_col and (ncas_col in df.columns) and (ncon_col in df.columns):
        return df.with_columns(
            (pl.col(ncas_col).cast(pl.Float64, strict=False) + pl.col(ncon_col).cast(pl.Float64, strict=False)).alias(n_out_col)
        )
    return df.with_columns(pl.lit(None, dtype=pl.Float64).alias(n_out_col))

def add_info(df: pl.DataFrame, info_out_col: str, info_col: str | None, require_info: bool) -> pl.DataFrame:
    if info_out_col in df.columns:
        return df
    if info_col and (info_col in df.columns):
        return df.with_columns(pl.col(info_col).cast(pl.Float64, strict=False).alias(info_out_col))
    if require_info:
        raise SystemExit("ERROR: INFO required but --info_col not provided or not found")
    return df.with_columns(pl.lit(None, dtype=pl.Float64).alias(info_out_col))

def filter_maf_info(df: pl.DataFrame, maf_min: float, info_min: float, info_col: str, require_info: bool) -> pl.DataFrame:
    if "MAF" in df.columns:
        df = df.filter(pl.col("MAF") >= maf_min)
    if require_info:
        df = df.filter(pl.col(info_col).is_not_null() & (pl.col(info_col) >= info_min))
    else:
        if info_col in df.columns:
            df = df.filter(pl.col(info_col).is_null() | (pl.col(info_col) >= info_min))
    return df

def drop_missing_required(df: pl.DataFrame, required_cols: list[str]) -> pl.DataFrame:
    have = [c for c in required_cols if c in df.columns]
    return df.drop_nulls(subset=have)

def make_ldsc_and_write(
    df: pl.DataFrame,
    out: Path,
    snp_col: str,
    a1_col: str,
    a2_col: str,
    frq_col: str,
    n_col_out: str,
    beta_col: str,
    se_col: str,
    p_col: str,
    chr_col: str,
    pos_col: str,
    info_col_out: str,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    def safe(colname: str) -> pl.Expr:
        if colname in df.columns:
            return pl.col(colname)
        return pl.lit(None)
    ldsc = df.select(
        safe(snp_col).alias("SNP"),
        safe(a1_col).alias("A1"),
        safe(a2_col).alias("A2"),
        safe(frq_col).alias("FRQ"),
        safe(n_col_out).alias("N"),
        safe(beta_col).alias("BETA"),
        safe(se_col).alias("SE"),
        safe(p_col).alias("P"),
        safe(chr_col).alias("CHR"),
        safe(pos_col).alias("POS"),
        safe(info_col_out).alias("INFO"),
    )
    ldsc.write_csv(out, separator="\t")

def main():

    ap = argparse.ArgumentParser(description="Master GWAS QC -> LDSC-ready (polars)")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--pheno_id", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--sep", default="\t")
    ap.add_argument("--snp_col", required=True)
    ap.add_argument("--chr_col", required=True)
    ap.add_argument("--pos_col", required=True)
    ap.add_argument("--a1_col", required=True)
    ap.add_argument("--a2_col", required=True)
    ap.add_argument("--beta_col", required=True)
    ap.add_argument("--se_col", required=True)
    ap.add_argument("--p_col", required=True)
    ap.add_argument("--eaf_col", default=None)
    ap.add_argument("--freq_case_col", default=None)
    ap.add_argument("--freq_ctrl_col", default=None)
    ap.add_argument("--n_case_col", default=None)
    ap.add_argument("--n_ctrl_col", default=None)
    ap.add_argument("--n_col", default=None)
    ap.add_argument("--info_col", default=None)
    ap.add_argument("--require_info", action="store_true")
    ap.add_argument("--maf_min", type=float, default=0.01)
    ap.add_argument("--info_min", type=float, default=0.90)
    ap.add_argument("--exclude_mhc", action="store_true")
    ap.add_argument("--exclude_apoe", action="store_true")
    ap.add_argument("--apoe_chr", type=int, default=19)
    ap.add_argument("--apoe_start", type=int, default=44_000_000)
    ap.add_argument("--apoe_end", type=int, default=46_500_000)
    ap.add_argument("--drop_palindromes", action="store_true")
    ap.add_argument("--keep_snps_only", action="store_true")
    a = ap.parse_args()
    src = Path(a.inp)
    outdir = Path(a.outdir)
    out = outdir / a.pheno_id / "post-qc" / f"{a.pheno_id}_ldsc_ready.tsv"
    df = read_table(src, a.sep)
    count("loaded", df)
    df = as_upper(df, a.a1_col)
    df = as_upper(df, a.a2_col)
    df = to_int(df, a.chr_col)
    df = to_int(df, a.pos_col)
    df = to_float(df, a.beta_col)
    df = to_float(df, a.se_col)
    df = to_float(df, a.p_col)
    if a.n_col:
        df = to_float(df, a.n_col)
    if a.info_col:
        df = to_float(df, a.info_col)
    if a.eaf_col:
        df = to_float(df, a.eaf_col)
    if a.freq_case_col:
        df = to_float(df, a.freq_case_col)
    if a.freq_ctrl_col:
        df = to_float(df, a.freq_ctrl_col)
    if a.n_case_col:
        df = to_float(df, a.n_case_col)
    if a.n_ctrl_col:
        df = to_float(df, a.n_ctrl_col)
    count("after_types", df)
    if a.exclude_mhc:
        df = exclude_region(df, a.chr_col, a.pos_col, 6, 25_000_000, 34_000_000)
        count("after_exclude_MHC", df)
    if a.exclude_apoe:
        df = exclude_region(df, a.chr_col, a.pos_col, a.apoe_chr, a.apoe_start, a.apoe_end)
        count("after_exclude_APOE", df)
    if a.keep_snps_only:
        df = keep_snps_only(df, a.a1_col, a.a2_col)
        count("after_remove_indels", df)
    if a.drop_palindromes:
        df = drop_palindromes(df, a.a1_col, a.a2_col)
        count("after_remove_palindromes", df)
    frq_col = "__frq__"
    n_out_col = "__n__"
    info_out_col = "__info__"
    df = add_frq_maf(
        df,
        frq_col=frq_col,
        eaf_col=a.eaf_col,
        fcas_col=a.freq_case_col,
        fcon_col=a.freq_ctrl_col,
        ncas_col=a.n_case_col,
        ncon_col=a.n_ctrl_col,
    )
    df = add_n(df, n_out_col=n_out_col, n_col=a.n_col, ncas_col=a.n_case_col, ncon_col=a.n_ctrl_col)
    df = add_info(df, info_out_col=info_out_col, info_col=a.info_col, require_info=a.require_info)
    df = filter_maf_info(df, maf_min=a.maf_min, info_min=a.info_min, info_col=info_out_col, require_info=a.require_info)
    count("after_maf_info", df)
    required = [
        a.snp_col, a.a1_col, a.a2_col, frq_col,
        a.beta_col, a.se_col, a.p_col, a.chr_col, a.pos_col
    ]
    df = drop_missing_required(df, required)
    count("after_dropna_required", df)
    make_ldsc_and_write(
        df=df,
        out=out,
        snp_col=a.snp_col,
        a1_col=a.a1_col,
        a2_col=a.a2_col,
        frq_col=frq_col,
        n_col_out=n_out_col,
        beta_col=a.beta_col,
        se_col=a.se_col,
        p_col=a.p_col,
        chr_col=a.chr_col,
        pos_col=a.pos_col,
        info_col_out=info_out_col,
    )
    print("wrote:", out)

if __name__ == "__main__":
    main()