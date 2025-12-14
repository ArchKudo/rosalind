import pysam
import polars as pl
import pyximport
from typing import List, Tuple, Any

pyximport.install(language_level=3)

from extract_reads import extract_read_info # type: ignore


def is_snv(record: Any) -> bool:
    """
    Ignore anything but lines containing single character REF,ALT
    """
    return len(record.ref) == len(record.alts) == len(record.alts[0]) == 1


def load_variants(vcf_path: str) -> List[Tuple[str, int, str, str]]:
    """
    Get the first two columns of output
    And REF/ALT for searching in bam file
    """
    vcf: pysam.VariantFile = pysam.VariantFile(vcf_path)
    return [
        (record.chrom, record.pos, record.ref, record.alts[0]) for record in vcf.fetch() if is_snv(record)
    ]


def build_dataframe(
    read_records: List[Tuple[str, int, str, str, int, str]],
) -> pl.DataFrame:
    """
    Build dataframe from the extracted reads
    """
    return pl.DataFrame(
        read_records, schema=["chrom", "pos", "ref", "alt", "hap", "base"], orient="row"
    )


def aggregate_counts(df: pl.DataFrame) -> pl.DataFrame:
    """
    Count the remaining 4 columns of output
    """
    return (
        df.with_columns(
            [
                (pl.col("base") == pl.col("ref")).alias("is_ref"),
                (pl.col("base") == pl.col("alt")).alias("is_alt"),
            ]
        )
        .group_by(["chrom", "pos"])
        .agg(
            [
                pl.col("hap")
                .filter((pl.col("hap") == 1) & pl.col("is_ref"))
                .count()
                .alias("h1_ref"),
                pl.col("hap")
                .filter((pl.col("hap") == 1) & pl.col("is_alt"))
                .count()
                .alias("h1_alt"),
                pl.col("hap")
                .filter((pl.col("hap") == 2) & pl.col("is_ref"))
                .count()
                .alias("h2_ref"),
                pl.col("hap")
                .filter((pl.col("hap") == 2) & pl.col("is_alt"))
                .count()
                .alias("h2_alt"),
            ]
        )
        .sort("pos")
    )


def main(bam_path: str, vcf_path: str, output_path: str) -> pl.DataFrame:
    variant_records: List[Tuple[str, int, str, str]] = load_variants(vcf_path)
    df: pl.DataFrame = extract_read_info(bam_path, variant_records)
    agg_df: pl.DataFrame = aggregate_counts(df)
    agg_df.write_csv(output_path, separator="\t")
    return agg_df


if __name__ == "__main__":
    main(
        "test_data/giab_2023.05.hg002.haplotagged.chr16_28000000_29000000.processed.30x.bam",
        "test_data/giab_2023.05.hg002.wf_snp.chr16_28000000_29000000.vcf.gz",
        "output_main.tsv",
    )
