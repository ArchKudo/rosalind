import polars as pl
import main


def test_is_snv_true():
    class DummyRec:
        ref = "A"
        alts = ["T"]

    assert main.is_snv(DummyRec())


def test_is_snv_false():
    class DummyRec:
        ref = "AT"
        alts = ["T"]

    assert not main.is_snv(DummyRec())


def test_build_dataframe():
    records = [
        ("chr1", 100, "A", "T", 1, "A"),
        ("chr1", 100, "A", "T", 2, "T"),
    ]
    df = main.build_dataframe(records)
    assert isinstance(df, pl.DataFrame)
    assert df.shape == (2, 6)
    assert set(df.columns) == {"chrom", "pos", "ref", "alt", "hap", "base"}


def test_aggregate_counts():
    records = [
        ("chr1", 100, "A", "T", 1, "A"),
        ("chr1", 100, "A", "T", 1, "T"),
        ("chr1", 100, "A", "T", 2, "A"),
        ("chr1", 100, "A", "T", 2, "T"),
        ("chr1", 100, "A", "T", 2, "T"),
    ]
    df = main.build_dataframe(records)
    agg = main.aggregate_counts(df)
    row = agg.row(0)
    # Should have 1 h1_ref, 1 h1_alt, 1 h2_ref, 2 h2_alt
    assert row[agg.columns.index("h1_ref")] == 1
    assert row[agg.columns.index("h1_alt")] == 1
    assert row[agg.columns.index("h2_ref")] == 1
    assert row[agg.columns.index("h2_alt")] == 2
