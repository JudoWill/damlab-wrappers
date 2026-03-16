import glob
import os

import pysam
import yaml


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def count_fasta_records(path: str) -> int:
    count = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def count_fastq_records(path: str) -> int:
    with open(path) as f:
        lines = [l for l in f if l.strip()]
    assert len(lines) % 4 == 0, f"{path}: line count {len(lines)} not divisible by 4"
    return len(lines) // 4


def count_bam_records(path: str) -> int:
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        return sum(1 for _ in bam)


def load_yaml(path: str):
    with open(path) as f:
        content = f.read()
    assert content.startswith("# Cigarmath Split"), \
        f"{path}: YAML should start with '# Cigarmath Split'"
    return yaml.safe_load(content)


REQUIRED_YAML_KEYS = {
    "sample_name",
    "split_by",
    "output_format",
    "output_directory",
    "is_name_sorted",
    "top_n",
    "required_reference_overlap",
    "total_reads",
    "required_overlap_filtered",
    "categories",
    "unclassified",
}


# ===========================================================================
# Group 1 — test.sam-based: structure and sanity
# ===========================================================================

# --- flag / BAM output ------------------------------------------------------

def test_flag_bam_files_exist():
    assert os.path.isdir("test_output/flag_bam/"), "flag_bam output dir missing"
    bam_files = glob.glob("test_output/flag_bam/*.bam")
    assert len(bam_files) >= 1, "flag_bam output dir has no .bam files"
    assert os.path.exists("test_output/flag_bam_summary.yaml"), "flag_bam summary YAML missing"


def test_flag_bam_summary_structure():
    summary = load_yaml("test_output/flag_bam_summary.yaml")
    assert REQUIRED_YAML_KEYS == set(summary.keys()), \
        f"Unexpected YAML keys: {set(summary.keys()) ^ REQUIRED_YAML_KEYS}"
    assert summary["split_by"] == "flag"
    assert summary["output_format"] == "bam"
    assert summary["total_reads"] > 0
    assert summary["required_overlap_filtered"] == 0
    assert 0 <= summary["unclassified"] <= summary["total_reads"]


def test_flag_bam_record_totals():
    """Sum of BAM records across all output files equals total_reads in summary."""
    summary = load_yaml("test_output/flag_bam_summary.yaml")
    bam_files = glob.glob("test_output/flag_bam/*.bam")
    total_in_files = sum(count_bam_records(p) for p in bam_files)
    assert total_in_files == summary["total_reads"], (
        f"Total BAM records ({total_in_files}) != summary total_reads ({summary['total_reads']})"
    )


def test_flag_bam_category_counts_match_summary():
    """Per-category counts in summary match actual BAM file record counts."""
    summary = load_yaml("test_output/flag_bam_summary.yaml")
    for cat, expected_count in summary["categories"].items():
        path = f"test_output/flag_bam/{cat}.bam"
        assert os.path.exists(path), f"Expected output file {path} not found"
        actual = count_bam_records(path)
        assert actual == expected_count, \
            f"Category {cat!r}: BAM has {actual} records, summary says {expected_count}"


# --- query_size / FASTQ output ----------------------------------------------

def test_size_fastq_files_exist():
    assert os.path.isdir("test_output/size_fastq/"), "size_fastq output dir missing"
    fq_files = glob.glob("test_output/size_fastq/*.fq")
    assert len(fq_files) >= 1, "size_fastq output dir has no .fq files"
    assert os.path.exists("test_output/size_fastq_summary.yaml"), "size_fastq summary YAML missing"


def test_size_fastq_format():
    """Every .fq file has valid 4-line FASTQ records."""
    for path in glob.glob("test_output/size_fastq/*.fq"):
        with open(path) as f:
            lines = [l.rstrip("\n") for l in f if l.strip()]
        assert len(lines) % 4 == 0, f"{path}: line count {len(lines)} not divisible by 4"
        for i in range(0, len(lines), 4):
            assert lines[i].startswith("@"), f"{path} record {i//4}: header should start with @"
            assert lines[i + 2] == "+", f"{path} record {i//4}: third line should be +"
            assert len(lines[i + 1]) == len(lines[i + 3]), \
                f"{path} record {i//4}: seq and quality lengths differ"


def test_size_fastq_summary_crosscheck():
    """Sum of FASTQ record counts equals total_reads - required_overlap_filtered."""
    summary = load_yaml("test_output/size_fastq_summary.yaml")
    fq_files = glob.glob("test_output/size_fastq/*.fq")
    total_in_files = sum(count_fastq_records(p) for p in fq_files)
    expected = summary["total_reads"] - summary["required_overlap_filtered"]
    assert total_in_files == expected, (
        f"Total FASTQ records ({total_in_files}) != total_reads - filtered ({expected})"
    )


# --- required_reference_overlap + flag / FASTA output -----------------------

def test_overlap_fasta_files_exist():
    assert os.path.isdir("test_output/overlap_fasta/"), "overlap_fasta output dir missing"
    assert os.path.exists("test_output/overlap_fasta_summary.yaml"), \
        "overlap_fasta summary YAML missing"


def test_overlap_fasta_filtered_count():
    """Some reads are filtered; remaining reads are positive."""
    summary = load_yaml("test_output/overlap_fasta_summary.yaml")
    assert summary["required_overlap_filtered"] >= 0
    passing = summary["total_reads"] - summary["required_overlap_filtered"]
    assert passing > 0, "No reads passed the required_reference_overlap filter"


def test_overlap_fasta_format():
    """Every .fa file has valid FASTA records."""
    for path in glob.glob("test_output/overlap_fasta/*.fa"):
        with open(path) as f:
            lines = [l.rstrip("\n") for l in f if l.strip()]
        assert len(lines) >= 2, f"{path}: FASTA file appears empty"
        for i, line in enumerate(lines):
            if line.startswith(">"):
                assert i + 1 < len(lines), f"{path}: header at end of file with no sequence"
                assert not lines[i + 1].startswith(">"), \
                    f"{path}: consecutive headers with no sequence"


def test_overlap_fasta_summary_crosscheck():
    """FASTA record counts match summary categories + unclassified."""
    summary = load_yaml("test_output/overlap_fasta_summary.yaml")
    fa_files = glob.glob("test_output/overlap_fasta/*.fa")
    total_in_files = sum(count_fasta_records(p) for p in fa_files)
    expected = sum(summary["categories"].values()) + summary["unclassified"]
    assert total_in_files == expected, \
        f"Total FASTA records ({total_in_files}) != categories sum + unclassified ({expected})"


# --- top_n / FASTQ output ---------------------------------------------------

def test_top_n_category_count():
    """With top_n=1, at most 2 output files (1 category + unclassified)."""
    summary = load_yaml("test_output/top_n_fastq_summary.yaml")
    assert summary["top_n"] == 1
    fq_files = glob.glob("test_output/top_n_fastq/*.fq")
    assert len(fq_files) <= 2, \
        f"top_n=1 should produce ≤2 .fq files, got {len(fq_files)}: {fq_files}"


def test_top_n_category_names():
    """The single kept category must correspond to the most common flag group."""
    summary = load_yaml("test_output/top_n_fastq_summary.yaml")
    assert len(summary["categories"]) <= 1, \
        f"top_n=1 should have ≤1 named categories, got {summary['categories']}"


# ===========================================================================
# Group 2 — toy SAM: exact value assertions
# ===========================================================================

# --- toy CSV / FASTA output -------------------------------------------------

def test_toy_csv_output_files():
    """Exactly 3 FASTA files: alpha.fa, beta.fa, unclassified.fa."""
    fa_files = {os.path.basename(p) for p in glob.glob("test_output/toy_csv_fasta/*.fa")}
    assert fa_files == {"alpha.fa", "beta.fa", "unclassified.fa"}, \
        f"Unexpected FASTA files: {fa_files}"


def test_toy_csv_record_counts():
    """alpha: 1, beta: 1, unclassified: 1."""
    assert count_fasta_records("test_output/toy_csv_fasta/alpha.fa") == 1
    assert count_fasta_records("test_output/toy_csv_fasta/beta.fa") == 1
    assert count_fasta_records("test_output/toy_csv_fasta/unclassified.fa") == 1


def test_toy_csv_sequences():
    """alpha.fa contains all-A sequence; beta.fa all-T; unclassified.fa all-G."""
    for path, expected_base in [
        ("test_output/toy_csv_fasta/alpha.fa", "A"),
        ("test_output/toy_csv_fasta/beta.fa", "T"),
        ("test_output/toy_csv_fasta/unclassified.fa", "G"),
    ]:
        with open(path) as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith(">")]
        seq = "".join(lines)
        assert seq, f"{path}: no sequence found"
        assert all(c == expected_base for c in seq), \
            f"{path}: expected all {expected_base!r}, got {seq!r}"


def test_toy_csv_summary_exact():
    summary = load_yaml("test_output/toy_csv_fasta_summary.yaml")
    assert summary["split_by"] == "csv"
    assert summary["is_name_sorted"] is True, "toy SAM should be detected as name-sorted"
    assert summary["total_reads"] == 3
    assert summary["required_overlap_filtered"] == 0
    assert summary["categories"] == {"alpha": 1, "beta": 1}
    assert summary["unclassified"] == 1


# --- toy HP tag + overlap filter / BAM output --------------------------------

def test_toy_tag_summary_exact():
    """Overlap filter removes readC; HP_0 and HP_1 each get 1 read."""
    summary = load_yaml("test_output/toy_tag_bam_summary.yaml")
    assert summary["split_by"] == "tag"
    assert summary["is_name_sorted"] is True, "toy SAM should be detected as name-sorted"
    assert summary["total_reads"] == 3
    assert summary["required_overlap_filtered"] == 1, \
        "readC (200-220) should be filtered by HXB2F:3000-3020"
    assert summary["categories"] == {"0": 1, "1": 1}, \
        f"Expected categories {{'0': 1, '1': 1}}, got {summary['categories']}"
    assert summary["unclassified"] == 0


def test_toy_tag_bam_output_files():
    """Exactly 2 BAM files: 0.bam and 1.bam (raw HP tag values)."""
    bam_files = {os.path.basename(p) for p in glob.glob("test_output/toy_tag_bam/*.bam")}
    assert bam_files == {"0.bam", "1.bam"}, f"Unexpected BAM files: {bam_files}"


def test_toy_tag_bam_record_counts():
    """0.bam has 2 records (readA primary + supp); 1.bam has 2 records (readB)."""
    assert count_bam_records("test_output/toy_tag_bam/0.bam") == 2, \
        "HP_0 (readA) should contribute 2 BAM records (primary + supplementary)"
    assert count_bam_records("test_output/toy_tag_bam/1.bam") == 2, \
        "HP_1 (readB) should contribute 2 BAM records (primary + supplementary)"


def test_toy_tag_bam_read_names():
    """0.bam records are all readA; 1.bam records are all readB."""
    with pysam.AlignmentFile("test_output/toy_tag_bam/0.bam", "rb", check_sq=False) as bam:
        names = {seg.query_name for seg in bam}
    assert names == {"readA"}, f"0.bam should contain only readA, got {names}"

    with pysam.AlignmentFile("test_output/toy_tag_bam/1.bam", "rb", check_sq=False) as bam:
        names = {seg.query_name for seg in bam}
    assert names == {"readB"}, f"1.bam should contain only readB, got {names}"


def test_toy_tag_bam_pg_header():
    """Each output BAM has a PG entry with ID 'cigarmath-split'."""
    for path in glob.glob("test_output/toy_tag_bam/*.bam"):
        with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
            pg_ids = [pg["ID"] for pg in bam.header.to_dict().get("PG", [])]
        assert "cigarmath-split" in pg_ids, \
            f"{path}: PG entry 'cigarmath-split' not found. PG IDs: {pg_ids}"
