import pytest
import os
import csv
import yaml


def test_multiqc_report_exists():
    assert os.path.exists("build/multiqc/multiqc_report.html")


def test_multiqc_report_general_stats_exists():
    assert os.path.exists("build/multiqc/multiqc_report_data/multiqc_general_stats.txt")

def get_general_stats_data():
    with open("build/multiqc/multiqc_report_data/multiqc_general_stats.txt", "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def test_strainline():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that the wanted columns are present
    wanted_cols = ["strainline-haplotype_count", "strainline-haplotype_freqs_max"]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check that the data is correct
    assert [row['strainline-haplotype_count'] for row in data] == ['3', '1', '2']
    assert [row['strainline-haplotype_freqs_max'] for row in data] == ['0.75', '1.0', '0.5']


def test_dorado():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that all expected columns are present
    wanted_cols = [
        "dorado-total_reads",
        "dorado-duplex_reads",
        "dorado-duplex_rate",
        "dorado-simplex_qscore",
        "dorado-duplex_qscore"
    ]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check sample1 data (5 reads: 2 duplex, 1 simplex with duplex, 1 simplex no duplex, 1 no tag)
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    assert float(sample1['dorado-total_reads']) == 5
    assert float(sample1['dorado-duplex_reads']) == 2
    assert float(sample1['dorado-duplex_rate']) == pytest.approx(0.4)  # 2/5
    assert float(sample1['dorado-simplex_qscore']) == pytest.approx(13.27, abs=0.1)  # Mean of 12.5, 14.2, 13.1
    assert float(sample1['dorado-duplex_qscore']) == pytest.approx(18.95, abs=0.1)  # Mean of 18.7, 19.2
    
    # Check sample2 data (6 reads: 2 duplex, 2 simplex with duplex, 1 simplex no duplex, 1 no tag)
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    assert float(sample2['dorado-total_reads']) == 6
    assert float(sample2['dorado-duplex_reads']) == 2
    assert float(sample2['dorado-duplex_rate']) == pytest.approx(0.333, abs=0.01)  # 2/6
    
    # Check sample3 data (7 reads: 3 duplex, 2 simplex with duplex, 1 simplex no duplex, 1 no tag)
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    assert float(sample3['dorado-total_reads']) == 7
    assert float(sample3['dorado-duplex_reads']) == 3
    assert float(sample3['dorado-duplex_rate']) == pytest.approx(0.429, abs=0.01)  # 3/7


def test_generic():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that the wanted column is present
    assert 'generic-gel_band' in data[0]
    
    # Check that the data is correct for each sample
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    
    assert sample1['generic-gel_band'] == 'yes'
    assert sample2['generic-gel_band'] == 'yes'
    assert sample3['generic-gel_band'] == 'no'


def test_primercheck():
    data = get_general_stats_data()

    assert len(data) == 3

    # Check that all expected columns are present
    wanted_cols = [
        "primercheck-total_sequences",
        "primercheck-primer_primer1",
        "primercheck-primer_primer2",
    ]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check sample1 data (4/5 reads with primer1, 1/5 with primer2/3)
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    assert float(sample1['primercheck-total_sequences']) == 5
    assert float(sample1['primercheck-primer_primer1']) == pytest.approx(0.8)  # 4/5 reads
    assert float(sample1['primercheck-primer_primer2']) == pytest.approx(0.2)  # 1/5 reads
    
    # Check sample2 data (2/5 with primer1, 3/5 with primer2/3)
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    assert float(sample2['primercheck-total_sequences']) == 5
    assert float(sample2['primercheck-primer_primer1']) == pytest.approx(0.4)  # 2/5 reads
    assert float(sample2['primercheck-primer_primer2']) == pytest.approx(0.6)  # 3/5 reads
    
    # Check sample3 data (all reads with primer2/3)
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    assert float(sample3['primercheck-total_sequences']) == 4
    assert float(sample3['primercheck-primer_primer1']) == pytest.approx(0.0)  # No reads
    assert float(sample3['primercheck-primer_primer2']) == pytest.approx(1.0)  # All reads


def test_intactness():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that all expected columns are present
    wanted_cols = [
        "intactness-total_reads",
        "intactness-countable_reads",
        "intactness-intact_reads",
        "intactness-percent_countable",
        "intactness-percent_intact"
    ]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check sample1 data (5 reads: 2 countable, 1 intact)
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    assert float(sample1['intactness-total_reads']) == 5
    assert float(sample1['intactness-countable_reads']) == 3
    assert float(sample1['intactness-intact_reads']) == 1
    assert float(sample1['intactness-percent_countable']) == pytest.approx(60.0)  # 3/5 * 100
    assert float(sample1['intactness-percent_intact']) == pytest.approx(33.33, abs=0.01)  # 1/3 * 100
    
    # Check sample2 data (5 reads: 3 countable, 2 intact)
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    assert float(sample2['intactness-total_reads']) == 5
    assert float(sample2['intactness-countable_reads']) == 3
    assert float(sample2['intactness-intact_reads']) == 2
    assert float(sample2['intactness-percent_countable']) == pytest.approx(60.0)  # 3/5 * 100
    assert float(sample2['intactness-percent_intact']) == pytest.approx(66.67, abs=0.01)  # 2/3 * 100
    
    # Check sample3 data (4 reads: 4 countable, 1 intact)
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    assert float(sample3['intactness-total_reads']) == 4
    assert float(sample3['intactness-countable_reads']) == 4
    assert float(sample3['intactness-intact_reads']) == 1
    assert float(sample3['intactness-percent_countable']) == pytest.approx(100.0)  # 4/4 * 100
    assert float(sample3['intactness-percent_intact']) == pytest.approx(25.0)  # 1/4 * 100


def test_hivbert():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that all expected columns are present
    wanted_cols = [
        "hivbert-total_sequences",
        "hivbert-processed_sequences",
        "hivbert-processed_rate"
    ]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check sample1 data (3 sequences: 2 CCR5, 1 dual)
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    assert float(sample1['hivbert-total_sequences']) == 3
    assert float(sample1['hivbert-processed_sequences']) == 3
    assert float(sample1['hivbert-processed_rate']) == 1.0
    
    # Check sample2 data (3 sequences: 2 CXCR4, 1 CCR5)
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    assert float(sample2['hivbert-total_sequences']) == 3
    assert float(sample2['hivbert-processed_sequences']) == 3
    assert float(sample2['hivbert-processed_rate']) == 1.0
    
    # Check sample3 data (3 sequences: 2 dual, 1 CXCR4)
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    assert float(sample3['hivbert-total_sequences']) == 3
    assert float(sample3['hivbert-processed_sequences']) == 3
    assert float(sample3['hivbert-processed_rate']) == 1.0


def test_hivbert_metrics_file():
    """Test that metrics files contain expected values"""
    # Check sample1 metrics
    with open('build/hivbert/sample1.hivbert.yaml', 'r') as f:
        metrics = yaml.safe_load(f)
        assert metrics['model_name'] == 'damlab/HIV_V3_coreceptor'
        assert metrics['model_type'] == 'classification'
        assert metrics['total_sequences'] == 3
        assert metrics['processed_sequences'] == 3
        # Check high confidence counts
        assert metrics['high_confidence_counts']['CCR5'] == 2  # 2 CCR5 sequences
        assert metrics['high_confidence_counts']['CCR5+CXCR4'] == 1  # 1 dual tropic
    
    # Check sample2 metrics
    with open('build/hivbert/sample2.hivbert.yaml', 'r') as f:
        metrics = yaml.safe_load(f)
        assert metrics['high_confidence_counts']['CXCR4'] == 2  # 2 CXCR4 sequences
        assert metrics['high_confidence_counts']['CCR5'] == 1  # 1 CCR5 sequence
    
    # Check sample3 metrics
    with open('build/hivbert/sample3.hivbert.yaml', 'r') as f:
        metrics = yaml.safe_load(f)
        assert metrics['high_confidence_counts']['CCR5+CXCR4'] == 2  # 2 dual tropic
        assert metrics['high_confidence_counts']['CXCR4'] == 1  # 1 CXCR4 sequence


def test_slice():
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that the wanted columns are present
    wanted_cols = [
        "slice-pct_overlapping",
        "slice-overlapping",
        "slice-total_segments"
    ]
    assert all(col in data[0] for col in wanted_cols)
    
    # Check sample1 data for gag region (100-200)
    # In sample1, 5 out of 7 reads overlap this region
    sample1 = next(row for row in data if row['Sample'] == 'sample1')
    gag_pct = float(sample1['slice-pct_overlapping'])
    assert pytest.approx(gag_pct, abs=1.0) == 71.4  # ~5/7 * 100 = 71.4%
    
    # Check sample2 data
    # In sample2, 8 out of 9 reads overlap either gag or pol regions
    sample2 = next(row for row in data if row['Sample'] == 'sample2')
    assert float(sample2['slice-overlapping']) > 0
    
    # Check sample3 data
    # In sample3, 3 reads overlap gag and 1 read overlaps pol
    sample3 = next(row for row in data if row['Sample'] == 'sample3')
    assert float(sample3['slice-overlapping']) > 0
    assert float(sample3['slice-total_segments']) == 5


def test_deletion_block_detection():
    """Test that deletion block detection columns appear in general stats."""
    data = get_general_stats_data()
    assert len(data) == 3

    # Check that all expected columns are present
    wanted_cols = [
        "deletion_block_detection-deletion_frequency",
        "deletion_block_detection-reads_with_deletions",
        "deletion_block_detection-unique_deletion_count",
    ]
    assert all(col in data[0] for col in wanted_cols), \
        f"Missing columns. Available: {list(data[0].keys())}"
    
    # Check that values are reasonable for each sample
    for row in data:
        sample = row['Sample']
        
        # Deletion frequency should be between 0 and 1
        del_freq = float(row['deletion_block_detection-deletion_frequency'])
        assert 0.0 <= del_freq <= 1.0, \
            f"Sample {sample}: deletion_frequency {del_freq} not in [0, 1]"
        
        # Reads with deletions should be non-negative
        reads_with_dels = float(row['deletion_block_detection-reads_with_deletions'])
        assert reads_with_dels >= 0, \
            f"Sample {sample}: reads_with_deletions should be non-negative"
        
        # Unique deletion count should be non-negative
        unique_dels = float(row['deletion_block_detection-unique_deletion_count'])
        assert unique_dels >= 0, \
            f"Sample {sample}: unique_deletion_count should be non-negative"


def test_deletion_block_detection_output_files():
    """Test that deletion block detection output files are created correctly."""
    for sample in ["sample1", "sample2", "sample3"]:
        summary_path = f'build/deletion_block_detection/{sample}.summary.yaml'
        reads_path = f'build/deletion_block_detection/{sample}.reads.csv'
        deletions_path = f'build/deletion_block_detection/{sample}.deletions.csv'
        
        # Check summary YAML
        assert os.path.exists(summary_path), f"Missing {summary_path}"
        with open(summary_path, 'r') as f:
            content = f.read()
            assert content.startswith('# Cigarmath Deletion Block Detection'), \
                f"YAML should start with MultiQC marker comment"
            metrics = yaml.safe_load(content)
        
        assert metrics['sample_name'] == sample
        assert metrics['total_reads'] > 0
        assert metrics['reads_with_deletions'] >= 0
        assert metrics['unique_deletion_count'] >= 0
        assert metrics['total_deletion_count'] >= 0
        assert 0.0 <= metrics['deletion_frequency'] <= 1.0
        
        # Check reads CSV exists and has correct structure
        assert os.path.exists(reads_path), f"Missing {reads_path}"
        with open(reads_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) > 0, f"reads.csv for {sample} should have data"
            expected_cols = {'read_name', 'reference_start', 'reference_end', 'deletions'}
            assert set(rows[0].keys()) == expected_cols, \
                f"reads.csv columns mismatch for {sample}"
        
        # Check deletions CSV exists and has correct structure
        assert os.path.exists(deletions_path), f"Missing {deletions_path}"
        with open(deletions_path, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            if len(rows) > 0:
                expected_cols = {'deletion_start', 'deletion_end', 'deletion_size', 
                               'read_count', 'coverage_count'}
                assert set(rows[0].keys()) == expected_cols, \
                    f"deletions.csv columns mismatch for {sample}"


def test_deletion_block_detection_plots():
    """Test that deletion block detection plot data is generated in MultiQC output."""
    # Check that the MultiQC report data directory contains plot data
    plot_data_dir = "build/multiqc/multiqc_report_data"
    assert os.path.exists(plot_data_dir), "MultiQC report data directory should exist"
    
    # Check for deletion_block_detection data file
    expected_data_file = os.path.join(plot_data_dir, "multiqc_deletion_block_detection.txt")
    assert os.path.exists(expected_data_file), \
        f"Expected deletion_block_detection data file at {expected_data_file}"
    
    # Verify the data file has content for all samples
    with open(expected_data_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    
    sample_names = {row['Sample'] for row in rows}
    assert 'sample1' in sample_names, "sample1 should be in plot data"
    assert 'sample2' in sample_names, "sample2 should be in plot data"
    assert 'sample3' in sample_names, "sample3 should be in plot data"

