import os
import pandas as pd
import yaml
import pytest

@pytest.fixture
def results_dir():
    return "results"

@pytest.fixture
def csv_file(results_dir):
    return os.path.join(results_dir, "matches.csv")

@pytest.fixture
def metrics_file(results_dir):
    return os.path.join(results_dir, "metrics.yaml")

def test_csv_output_exists(csv_file):
    """Test that CSV output file exists"""
    assert os.path.exists(csv_file), "CSV output file not found"

def test_metrics_output_exists(metrics_file):
    """Test that metrics output file exists"""
    assert os.path.exists(metrics_file), "Metrics output file not found"

def test_csv_content(csv_file):
    """Test CSV output content"""
    df = pd.read_csv(csv_file)
    
    # Check column names
    expected_columns = ['read_name', 'start_pattern', 'end_pattern']
    assert all(col in df.columns for col in expected_columns), "Missing expected columns"
    
    # Check number of rows
    assert len(df) == 4, "Expected 4 sequences in output"
    
    # Check specific matches
    seq1_row = df[df['read_name'] == 'seq1']
    assert not seq1_row['start_pattern'].isna().iloc[0], "Pattern 1 should match seq1"
    assert not seq1_row['end_pattern'].isna().iloc[0], "Pattern 2 should match seq1"
    
    seq2_row = df[df['read_name'] == 'seq2']
    assert not seq2_row['start_pattern'].isna().iloc[0], "Pattern 1 should match seq2"
    assert seq2_row['end_pattern'].isna().iloc[0], "Pattern 2 should not match seq2"
    
    seq3_row = df[df['read_name'] == 'seq3']
    assert seq3_row['start_pattern'].isna().iloc[0], "Pattern 1 should not match seq3"
    assert not seq3_row['end_pattern'].isna().iloc[0], "Pattern 2 should match seq3"
    
    seq4_row = df[df['read_name'] == 'seq4']
    assert seq4_row['start_pattern'].isna().iloc[0], "Pattern 1 should not match seq4"
    assert seq4_row['end_pattern'].isna().iloc[0], "Pattern 2 should not match seq4"

def test_metrics_content(metrics_file):
    """Test metrics output content"""
    with open(metrics_file) as f:
        metrics = yaml.safe_load(f)
    
    # Check metrics structure
    assert 'total_reads' in metrics, "Missing total_reads in metrics"
    assert 'pattern_matches' in metrics, "Missing pattern_matches in metrics"
    assert 'start_pattern' in metrics['pattern_matches'], "Missing pattern_1 count in metrics"
    assert 'end_pattern' in metrics['pattern_matches'], "Missing pattern_2 count in metrics"
    
    # Check metrics values
    assert metrics['total_reads'] == 4, "Expected 4 total reads"
    assert metrics['pattern_matches']['start_pattern'] == 2, "Expected 2 matches for start_pattern"
    assert metrics['pattern_matches']['end_pattern'] == 2, "Expected 2 matches for end_pattern" 