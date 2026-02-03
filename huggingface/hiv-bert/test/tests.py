import pytest
import os
import pandas as pd
import yaml

# Define build directory path
BUILD_DIR = os.path.join(os.path.dirname(__file__), 'test_output')

def test_output_exists():
    """Test that output files exist"""
    assert os.path.exists(os.path.join(BUILD_DIR, 'aa_v3_results.csv'))
    assert os.path.exists(os.path.join(BUILD_DIR, 'dna_v3_results.csv'))
    assert os.path.exists(os.path.join(BUILD_DIR, 'protein_results.csv'))
    assert os.path.exists(os.path.join(BUILD_DIR, 'edge_case_results.csv'))

def test_v3_sequence_detection():
    """Test that V3 sequences are processed correctly regardless of input type"""
    aa_df = pd.read_csv(os.path.join(BUILD_DIR, 'aa_v3_results.csv'))
    dna_df = pd.read_csv(os.path.join(BUILD_DIR, 'dna_v3_results.csv'))
    
    # Get predictions for equivalent sequences
    aa_pred = aa_df[aa_df['id'] == 'V3_example'].iloc[0]
    dna_pred = dna_df[dna_df['id'] == 'V3_example_DNA'].iloc[0]
    
    # Predictions should be similar for equivalent sequences
    prob_cols = [col for col in aa_df.columns if col != 'id']
    for col in prob_cols:
        assert abs(aa_pred[col] - dna_pred[col]) < 0.01

def test_protein_processing():
    """Test processing of protein sequences"""
    df = pd.read_csv(os.path.join(BUILD_DIR, 'protein_results.csv'))
    
    # Check that protein sequences are present
    assert 'Tat_example' in df['id'].values
    assert 'PR_example' in df['id'].values
    
    # Check embedding dimensions for protein sequences
    if 'dim_0' in df.columns:  # If using embedding model
        embedding_cols = [col for col in df.columns if col.startswith('dim_')]
        assert len(embedding_cols) == 1024
        assert df[embedding_cols].notnull().all().all()

def test_edge_cases():
    """Test handling of edge cases"""
    df = pd.read_csv(os.path.join(BUILD_DIR, 'edge_case_results.csv'))
    
    # Short sequence should be filtered out
    assert 'AA_short' not in df['id'].values
    
    # Invalid DNA should still be processed
    assert 'DNA_invalid_chars' in df['id'].values
    
    # Mixed sequence should be processed
    assert 'Mixed_sequence' in df['id'].values

def test_case_sensitivity():
    """Test that case doesn't affect predictions"""
    df = pd.read_csv(os.path.join(BUILD_DIR, 'dna_v3_results.csv'))
    
    # Get predictions for normal and mixed case sequences
    normal_pred = df[df['id'] == 'V3_Tcell_DNA_frame0'].iloc[0]
    mixed_pred = df[df['id'] == 'V3_Mixed_Case_DNA'].iloc[0]
    
    # Predictions should be identical regardless of case
    prob_cols = [col for col in df.columns if col != 'id']
    for col in prob_cols:
        assert abs(normal_pred[col] - mixed_pred[col]) < 0.01

def test_metrics_output():
    """Test metrics file format and content"""
    # Test bodysite metrics
    with open(os.path.join(BUILD_DIR, 'aa_v3_metrics.yaml'), 'r') as f:
        metrics = yaml.safe_load(f)
        
        # Check structure
        assert 'total_sequences' in metrics
        assert 'processed_sequences' in metrics
        assert 'mean_values' in metrics
        assert 'high_confidence_counts' in metrics
        assert 'model_name' in metrics
        assert 'model_type' in metrics
        assert 'parameters' in metrics
        
        # Check content
        assert metrics['total_sequences'] > 0
        assert metrics['processed_sequences'] > 0
        assert metrics['model_type'] == 'classification'
        
        # Check mean values
        for label in ['periphery-tcell', 'CNS', 'lung']:  # Example labels
            assert label in metrics['mean_values']
            assert 0 <= metrics['mean_values'][label] <= 1
            
        # Check high confidence counts
        for label in ['periphery-tcell', 'CNS', 'lung']:
            assert label in metrics['high_confidence_counts']
            assert isinstance(metrics['high_confidence_counts'][label], int)
            assert metrics['high_confidence_counts'][label] >= 0

