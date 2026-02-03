"""Unit tests for pandas merge wrapper"""

import os
import pytest
import pandas as pd
from pathlib import Path

def test_output_exists():
    """Test that output files were created"""
    assert os.path.exists('test_output/inner_merge.csv'), "Inner merge output not found"
    assert os.path.exists('test_output/outer_merge.csv'), "Outer merge output not found"
    assert os.path.exists('test_output/no_on_merge.csv'), "No on merge output not found"
    assert os.path.exists('test_output/three_file_merge.csv'), "Three file merge output not found"
    assert os.path.exists('test_output/different_columns_merge.csv'), "Different columns merge output not found"
    assert os.path.exists('test_output/three_file_different_columns_merge.csv'), "Three file different columns merge output not found"

def test_inner_merge():
    """Test inner merge with common column"""
    df = pd.read_csv('test_output/inner_merge.csv')
    assert len(df) > 0, "Inner merge produced empty dataframe"
    assert 'id' in df.columns, "Common column 'id' not found in output"
    
    # Check that only matching rows are included
    df1 = pd.read_csv('test_data1.csv')
    df2 = pd.read_csv('test_data2.csv')
    expected_rows = len(pd.merge(df1, df2, on='id', how='inner'))
    assert len(df) == expected_rows, "Inner merge row count mismatch"

def test_outer_merge():
    """Test outer merge with common column"""
    df = pd.read_csv('test_output/outer_merge.csv')
    assert len(df) > 0, "Outer merge produced empty dataframe"
    assert 'id' in df.columns, "Common column 'id' not found in output"
    
    # Check that all rows are included
    df1 = pd.read_csv('test_data1.csv')
    df2 = pd.read_csv('test_data2.csv')
    expected_rows = len(pd.merge(df1, df2, on='id', how='outer'))
    assert len(df) == expected_rows, "Outer merge row count mismatch"

def test_no_on_merge():
    """Test merge without specifying 'on' parameter"""
    df = pd.read_csv('test_output/no_on_merge.csv')
    assert len(df) > 0, "No on merge produced empty dataframe"
    
    # Check that all columns from both dataframes are present
    df1 = pd.read_csv('test_data1.csv')
    df2 = pd.read_csv('test_data2.csv')
    expected_columns = set(df1.columns) | set(df2.columns)
    assert set(df.columns) == expected_columns, "Column mismatch in no on merge"

def test_three_file_merge():
    """Test merging three files with suffixes"""
    df = pd.read_csv('test_output/three_file_merge.csv')
    assert len(df) > 0, "Three file merge produced empty dataframe"
    
    # Check that all rows are included
    df1 = pd.read_csv('test_data1.csv')
    df2 = pd.read_csv('test_data2.csv')
    df3 = pd.read_csv('test_data3.csv')
    
    # Merge manually to get expected result
    temp = pd.merge(df1, df2, on='id', how='outer', suffixes=('_first', '_second'))
    expected = pd.merge(temp, df3, on='id', how='outer', suffixes=('', '_third'))
    expected_rows = len(expected)
    assert len(df) == expected_rows, "Three file merge row count mismatch"
    
    # Check that overlapping columns have correct suffixes
    overlapping_cols = set(df1.columns) & set(df2.columns) & set(df3.columns)
    for col in overlapping_cols:
        if col != 'id':  # Skip the merge column
            assert f"{col}_first" in df.columns, f"Expected column {col}_first not found"
            assert f"{col}_second" in df.columns, f"Expected column {col}_second not found"
            assert f"{col}_third" in df.columns, f"Expected column {col}_third not found"

def test_different_columns_merge():
    """Test merging files with different column names using list of lists"""
    df = pd.read_csv('test_output/different_columns_merge.csv')
    assert len(df) > 0, "Different columns merge produced empty dataframe"
    
    # Check that all rows are included
    df1 = pd.read_csv('test_data4.csv')
    df2 = pd.read_csv('test_data5.csv')
    
    # Merge manually to get expected result
    expected = pd.merge(df1, df2, left_on='employee_id', right_on='staff_id', how='outer')
    expected_rows = len(expected)
    assert len(df) == expected_rows, "Different columns merge row count mismatch"
    
    # Check that both ID columns are present
    assert 'employee_id' in df.columns, "employee_id column not found"
    assert 'staff_id' in df.columns, "staff_id column not found"
    
    # Check that the merge worked correctly by comparing a few rows
    sample_row = df[df['employee_id'] == 1].iloc[0]
    assert sample_row['name'] == 'John', "Name mismatch in merged data"
    assert sample_row['salary'] == 75000, "Salary mismatch in merged data"

def test_three_file_different_columns_merge():
    """Test merging three files with different column names using list of lists"""
    df = pd.read_csv('test_output/three_file_different_columns_merge.csv')
    assert len(df) > 0, "Three file different columns merge produced empty dataframe"
    
    # Check that all rows are included
    df1 = pd.read_csv('test_data4.csv')
    df2 = pd.read_csv('test_data5.csv')
    df3 = pd.read_csv('test_data6.csv')
    
    # Merge manually to get expected result
    temp = pd.merge(df1, df2, left_on='employee_id', right_on='staff_id', how='outer')
    expected = pd.merge(temp, df3, left_on='employee_id', right_on='person_id', how='outer')
    expected_rows = len(expected)
    assert len(df) == expected_rows, "Three file different columns merge row count mismatch"
    
    # Check that all ID columns are present
    assert 'employee_id' in df.columns, "employee_id column not found"
    assert 'staff_id' in df.columns, "staff_id column not found"
    assert 'person_id' in df.columns, "person_id column not found"
    
    # Check that the merge worked correctly by comparing a few rows
    sample_row = df[df['employee_id'] == 1].iloc[0]
    assert sample_row['name'] == 'John', "Name mismatch in merged data"
    assert sample_row['salary'] == 75000, "Salary mismatch in merged data"
    assert sample_row['position'] == 'Senior Engineer', "Position mismatch in merged data"

def test_three_file_same_key_merge():
    """Test merging three files where the middle table has a different key name"""
    df = pd.read_csv('test_output/three_file_same_key_merge.csv')
    print(df)
    assert len(df) > 0, "Three file same key merge produced empty dataframe"
    
    # Check that the key columns are preserved
    assert 'key' in df.columns, "Key column not found in output"
    assert 'key2' in df.columns, "Key2 column not found in output"
    
    # Check that all value columns are present
    assert 'value1' in df.columns, "value1 column not found"
    assert 'value2' in df.columns, "value2 column not found"
    assert 'value3' in df.columns, "value3 column not found"
    
    # Check that we have the expected number of rows
    assert len(df) == 3, "Expected 3 rows in output"
    
    # Check the values directly
    expected_data = [
        {'key': 1, 'key2': 1, 'value1': 'A', 'value2': 'X', 'value3': 'M'},
        {'key': 2, 'key2': 2, 'value1': 'B', 'value2': 'Y', 'value3': 'N'},
        {'key': 3, 'key2': 3, 'value1': 'C', 'value2': 'Z', 'value3': 'O'}
    ]
    
    for expected_row in expected_data:
        actual_row = df[df['key'] == expected_row['key']].iloc[0]
        for col, expected_value in expected_row.items():
            actual_value = actual_row[col]
            if isinstance(expected_value, int):
                # For integer keys, compare as integers
                assert int(actual_value) == expected_value, f"Value mismatch for {col}: expected {expected_value}, got {actual_value}"
            else:
                # For string values, compare as strings
                assert str(actual_value).strip() == expected_value, f"Value mismatch for {col}: expected '{expected_value}', got '{actual_value}'" 