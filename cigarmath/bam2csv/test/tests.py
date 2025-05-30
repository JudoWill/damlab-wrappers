"""Unit tests for bam2csv functions"""

import csv

def test_data_extracted():
    """Test that the data is extracted correctly"""

    field_order = ["query_name", "reference_start", "CR", "OX"]
    correct_data = [
        ("568f3082-407b-42d7-a287-94565d242d4f", "1134", "ACGCGCATCGAGCCCATGTGCACAGCTCGCGTAT", "TCTTGCATCATAGCGTTC"),
        ("0c2a9d1a-1943-47bf-b43c-3a440d31a61d", "1132", "ACGCGGAAGCCATGACTATAGTATCGATGCGCGC", "ACGTGAAACAAATTGACA"),
        ("0c2a9d1a-1943-47bf-b43c-3a440d31a61d", "6154", "ACGCGGAAGCCATGACTATAGTATCGATGCGCGC", ""),
        ("f64d3b26-b81b-4990-a1cf-5e57f1bd0fac", "1100", "ACGCGCGCATCGGGTATATGCATGGCTCGCGTAT", ""),
        ("a948acbf-a190-4ad3-bd2a-19a731c13c72", "1131", "ACGCGAGCCATGACTATAGTATCGATGCGCGCGT", ""),
        ("a948acbf-a190-4ad3-bd2a-19a731c13c72", "6154", "ACGCGAGCCATGACTATAGTATCGATGCGCGCGT", ""),
        ("638f8412-9031-4cf5-a074-db6e5a80c146", "6093", "ACGCGCGCATCGAGCCCATGTGCATGGCTCGCGT", ""),
        ("4c6dae94-a7a5-47ad-a5c5-948847d6ecff", "1132", "ACGCGAGCCATGACTATAGTATCGATGCGCGTAT", "CACTACTACGCTACAAGC"),
        ("4c6dae94-a7a5-47ad-a5c5-948847d6ecff", "6154", "ACGCGAGCCATGACTATAGTATCGATGCGCGTAT", ""),
        ("9d3abf22-fa9b-4e6b-8447-5d5deea42a52", "1132", "ACGCGAGCCATGACTATAGTATCGATGCGCGCGT", "")
    ]

    with open('test_output/test_output.csv', newline='') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            print(row)
            print(correct_data[i])
            print()

            for j, field in enumerate(field_order):
                assert row[field] == correct_data[i][j]

