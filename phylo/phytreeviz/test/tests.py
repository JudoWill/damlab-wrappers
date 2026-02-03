import pytest # type: ignore
import os
from PIL import Image # type: ignore

def test_plot_output_exists():
    """Test plot exists"""
    assert os.path.exists('test_output/test.png')

def test_plot_is_valid():
    """Test plot is valid image"""
    try:
        img = Image.open('test_output/test.png')
        img.verify()
        assert True
    except:
        assert False, "Invalid image file" 