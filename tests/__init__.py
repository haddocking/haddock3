"""Test module."""
from pathlib import Path


data_path = Path(__file__).resolve().parents[0]
golden_data = Path(data_path, 'golden_data')
