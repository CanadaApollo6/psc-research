"""Load a local research API key without printing or putting it in Git."""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
KEY_FILE = ROOT / 'data/private/alphagenome_api_key'


def load_key():
    key = os.environ.get('ALPHAGENOME_API_KEY', '').strip()
    if not key and KEY_FILE.is_file():
        key = KEY_FILE.read_text().strip()
    if not key:
        raise RuntimeError('No AlphaGenome API key configured. Set ALPHAGENOME_API_KEY or use data/private/alphagenome_api_key.')
    if any(char.isspace() for char in key):
        raise ValueError('The key must be a single value, with no prefix or internal whitespace.')
    return key


def key_is_configured():
    try:
        load_key()
        return True
    except (RuntimeError, ValueError):
        return False
