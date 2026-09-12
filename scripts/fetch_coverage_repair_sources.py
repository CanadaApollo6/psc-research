"""Fetch and pin supporting sources for the separate variant-coverage follow-up."""

import argparse
import json
from pathlib import Path

import fetch_comparison_sources as shared

ROOT=Path(__file__).resolve().parents[1]
MANIFEST=ROOT/'config/coverage-repair-sources.json'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--requests',type=Path)
    parser.add_argument('--offline',action='store_true')
    args=parser.parse_args()
    shared.MANIFEST=MANIFEST
    specs=json.loads(args.requests.read_text()) if args.requests else json.loads(MANIFEST.read_text())['sources']
    shared.acquire(specs,offline=args.offline)


if __name__=='__main__':main()
