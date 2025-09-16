import pandas as pd
from Bio import SeqIO
from pathlib import Path
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chimeras", type=Path, required=True
    )
    parser.add_argument(
        "--db", type=Path, required=True
    )
    parser.add_argument(
        "--out", type=Path, required=True
    )
    return parser.parse_args()

args = parse_args()

df = pd.read_csv(args.chimeras)

db_targets = df['subject'].tolist()

chimeras_records = []
with open(args.db, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        if record.id in db_targets:
            chimeras_records.append(record)

with open(args.out, "w") as f:
    SeqIO.write(chimeras_records, f, "fasta")
