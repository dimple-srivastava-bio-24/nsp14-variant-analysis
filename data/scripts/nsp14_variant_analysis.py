"""
Project: SARS-CoV-2 nsp14 Variant Analysis
Author: Dimple Srivastava
"""

from Bio import SeqIO
import pandas as pd

# Load sequences
wuhan_seq = SeqIO.read("data/wuhan_extracted.fasta", "fasta").seq
delta_seq = SeqIO.read("data/delta_extracted.fasta", "fasta").seq
omicron_seq = SeqIO.read("data/omicron_extracted.fasta", "fasta").seq

# Compare mutations
mutations = []

for i in range(min(len(wuhan_seq), len(delta_seq), len(omicron_seq))):
    if wuhan_seq[i] != delta_seq[i] or wuhan_seq[i] != omicron_seq[i]:
        mutations.append([
            i + 1,
            wuhan_seq[i],
            delta_seq[i],
            omicron_seq[i]
        ])

# Create dataframe
df = pd.DataFrame(
    mutations,
    columns=["Position", "Wuhan", "Delta", "Omicron"])

# Save results
df.to_csv("results/nsp14_variant_comparison.csv", index=False)

print("nsp14 variant analysis completed")


