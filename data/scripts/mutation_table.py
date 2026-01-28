from Bio import SeqIO
import pandas as pd 

wuhan = SeqIO.read("wuhan_extracted.fasta", "fasta").seq
delta = SeqIO.read("delta_extracted.fasta", "fasta").seq
omicron = SeqIO.read("omicron_extracted.fasta", "fasta").seq

mutations = []

for i in range(min(len(wuhan), len(delta), len(omicron))):
    if wuhan[i] != delta[i] or wuhan[i] != omicron[i]:
        mutations.append([
            i + 1,        
            wuhan[i],
            delta[i],
            omicron[i],
            f"{wuhan[i]}→{delta[i]}/{omicron[i]}"])
        
        print(len(wuhan), len(delta), len(omicron))

df = pd.DataFrame(mutations,
    columns=["Position", "Wuhan", "Delta", "Omicron", "Mutation"])

print(df.head())


