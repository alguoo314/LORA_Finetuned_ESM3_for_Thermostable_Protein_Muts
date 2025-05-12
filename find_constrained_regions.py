#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd


# In[3]:


#!pip3 install openpyxl


# In[4]:


df = pd.read_excel("ireProt_HotMusic_combined.xlsx")

# Drop rows where 'dTm' is missing
df_filtered = df.dropna(subset=["dTm"])

# Keep only the desired columns
columns_to_keep = [
    "experiment_id", "protein_name", "uniprot_id", "pdb_id", "chain",
    "position", "wild_type", "mutation", "ddG", "dTm", "sequence"
]
df_filtered = df_filtered[columns_to_keep]


# In[5]:


df_filtered.shape


# In[6]:


len(set(df_filtered['uniprot_id'].dropna()))


# In[7]:


from Bio.Blast import NCBIWWW


# In[8]:


sequences = list(set(df_filtered["sequence"].dropna().tolist()))
len(sequences)


# In[9]:


grouped = df_filtered.groupby("uniprot_id")["sequence"].nunique()

pdb_ids_with_diff_sequences = grouped[grouped > 1].index

# Filter the dataframe to show only those rows
df_conflicting = df_filtered[df_filtered["uniprot_id"].isin(pdb_ids_with_diff_sequences)]

# If you just want two such rows (first example)
example = df_conflicting.groupby("uniprot_id").head(2)

print(example)


# In[10]:


grouped


# In[11]:


df_filtered_2 = df_filtered[["protein_name","uniprot_id","pdb_id","chain","sequence"]]
df_filtered_2


# In[12]:


df=df_filtered_2.drop_duplicates()
df


# In[23]:


df=df.reset_index()
df


# In[140]:


import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
import os
from collections import Counter
from pathlib import Path

Entrez.email = "alguo@mit.edu"  # ← required by NCBI

# Prepare output
output_rows = []

# Helper: calculate Shannon entropy per position
def shannon_entropy(column):
    freqs = Counter(column)
    probs = [v / len(column) for v in freqs.values()]
    return -sum(p * np.log2(p) for p in probs if p > 0)

# Filter function to check for mutant keywords in titles
def is_mutant(title):
    mutant_keywords = ['mutant', 'variant', 'mutation', 'substitution', 'mutagenesis']
    return any(keyword in title.lower() for keyword in mutant_keywords)

# Function to calculate sequence identity between two sequences
def calculate_identity(seq1, seq2):
    aligned = 0
    identical = 0
    for a, b in zip(seq1, seq2):
        if a != '-' and b != '-':  # Both are not gaps
            aligned += 1
            if a == b:
                identical += 1
    return identical / aligned if aligned > 0 else 0

# Loop through sequences
for i, row in df.iterrows():
    print(i)
    seq = row["sequence"]
    if not isinstance(seq, str) or len(seq) < 30:
        output_rows.append({**row, "conserved residues indices": None})
        print("failed", i, row["pdb_id"])
        continue
    
    try:
        # Run BLAST with modified parameters
        # Add entrez_query to exclude artificial mutants and use a stricter E-value
        result_handle = NCBIWWW.qblast(
            "blastp", "swissprot", seq, 
            hitlist_size=100,  # Fetch more hits initially for filtering
            expect= 0.001,     # Stricter E-value threshold
            entrez_query="NOT (mutant[Title] OR mutation[Title] OR variant[Title])"
        )
        
        blast_xml = f"blast_{i}.xml"
        with open(blast_xml, "w") as out_handle:
            out_handle.write(result_handle.read())
        
        with open(blast_xml) as handle:
            blast_record = NCBIXML.read(handle)
        
        # Filter alignments to exclude likely mutants based on titles
        filtered_alignments = []
        for alignment in blast_record.alignments:
            if not is_mutant(alignment.title):
                filtered_alignments.append(alignment)
        
        # If we have too few sequences after filtering, relax constraints
        if len(filtered_alignments) < 10:
            print("len(filtered_alignments) < 10")
            with open(blast_xml) as handle:
                blast_record = NCBIXML.read(handle)
            filtered_alignments = blast_record.alignments
            
        
        ids = [hit.accession for hit in filtered_alignments]
        
        # Get sequences
        records = []
        query_record = SeqRecord(Seq(seq), id="query", description="")
        records.append(query_record)  # Add query sequence first
        
        # Download the sequences
        for acc in ids:
            try:
                fetch = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
                record = SeqIO.read(fetch, "fasta")
                records.append(record)
                fetch.close()
            except Exception as e:
                print(f"Failed to fetch {acc}: {e}")
        
        # Filter sequences by identity to ensure diversity (if we have enough sequences)
        if len(records) > 10:
            diverse_records = [records[0]]  # Always keep query
            for record in records[1:]:
                # Only add if not too similar to query (avoid close homologs & mutants)
                seq_identity = calculate_identity(str(records[0].seq), str(record.seq))
                if 0.3 <= seq_identity <= 0.9:  # Keep sequences with 30-90% identity
                    diverse_records.append(record)
                if len(diverse_records) >= 20:  # Limit to 20 sequences for efficiency
                    break
            
            records = diverse_records if len(diverse_records) >= 5 else records
            if len(diverse_records) < 5:
                print("not enough 0.3 <= seq_identity <= 0.9")
                
        
        # Write sequences to FASTA file
        fasta_file = f"seq_{i}.fasta"
        aligned_file = f"aligned_{i}.fasta"
        SeqIO.write(records, fasta_file, "fasta")
        
        # Align
        clustalo_path = "/home/ubuntu/anaconda3/bin/clustalo"
        cline = ClustalOmegaCommandline(
            clustalo_path, 
            infile=fasta_file, 
            outfile=aligned_file, 
            auto=True, 
            verbose=False
        )
        cline()
        
        # Compute conservation scores
        alignment = AlignIO.read(aligned_file, "fasta")
        columns = list(zip(*[str(rec.seq) for rec in alignment]))
        entropies = [shannon_entropy(col) for col in columns]
        max_entropy = np.log2(20)
        conservation_scores = [1 - e / max_entropy for e in entropies]
        
        # Map conserved indices from aligned seq to original seq
        query_seq = str(alignment[0].seq)
        orig_idx = 0
        conserved_positions = []
        for aln_pos, aa in enumerate(query_seq):
            if aa == "-":
                continue  # gap
            if conservation_scores[aln_pos] >= 0.8:
                conserved_positions.append(orig_idx)
            orig_idx += 1
        
        # Save result
        print("len of conserved_positions", len(conserved_positions))
        output_rows.append({
            "protein_name": row["protein_name"],
            "uniprot_id": row["uniprot_id"],
            "pdb_id": row["pdb_id"],
            "sequence": seq,
            "chain": row["chain"],
            "conserved residues indices": conserved_positions
        })
        out_df = pd.DataFrame(output_rows)
        out_df.to_excel(str(i)+"_conserved_residues_summary.xlsx", index=False)
        
        # Cleanup
        # Uncomment these if you want to delete intermediate files
        # os.remove(fasta_file)
        # os.remove(aligned_file)
        # os.remove(blast_xml)
        
    except Exception as e:
        print(f"Row {i} failed: {e}")
        output_rows.append({**row, "conserved residues indices": None})

# Final table
out_df = pd.DataFrame(output_rows)
out_df.to_excel("conserved_residues_summary.xlsx", index=False)
print("Saved to conserved_residues_summary.xlsx")

