#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd


# In[3]:


#!pip3 install openpyxl


# In[4]:


df = pd.read_excel("fireprotdb_results_unique_entries.xlsx")  # ← Replace with your actual filename

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


# In[5]:


len(set(df_filtered['uniprot_id'].dropna()))


# In[6]:


df_filtered["pdb_id"] = df_filtered["pdb_id"].astype(str)


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


# In[ ]:





# In[12]:


df=df_filtered_2.drop_duplicates()
df


# In[13]:


df=df.reset_index()
df


# In[14]:


#average duplicated mutations
dtm_avg = df_filtered.groupby(["sequence", "uniprot_id", "position","mutation"], as_index=False)["dTm"].mean()

df_filtered2 = df_filtered.drop(columns=["dTm"]).drop_duplicates(subset=["sequence", "uniprot_id", "position","mutation"])
df_filtered2 = df_filtered2.merge(dtm_avg, on=["sequence", "uniprot_id", "position","mutation"])
df_filtered2


# In[15]:


duplicates = df_filtered[df_filtered.duplicated(subset=["sequence", "uniprot_id", "position","mutation"], keep=False)]
duplicates


# In[16]:


grouped = df_filtered2.groupby(["pdb_id","chain", "sequence", "position"])



# In[17]:


triplets = []

for _, group in grouped:
    if len(group) < 2:
        continue

    wild_types = group["wild_type"].unique()
    if len(wild_types) != 1:
        print('inconsistent wild types')
        #print(wild_types)
        continue  # skip inconsistent wild types

    wt = wild_types[0]
    pos_examples = group[group["dTm"] > 0]
    neg_examples = group[group["dTm"] < 0]

    # Create all combinations of positive and negative examples
    for _, pos in pos_examples.iterrows():
        for _, neg in neg_examples.iterrows():
            triplets.append({
                "protein_name": pos["protein_name"],
                "uniprot_id": pos["uniprot_id"],
                "pdb_id": pos["pdb_id"],
                "chain": pos["chain"],
                "sequence": pos["sequence"],
                "position": pos["position"],
                "wild_type": wt,
                "positive_example_mutation": pos["mutation"],
                "positive_example_dTm": pos["dTm"],
                "negative_example_mutation": neg["mutation"],
                "negative_example_dTm": neg["dTm"],
            })

# Create the new DataFrame
triplet_df = pd.DataFrame(triplets)


# In[ ]:





# In[18]:


triplet_df


# In[81]:


triplet_df.to_excel("triplets_positive_example_has_to_be_positive_dTm_and_vice_versa.xlsx", index=False)


# In[84]:


len(set(triplet_df['pdb_id']))


# In[85]:


#the version where positive example just needs to have higher dTm than negative example
import pandas as pd
from itertools import combinations

triplets = []

for _, group in grouped:
    if len(group) < 2:
        continue  # Need at least 2 mutations to form a pair

    wild_types = group["wild_type"].unique()
    if len(wild_types) != 1:
        continue

    wt = wild_types[0]
    
    # Generate all mutation pairs
    for idx1, idx2 in combinations(group.index, 2):
        row1 = group.loc[idx1]
        row2 = group.loc[idx2]

        if row1["dTm"] == row2["dTm"]:
            continue  # skip if same dTm, no clear positive/negative

        if row1["dTm"] > row2["dTm"]:
            pos, neg = row1, row2
        else:
            pos, neg = row2, row1

        triplets.append({
            "protein_name": pos["protein_name"],
            "uniprot_id": pos["uniprot_id"],
            "pdb_id": pos["pdb_id"],
            "chain": pos["chain"],
            "sequence": pos["sequence"],
            "position": pos["position"],
            "wild_type": wt,
            "positive_example_mutation": pos["mutation"],
            "positive_example_dTm": pos["dTm"],
            "negative_example_mutation": neg["mutation"],
            "negative_example_dTm": neg["dTm"],
        })

# Create the output DataFrame
triplet_df = pd.DataFrame(triplets)


# In[88]:


triplet_df


# In[89]:


len(set(triplet_df['pdb_id']))


# In[90]:


triplet_df.loc[triplet_df['uniprot_id']=='P00883',]


# In[91]:


triplet_df.to_excel("triplets_positive_example_higher_dTm_than_negative_example.xlsx", index=False)


# In[10]:


#if a group has only one mutation, you create a synthetic triplet using the wild type as the opposite of the mutation (depending on the sign of dTm)
import pandas as pd
from itertools import combinations

triplets = []

for _, group in grouped:
    group = group.dropna(subset=["dTm"])  # Ensure dTm is present

    wild_types = group["wild_type"].unique()
    if len(wild_types) != 1:
        continue  # Skip if inconsistent wild types

    wt = wild_types[0]

    if len(group) >= 2:
        # Generate all mutation pairs
        for idx1, idx2 in combinations(group.index, 2):
            row1 = group.loc[idx1]
            row2 = group.loc[idx2]

            if row1["dTm"] == row2["dTm"]:
                continue  # skip if same dTm

            if row1["dTm"] > row2["dTm"]:
                pos, neg = row1, row2
            else:
                pos, neg = row2, row1

            triplets.append({
                "protein_name": pos["protein_name"],
                "uniprot_id": pos["uniprot_id"],
                "pdb_id": pos["pdb_id"],
                "chain": pos["chain"],
                "sequence": pos["sequence"],
                "position": pos["position"],
                "wild_type": wt,
                "positive_example_mutation": pos["mutation"],
                "positive_example_dTm": pos["dTm"],
                "negative_example_mutation": neg["mutation"],
                "negative_example_dTm": neg["dTm"],
            })

    elif len(group) == 1:
        row = group.iloc[0]
        # Compare mutation to wild type based on dTm
        if row["dTm"] > 0:
            pos_mut, pos_dTm = row["mutation"], row["dTm"]
            neg_mut, neg_dTm = wt, 0.0
        elif row["dTm"] < 0:
            pos_mut, pos_dTm = wt, 0.0
            neg_mut, neg_dTm = row["mutation"], row["dTm"]
        else:
            continue  # dTm == 0, skip

        triplets.append({
            "protein_name": row["protein_name"],
            "uniprot_id": row["uniprot_id"],
            "pdb_id": row["pdb_id"],
            "chain": row["chain"],
            "sequence": row["sequence"],
            "position": row["position"],
            "wild_type": wt,
            "positive_example_mutation": pos_mut,
            "positive_example_dTm": pos_dTm,
            "negative_example_mutation": neg_mut,
            "negative_example_dTm": neg_dTm,
        })

# Create the output DataFrame
triplet_df = pd.DataFrame(triplets)


# In[11]:


triplet_df


# In[12]:


len(set(triplet_df['pdb_id']))


# In[13]:


triplet_df.to_excel("triplets_include_wild_type_as_pos_or_neg_example_if_single_mutation_entry.xlsx", index=False)


# In[19]:


hot_music_df = pd.read_csv('HotMuSic_clean_with_sequences_tab_sep.csv',sep='\t')
hot_music_df


# In[20]:


hot_music_df.columns


# In[192]:


hot_music_df_2=hot_music_df[['PDB_id','Chain','Protein','DEFAULT_SEQUENCE','MUTATED_SEQUENCE', 'ΔTm_exp']]
hot_music_df_2 = hot_music_df_2.dropna(subset=['ΔTm_exp'])
hot_music_df_2 = hot_music_df_2.dropna(subset=['DEFAULT_SEQUENCE'])

hot_music_df_2 = hot_music_df_2.dropna(subset=['MUTATED_SEQUENCE'])



# In[193]:


hot_music_df_2


# In[194]:


def find_mutation_info(row):
    default = row['DEFAULT_SEQUENCE']
    mutated = row['MUTATED_SEQUENCE']
    for i, (a, b) in enumerate(zip(default, mutated)):
        if a != b:
            return pd.Series({'position': i+1, 'wild_type': a, 'mutation': b})
    return pd.Series({'position': None, 'wild_type': None, 'mutation': None})

# Apply to your DataFrame
hot_music_df_2[['position', 'wild_type', 'mutation']] = hot_music_df_2.apply(find_mutation_info, axis=1)


# In[195]:


hot_music_df_2


# In[196]:


overlap_count = len(
    set(hot_music_df_2['PDB_id'].str.upper()) & set(df_filtered2['pdb_id'].str.upper())
)
overlap_count


# In[197]:


len(set(hot_music_df_2['PDB_id'].str.upper()))


# In[198]:


len(set(df_filtered2['pdb_id'].str.upper()))


# In[199]:


hot_music_df_2['PDB_id'] = hot_music_df_2['PDB_id'].str.upper()
hot_music_df_2


# In[207]:


hot_music_df_2=hot_music_df_2.drop('MUTATED_SEQUENCE',axis=1)
hot_music_df_2.columns = ['pdb_id','chain','protein_name','sequence','dTm','position','wild_type','mutation']
hot_music_df_2


# In[44]:


#hot_music_df_2.to_excel('cleaned_up_HotMusic_dataset.xlsx')


# In[208]:


from difflib import SequenceMatcher

modified_rows = []
unmatched_proteins = []

sequence_lookup = {}

# Helper: check if two sequences match with ≤3 mismatches in the middle
def is_fuzzy_match(short_seq, long_seq, max_mismatches=3):
    len_short = len(short_seq)
    for i in range(len(long_seq) - len_short + 1):
        window = long_seq[i:i+len_short]
        if window[0] == short_seq[0] and window[-1] == short_seq[-1]:
            mismatches = sum(1 for a, b in zip(window[1:-1], short_seq[1:-1]) if a != b)
            if mismatches <= max_mismatches:
                return i  # Return start index of match
    return -1

# Build lookup
for _, row in df_filtered2.iterrows():
    full_id = row['pdb_id']
    chain = row['chain']
    sequence = row['sequence']
    
    variants = [full_id]
    if '|' in full_id:
        before, after = full_id.split('|', 1)
        variants.extend([before, after])
    
    for variant in variants:
        sequence_lookup[(variant.strip(), chain)] = sequence

# Process each row in hot_music_df_2
for _, row in hot_music_df_2.iterrows():
    pdb_id = row['pdb_id']
    chain = row['chain']
    hot_music_seq = row['sequence']
    
    new_row = row.copy()
    new_row['original_position'] = row['position']
    
    key = (pdb_id, chain)
    full_seq = sequence_lookup.get(key)

    if full_seq:
        if hot_music_seq in full_seq:
            start_index = full_seq.index(hot_music_seq)
            new_row['position'] = row['position'] + start_index
            new_row['sequence'] = full_seq
            new_row['status'] = 'matched'
        else:
            fuzzy_start = is_fuzzy_match(hot_music_seq, full_seq)
            if fuzzy_start != -1:
                new_row['position'] = row['position'] + fuzzy_start
                new_row['sequence'] = full_seq
                new_row['status'] = 'matched'
            else:
                new_row['status'] = 'alignment failed'
                unmatched_proteins.append(row['pdb_id'])
    else:
        new_row['status'] = 'original'

    modified_rows.append(new_row)

# Final trimmed DataFrame
hot_music_df_2_modified_df = pd.DataFrame(modified_rows)

# Print unmatched proteins
if unmatched_proteins:
    print("Unmatched proteins (alignment failed):")
    for protein in set(unmatched_proteins):
        print(protein)


# In[209]:


from difflib import SequenceMatcher

modified_rows2 = []
unmatched_proteins = []

sequence_lookup = {}

# Helper: check for middle-aligned match with ≤3 mismatches
def is_fuzzy_match(short_seq, long_seq, max_mismatches=3):
    len_short = len(short_seq)
    for i in range(len(long_seq) - len_short + 1):
        window = long_seq[i:i+len_short]
        if window[0] == short_seq[0] and window[-1] == short_seq[-1]:
            mismatches = sum(1 for a, b in zip(window[1:-1], short_seq[1:-1]) if a != b)
            if mismatches <= max_mismatches:
                return i  # Return start index of the match
    return -1

# Build the lookup from hot_music_df_2_modified_df
for _, row in hot_music_df_2_modified_df.iterrows():
    full_id = row['pdb_id']
    chain = row['chain']
    sequence = row['sequence']
    
    variants = [full_id]
    if '|' in full_id:
        before, after = full_id.split('|', 1)
        variants.extend([before, after])
    
    for variant in variants:
        sequence_lookup[(variant.strip(), chain)] = sequence

# Align df_filtered2 sequences to those in the lookup
for _, row in df_filtered2.iterrows():
    pdb_id = row['pdb_id'].split('|')[0]
    chain = row['chain']
    df_filtered2_seq = row['sequence']
    
    new_row = row.copy()
    new_row['original_position'] = row['position']
    
    key = (pdb_id, chain)
    full_seq = sequence_lookup.get(key)
    if '|' in row['pdb_id'] and not full_seq:
        key = (row['pdb_id'].split('|')[1], chain)
        full_seq = sequence_lookup.get(key)

    if full_seq:
        if df_filtered2_seq in full_seq:
            start_index = full_seq.index(df_filtered2_seq)
            new_row['position'] = row['position'] - start_index
            new_row['status'] = 'matched'
            new_row['sequence'] = full_seq
        else:
            fuzzy_start = is_fuzzy_match(df_filtered2_seq, full_seq)
            if fuzzy_start != -1:
                new_row['position'] = row['position'] - fuzzy_start
                new_row['status'] = 'matched'
                new_row['sequence'] = full_seq
            else:
                new_row['status'] = 'alignment failed'
                unmatched_proteins.append(row['pdb_id'])
    else:
        new_row['status'] = 'original'

    modified_rows2.append(new_row)

# Final DataFrame
modified_df_filtered2_seq_df = pd.DataFrame(modified_rows2)

# Print unmatched proteins
if unmatched_proteins:
    print("Unmatched proteins (alignment failed):")
    for protein in set(unmatched_proteins):
        print(protein)


# In[217]:


#last pass
from difflib import SequenceMatcher

modified_rows = []
unmatched_proteins = []

sequence_lookup = {}

# Helper: check if two sequences match with ≤3 mismatches in the middle
def is_fuzzy_match(short_seq, long_seq, max_mismatches=3):
    len_short = len(short_seq)
    for i in range(len(long_seq) - len_short + 1):
        window = long_seq[i:i+len_short]
        if window[0] == short_seq[0] and window[-1] == short_seq[-1]:
            mismatches = sum(1 for a, b in zip(window[1:-1], short_seq[1:-1]) if a != b)
            if mismatches <= max_mismatches:
                return i  # Return start index of match
    return -1

# Build lookup
for _, row in modified_df_filtered2_seq_df.iterrows():
    full_id = row['pdb_id']
    chain = row['chain']
    sequence = row['sequence']
    
    variants = [full_id]
    if '|' in full_id:
        before, after = full_id.split('|', 1)
        variants.extend([before, after])
    
    for variant in variants:
        sequence_lookup[(variant.strip(), chain)] = sequence

# Process each row in hot_music_df_2
for _, row in hot_music_df_2.iterrows():
    pdb_id = row['pdb_id']
    chain = row['chain']
    hot_music_seq = row['sequence']
    
    new_row = row.copy()
    new_row['original_position'] = row['position']
    
    key = (pdb_id, chain)
    full_seq = sequence_lookup.get(key)

    if full_seq:
        if hot_music_seq in full_seq:
            start_index = full_seq.index(hot_music_seq)
            new_row['position'] = row['position'] + start_index
            new_row['sequence'] = full_seq
            new_row['status'] = 'matched'
        else:
            fuzzy_start = is_fuzzy_match(hot_music_seq, full_seq)
            if fuzzy_start != -1:
                new_row['position'] = row['position'] + fuzzy_start
                new_row['sequence'] = full_seq
                new_row['status'] = 'matched'
            else:
                new_row['status'] = 'alignment failed'
                unmatched_proteins.append(row['pdb_id'])
    else:
        new_row['status'] = 'original'

    modified_rows.append(new_row)

# Final trimmed DataFrame
final_hot_music_df_2_modified_df = pd.DataFrame(modified_rows)

# Print unmatched proteins
if unmatched_proteins:
    print("Unmatched proteins (alignment failed):")
    for protein in set(unmatched_proteins):
        print(protein)


# In[ ]:


final_hot_music_df_2_modified_df=final_hot_music_df_2_modified_df.drop('ddG',axis=1)
final_hot_music_df_2_modified_df


# In[261]:


final_hot_music_df_2_modified_df['source']='HotMusic'
modified_df_filtered2_seq_df['source']='FireProtDB'


# In[262]:


all_columns = list(set(modified_df_filtered2_seq_df.columns).union(final_hot_music_df_2_modified_df.columns))

# Reindex both DataFrames to have the same columns in the same order
df1_aligned = modified_df_filtered2_seq_df.reindex(columns=all_columns)
df2_aligned = final_hot_music_df_2_modified_df.reindex(columns=all_columns)


stacked_df = pd.concat([df1_aligned, df2_aligned], ignore_index=True)
stacked_df


# In[263]:


stacked_df=stacked_df[['protein_name', 'pdb_id','uniprot_id','source','chain','position','wild_type', 'original_position',
       'mutation', 'dTm', 'sequence','status', 'experiment_id']]
stacked_df


# In[264]:


# Group by 'sequence' and fill missing values in each group
stacked_df_filled = (
    stacked_df
    .groupby('sequence', group_keys=False)
    .apply(lambda group: group.ffill().bfill())
    .reset_index(drop=True)
)
stacked_df_filled


# In[265]:


stacked_df_filled['pdb_id'].isna().sum()


# In[266]:


# Define a function that replaces shorter pdb_ids with the longest one
def replace_with_longest_pdb_id(group):
    # Ensure 'pdb_id' is a string and fill NaN with an empty string
    group['pdb_id'] = group['pdb_id'].fillna('').astype(str)
    
    # Find the longest pdb_id in the group
    
    longest_pdb_id = group['pdb_id'][group['pdb_id'].str.len().idxmax()]
    
    # Replace all pdb_id values in the group with the longest pdb_id
    group['pdb_id'] = longest_pdb_id
    return group

# Apply the function to each group of rows with identical sequences
stacked_df_filled_updated = stacked_df_filled.groupby('sequence', group_keys=False).apply(replace_with_longest_pdb_id)
stacked_df_filled_updated


# In[267]:


stacked_df_filled_updated_sorted = stacked_df_filled_updated.sort_values(by=['pdb_id', 'chain','position'])

stacked_df_filled_updated_sorted.to_excel('fireProt_HotMusic_combined.xlsx')


# In[213]:


#check good
modified_df_filtered2_seq_df[modified_df_filtered2_seq_df['pdb_id'].str.contains('1MJ5',
                                                                             na=False)]


# In[216]:


modified_df_filtered2_seq_df['sequence'][22][177-1]


# In[210]:


#check good
hot_music_df_2_modified_df[hot_music_df_2_modified_df['pdb_id'].str.contains('1GV5',
                                                                             na=False)]


# In[211]:


ori_seq=hot_music_df_2[hot_music_df_2['pdb_id'].str.contains('1GV5',na=False)]['sequence'][415]
ori_seq


# In[212]:


ori_seq[14-1]


# In[169]:


strr=hot_music_df_2_modified_df[hot_music_df_2_modified_df['pdb_id'].str.contains('1GV5',
                                                                  na=False)]['sequence'][415]
strr


# In[173]:


strr[103-1]


# In[147]:


#check: good
hot_music_df_2_modified_df[hot_music_df_2_modified_df['pdb_id'].str.contains('5PTI', na=False)]


# In[ ]:





# In[140]:


hot_music_df_2_modified_df[hot_music_df_2_modified_df['pdb_id'].str.contains('5PTI', na=False)]['sequence'][1591][46-1]


# In[142]:


modified_df_filtered2_seq_df[modified_df_filtered2_seq_df['pdb_id'].str.contains('1L63', na=False)]


# In[273]:


modified_df_filtered2_seq_df[modified_df_filtered2_seq_df['pdb_id'].str.contains('1CUN', na=False)]


# In[277]:


hot_music_df_2_modified_df[hot_music_df_2_modified_df['pdb_id'].str.contains('1KF', na=False)]


# In[272]:


hot_music_df_2_modified_df[hot_music_df_2_modified_df['pdb_id'].str.contains('1L63', na=False)]


# In[279]:


modified = pd.read_excel('manual_modifed_fireProt_HotMusic_combined.xlsx').iloc[:,1:]
modified


# In[283]:


modified_unique = modified.drop_duplicates(subset=["chain","position", "wild_type", "mutation", "dTm"], keep='first')
modified_unique


# In[292]:


duplicates = modified[modified.duplicated(subset=["position", "wild_type", "mutation", "dTm"], keep=False)]
duplicates_sorted = duplicates.sort_values(by=["pdb_id", "position","mutation"])# Display the duplicated rows
duplicates_sorted


# In[294]:


duplicates_sorted = duplicates.sort_values(by=[ "position","mutation","pdb_id"])# Display the duplicated rows
duplicates_sorted


# In[299]:


dtm_avg = modified_unique.groupby(["sequence", "pdb_id", "position","mutation"], as_index=False)["dTm"].mean()

modified_unique2 = modified_unique.drop(columns=["dTm"]).drop_duplicates(subset=["sequence", "pdb_id", "position","mutation"])
modified_unique2 = modified_unique2.merge(dtm_avg, on=["sequence", "pdb_id", "position","mutation"])
modified_unique2


# In[300]:


modified_unique2.to_excel('final_fireProt_HotMusic_combined.xlsx')


# In[20]:


modified_unique2=pd.read_excel('final_fireProt_HotMusic_combined.xlsx')
modified_unique2


# In[21]:


for _, row in modified_unique2.iterrows():
    seq = row['sequence']
    pos = row['position']
    wt = row['wild_type']

    # Check: if the letter at (position - 1) doesn't match the wild type residue
    if pos <= len(seq) and seq[pos - 1] != wt:
        print(row[['sequence', 'position', 'wild_type']])


# In[22]:


correct= pd.read_excel('all_aug_hotmusic_data_conserved_residues_summary.xlsx')
correct


# In[23]:


# Merge on both 'protein_name' and 'pdb_id'
merged = modified_unique2.merge(
    correct[['protein_name', 'pdb_id', 'sequence']], 
    on=['protein_name', 'pdb_id'], 
    how='left', 
    suffixes=('', '_correct')
)

# Replace 'sequence' with 'sequence_correct' where available
merged['sequence'] = merged['sequence_correct'].combine_first(merged['sequence'])

# Drop the extra column used for replacement
modified_unique3 = merged.drop(columns=['sequence_correct'])
modified_unique3


# In[24]:


for _, row in modified_unique3.iterrows():
    seq = row['sequence']
    pos = row['position']
    wt = row['wild_type']

    # Check: if the letter at (position - 1) doesn't match the wild type residue
    if pos <= len(seq) and seq[pos - 1] != wt:
        print(row[['sequence', 'position', 'wild_type']])
#good


# In[43]:


# Drop rows where pdb_id is NaN
filtered_df = modified_unique3.dropna(subset=['pdb_id'])

# Group and find pdb_ids with conflicting sequences
duplicates = filtered_df.groupby(['pdb_id','chain'])['protein_name'].nunique()
conflicting = duplicates[duplicates > 1]

# Extract the conflicting rows
conflicting_rows = filtered_df[filtered_df['pdb_id'].isin(conflicting.index)]
conflicting_rows


# In[29]:


grouped = modified_unique3.groupby(["pdb_id","chain", "sequence", "position"])



# In[30]:


triplets = []

for _, group in grouped:
    if len(group) < 2:
        continue

    wild_types = group["wild_type"].unique()
    if len(wild_types) != 1:
        print('inconsistent wild types')
        #print(wild_types)
        continue  # skip inconsistent wild types

    wt = wild_types[0]
    pos_examples = group[group["dTm"] > 0]
    neg_examples = group[group["dTm"] < 0]

    # Create all combinations of positive and negative examples
    for _, pos in pos_examples.iterrows():
        for _, neg in neg_examples.iterrows():
            triplets.append({
                "protein_name": pos["protein_name"],
                "uniprot_id": pos["uniprot_id"],
                "pdb_id": pos["pdb_id"],
                "chain": pos["chain"],
                "sequence": pos["sequence"],
                "position": pos["position"],
                "wild_type": wt,
                "positive_example_mutation": pos["mutation"],
                "positive_example_dTm": pos["dTm"],
                "negative_example_mutation": neg["mutation"],
                "negative_example_dTm": neg["dTm"],
            })

# Create the new DataFrame
triplet_df = pd.DataFrame(triplets)


# In[31]:


triplet_df


# In[32]:


duplicates = triplet_df[triplet_df['positive_example_mutation'] == triplet_df['negative_example_mutation']]
print(duplicates)


# In[33]:


triplet_df.to_excel("fireProt_hotmusic_triplets_positive_example_has_to_be_positive_dTm_and_vice_versa.xlsx", index=False)


# In[34]:


len(set(triplet_df['pdb_id']))


# In[35]:


#the version where positive example just needs to have higher dTm than negative example
import pandas as pd
from itertools import combinations

triplets = []

for _, group in grouped:
    if len(group) < 2:
        continue  # Need at least 2 mutations to form a pair

    wild_types = group["wild_type"].unique()
    if len(wild_types) != 1:
        continue

    wt = wild_types[0]
    
    # Generate all mutation pairs
    for idx1, idx2 in combinations(group.index, 2):
        row1 = group.loc[idx1]
        row2 = group.loc[idx2]

        if row1["dTm"] == row2["dTm"]:
            continue  # skip if same dTm, no clear positive/negative

        if row1["dTm"] > row2["dTm"]:
            pos, neg = row1, row2
        else:
            pos, neg = row2, row1

        triplets.append({
            "protein_name": pos["protein_name"],
            "uniprot_id": pos["uniprot_id"],
            "pdb_id": pos["pdb_id"],
            "chain": pos["chain"],
            "sequence": pos["sequence"],
            "position": pos["position"],
            "wild_type": wt,
            "positive_example_mutation": pos["mutation"],
            "positive_example_dTm": pos["dTm"],
            "negative_example_mutation": neg["mutation"],
            "negative_example_dTm": neg["dTm"],
        })

# Create the output DataFrame
triplet_df = pd.DataFrame(triplets)


# In[36]:


triplet_df


# In[37]:


len(set(triplet_df['pdb_id']))


# In[38]:


triplet_df.to_excel("fireProt_hotmusic_triplets_positive_example_higher_dTm_than_negative_example.xlsx", index=False)


# In[39]:


#if a group has only one mutation, you create a synthetic triplet using the wild type as the opposite of the mutation (depending on the sign of dTm)
import pandas as pd
from itertools import combinations

triplets = []

for _, group in grouped:
    group = group.dropna(subset=["dTm"])  # Ensure dTm is present

    wild_types = group["wild_type"].unique()
    if len(wild_types) != 1:
        continue  # Skip if inconsistent wild types

    wt = wild_types[0]

    if len(group) >= 2:
        # Generate all mutation pairs
        for idx1, idx2 in combinations(group.index, 2):
            row1 = group.loc[idx1]
            row2 = group.loc[idx2]

            if row1["dTm"] == row2["dTm"]:
                continue  # skip if same dTm

            if row1["dTm"] > row2["dTm"]:
                pos, neg = row1, row2
            else:
                pos, neg = row2, row1

            triplets.append({
                "protein_name": pos["protein_name"],
                "uniprot_id": pos["uniprot_id"],
                "pdb_id": pos["pdb_id"],
                "chain": pos["chain"],
                "sequence": pos["sequence"],
                "position": pos["position"],
                "wild_type": wt,
                "positive_example_mutation": pos["mutation"],
                "positive_example_dTm": pos["dTm"],
                "negative_example_mutation": neg["mutation"],
                "negative_example_dTm": neg["dTm"],
            })

    elif len(group) == 1:
        row = group.iloc[0]
        # Compare mutation to wild type based on dTm
        if row["dTm"] > 0:
            pos_mut, pos_dTm = row["mutation"], row["dTm"]
            neg_mut, neg_dTm = wt, 0.0
        elif row["dTm"] < 0:
            pos_mut, pos_dTm = wt, 0.0
            neg_mut, neg_dTm = row["mutation"], row["dTm"]
        else:
            continue  # dTm == 0, skip

        triplets.append({
            "protein_name": row["protein_name"],
            "uniprot_id": row["uniprot_id"],
            "pdb_id": row["pdb_id"],
            "chain": row["chain"],
            "sequence": row["sequence"],
            "position": row["position"],
            "wild_type": wt,
            "positive_example_mutation": pos_mut,
            "positive_example_dTm": pos_dTm,
            "negative_example_mutation": neg_mut,
            "negative_example_dTm": neg_dTm,
        })

# Create the output DataFrame
triplet_df = pd.DataFrame(triplets)


# In[40]:


triplet_df


# In[41]:


len(set(triplet_df['pdb_id']))


# In[42]:


triplet_df.to_excel("fireProt_hotmusic_triplets_include_wild_type_as_pos_or_neg_example_if_single_mutation_entry.xlsx", index=False)


# In[ ]:





# In[ ]:





# In[ ]:





# In[48]:


sequence_to_mut_positions = triplet_df.groupby("sequence")["position"].apply(list).to_dict()


# In[45]:


modified_unique2=pd.read_excel('all_aug_hotmusic_data_conserved_residues_summary.xlsx')
modified_unique2


# In[60]:


import ast
sequence_to_conserved = modified_unique2.set_index("sequence")["conserved residues indices"].to_dict()

# convert string representations of lists to actual lists
sequence_to_conserved = {k: ast.literal_eval(v) for k, v in sequence_to_conserved.items()}
sequence_to_conserved = {k: [x + 1 for x in v] for k, v in sequence_to_conserved.items()}

sequence_to_conserved


# In[61]:


def compute_percent_overlap(pos_dict, conserved_dict):
    overlap_percent = {}
    for seq in pos_dict:
        if seq in conserved_dict:
            positions = set(pos_dict[seq])
            conserved = set(conserved_dict[seq])
            if positions:
                overlap = positions & conserved
                percent = len(overlap) / len(positions) * 100
                overlap_percent[seq] = percent
            else:
                overlap_percent[seq] = 0.0
    return overlap_percent
percents=compute_percent_overlap(sequence_to_conserved, sequence_to_mut_positions)
percents


# In[65]:


import matplotlib.pyplot as plt
plt.hist(percents.values())


# In[69]:


import numpy as np
np.mean(list(percents.values()))

