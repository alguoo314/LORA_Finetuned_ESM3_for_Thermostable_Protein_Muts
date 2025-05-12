import torch
import pandas as pd
import random
from test import Experiment  # assumes test (1).py has been renamed to test.py
from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer
from peft import set_peft_model_state_dict

# === 1. Initialize Experiment and load model + test data ===
exp = Experiment()
exp.load_data('/om/user/oliviat/18.S997_project/dpo_triplets.xlsx')
# Load LoRA weights
lora_path = '/om/user/oliviat/18.S997_project/lora_state_dict.pt'
state_dict = torch.load(lora_path, map_location=exp.device)
set_peft_model_state_dict(exp.model, state_dict)
exp.model.eval()

# === 2. Load conserved site info ===
conserved_df = pd.read_excel("/om/user/oliviat/18.S997_project/all_aug_hotmusic_data_conserved_residues_summary.xlsx")
conserved_map = dict(zip(conserved_df['sequence'], conserved_df['conserved residues indices']))

# === 3. Get first 5 unique test sequences ===
test_df = exp.test_loader.dataset.triplet_df
sequences = test_df['sequence'].drop_duplicates().head(100).tolist()
unique_test_seqs = [sequences[i] for i in [0,1,4,6,7,9,10,14,15,17,18]] #single-chain-entries in the testing set
tokenizer = EsmSequenceTokenizer()


results = []
id_to_token = exp.tokenizer.convert_ids_to_tokens
pdbs = test_df['pdb_id'].drop_duplicates().head(100).tolist()
selected = [pdbs[i] for i in [0,1,4,6,7,9,10,14,15,17,18]] #single-chain-entries in the testing set


# === 4. Predict masked positions ===
j=-1
for seq in unique_test_seqs:
    j+=1
    results = []
    conserved = conserved_map.get(seq, "")
    conserved_indices = set(int(i) - 1 for i in str(conserved).split(',') if i.isdigit())

    non_conserved_positions = [i for i in range(len(seq)) if i not in conserved_indices]
    if len(non_conserved_positions) < 5:
        print(f"Skipping {seq} — not enough non-conserved sites.")
        continue
    pdb_id = selected[j].split("|")[0]
    #follow the pdb dataset sequence starting and ending index
    if pdb_id == '1BYW':
       non_conserved_positions= [i for i in non_conserved_positions if i< 135 and i > 24]
    if pdb_id == '1M21':
       non_conserved_positions= [i for i in non_conserved_positions if i< 41+487 and i > 41]


    for i in range(60):
        #print(i)
        pos = random.choice(non_conserved_positions)
        masked_seq = seq[:pos] + tokenizer.mask_token + seq[pos+1:]

        # Tokenize and run model
        tokens = tokenizer(masked_seq, return_tensors='pt', truncation=True)
        input_ids = tokens['input_ids'].to(exp.device)
        with torch.no_grad():
            out = exp.model(sequence_tokens=input_ids)
            logits = out.sequence_logits[0]  # shape: (seq_len, vocab_size)

        mask_idx = input_ids[0].tolist().index(tokenizer.mask_token_id)
        pred_token_id = logits[mask_idx].argmax().item()
        pred_aa  = id_to_token(pred_token_id)

        results.append({
            'pdb_id':pdb_id,
            'sequence': seq,
            'wt':seq[pos],
            'masked_position': pos + 1,  # back to 1-indexed
            'predicted_aa': pred_aa
        })
    pred_df = pd.DataFrame(results)
    pred_df.to_csv(pdb_id+"_masked_position_predictions.csv", index=False)