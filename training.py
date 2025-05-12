import os
import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
import numpy as np

from huggingface_hub import login
from esm.models.esm3 import ESM3
from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer
from peft import LoraConfig, get_peft_model, get_peft_model_state_dict
from torch.utils.data import Dataset, DataLoader
from collections import defaultdict


class DPOTripletDataset(Dataset):
    def __init__(self, triplet_df):
        self.triplet_df = triplet_df

    def __len__(self):
        return len(self.triplet_df)

    def __getitem__(self, idx):
        row = self.triplet_df.iloc[idx]
        wild_type = row['sequence']
        position = int(row['position'])
        pos_mutation = row['positive_example_mutation']
        neg_mutation = row['negative_example_mutation']

        pos_sample, neg_sample = list(wild_type), list(wild_type)
        pos_sample[position-1] = pos_mutation
        neg_sample[position-1] = neg_mutation

        return {
            'wild_type': wild_type,
            'position': position - 1, # push to 0-indexing
            'pos_sample': ''.join(pos_sample),
            'neg_sample': ''.join(neg_sample)
        }

class Experiment:
    def __init__(self):
        self.set_seed()
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.beta = 1.0

        login(token='placeholder')

        og_model = ESM3.from_pretrained('esm3_sm_open_v1')
        self.ref_model = og_model.eval().to(self.device).float()
        for p in self.ref_model.parameters():
            p.requires_grad = False

        target_modules = self.generate_target_modules(num_layers=3)
        lora_config = LoraConfig(
            r=4,
            lora_alpha=32,
            target_modules=target_modules
        )

        train_base = ESM3.from_pretrained('esm3_sm_open_v1')
        self.model = get_peft_model(train_base, lora_config).to(self.device).float()

        lora_params = [p for p in self.model.parameters() if p.requires_grad]
        print(f"new Trainable parameters: {sum(p.numel() for p in lora_params)}")

        lora_params_2 = [p for p in self.ref_model.parameters() if p.requires_grad]
        print(f'old trainable parameters: {sum(p.numel() for p in lora_params_2)}')
        self.optimizer = torch.optim.AdamW(lora_params, lr=1e-4)
    

        self.tokenizer = EsmSequenceTokenizer()
        self.batch_size = 16
        self.epochs = 2

    def generate_target_modules(self, num_layers):
        blocks = sorted(47 - i for i in range(num_layers))
        modules = ['attn.layernorm_qkv.1', 'attn.out_proj', 'ffn.1', 'ffn.3']
        return [f'transformer.blocks.{b}.{m}' for b in blocks for m in modules]

    def set_seed(self, seed=42):
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    def protein_greedy_split(self, df):
        prot_counts = df['sequence'].value_counts().to_dict()
        proteins = list(prot_counts.keys())
        N = len(df)
        val_test = int(0.1 * N)
        targets = [N - 2 * val_test, val_test, val_test]

        proteins.sort(key=lambda p: prot_counts[p], reverse=True)
        splits = defaultdict(list)
        current = [0, 0, 0]
        names = ['train', 'val', 'test']

        for p in proteins:
            count = prot_counts[p]
            ratios = [(current[i] + count) / targets[i] for i in range(3)]
            i = min(range(3), key=lambda i: ratios[i])
            splits[names[i]].append(p)
            current[i] += count

        return (
            df[df['sequence'].isin(splits['train'])].reset_index(drop=True),
            df[df['sequence'].isin(splits['val'])].reset_index(drop=True),
            df[df['sequence'].isin(splits['test'])].reset_index(drop=True)
        )

    def load_data(self, path):
        df = pd.read_excel(path)
        train_df, val_df, test_df = self.protein_greedy_split(df)
        self.train_loader = DataLoader(DPOTripletDataset(train_df), batch_size=self.batch_size, shuffle=True)
        self.val_loader = DataLoader(DPOTripletDataset(val_df), batch_size=self.batch_size, shuffle=False)
        self.test_loader = DataLoader(DPOTripletDataset(test_df), batch_size=self.batch_size, shuffle=False)

    def extract_ref_logits(self, seqs):
        self.ref_model.eval()
        tokens = self.tokenizer(seqs, return_tensors='pt', padding=True, truncation=True)
        input_ids = tokens['input_ids'].to(self.device)
        with torch.no_grad():
            out = self.ref_model(sequence_tokens=input_ids)
            logits = out.sequence_logits[..., :29]
        return F.log_softmax(logits, dim=-1)

    def extract_train_logits(self, seqs):
        self.model.train()  
        tokens = self.tokenizer(seqs, return_tensors='pt', padding=True, truncation=True)
        input_ids = tokens['input_ids'].to(self.device)
        out = self.model(sequence_tokens=input_ids)
        logits = out.sequence_logits[..., :29]
        return F.log_softmax(logits, dim=-1)

    def score_sequence(self, logp, seq):
        tokens = self.tokenizer.encode(seq)[1:-1]
        idx = torch.tensor(tokens, dtype=torch.long, device=self.device)
        return logp[0, torch.arange(idx.size(0)), idx].sum() 
    
    def dpo_loss(self, wt, positions, y_pos, y_neg):
        masked_inputs = [seq[:pos] + self.tokenizer.mask_token + seq[pos+1:]
                         for seq, pos in zip(wt, positions)]
        logp_ref = self.extract_ref_logits(masked_inputs)
        logp_train = self.extract_train_logits(masked_inputs)

        losses = []
        for i in range(len(wt)):
            mask_idx = (tokens := self.tokenizer(masked_inputs[i], return_tensors='pt'))['input_ids'][0].tolist().index(self.tokenizer.mask_token_id)
            pos_ref = logp_ref[i, mask_idx, self.tokenizer._get_token_id(y_pos[i][positions[i]])]
            neg_ref = logp_ref[i, mask_idx, self.tokenizer._get_token_id(y_neg[i][positions[i]])]
            pos_train = logp_train[i, mask_idx, self.tokenizer._get_token_id(y_pos[i][positions[i]])]
            neg_train = logp_train[i, mask_idx, self.tokenizer._get_token_id(y_neg[i][positions[i]])]
            losses.append(-F.logsigmoid(self.beta * (pos_train - pos_ref - neg_train + neg_ref)))
        return torch.stack(losses).mean()
    
    def train_and_val(self):
        for epoch in range(self.epochs):
            self.model.train()
            total_loss = 0.0
            for batch in self.train_loader:
                wts, positions, y_pos, y_neg = batch['wild_type'], batch['position'], batch['pos_sample'], batch['neg_sample']
                loss = self.dpo_loss(wts, positions, y_pos, y_neg)
                
                self.optimizer.zero_grad()
                loss.backward()
                
                total_grad_norm = 0.0
                for _, p in self.model.named_parameters():
                    if p.requires_grad and p.grad is not None:
                        total_grad_norm += p.grad.norm().item()
                # print(f'gradient norm after backward(): {total_grad_norm}')
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)

                lora_name, lora_param = next((n, p) for n, p in self.model.named_parameters() if 'lora' in n)
                before = lora_param.detach().clone()
                self.optimizer.step()
                after = lora_param.detach()
                # print(f'grad norm of {lora_name}: {(after - before).norm()}')
                total_loss += loss.item()
            avg_loss = total_loss / len(self.train_loader)
            val_loss = self.eval()
            print(f'epoch {epoch+1}: training loss={avg_loss:.4f}, val loss={val_loss:.4f}')

    def eval(self, test=False):
        self.model.eval()
        loader = self.test_loader if test else self.val_loader
        total_loss = 0.0
        with torch.no_grad():
            for batch in loader:
                total_loss += self.dpo_loss(batch['wild_type'], batch['position'], batch['pos_sample'], batch['neg_sample']).item()
        return total_loss / len(loader)

    def save(self, path):
        sd = get_peft_model_state_dict(self.model)
        torch.save(sd, path)

if __name__ == "__main__":
    exp = Experiment()
    exp.load_data('/om/user/oliviat/18.S997_project/dpo_triplets.xlsx')
    print("Starting training...")
    exp.train_and_val()
    print("Testing...")
    print(f"Test loss: {exp.eval(test=True):.4f}")
    exp.save('/om/user/oliviat/18.S997_project/lora_state_dict.pt')
    print("Done!")
