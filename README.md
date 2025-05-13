This is the ThermoDPO repository for the Spring 2025 6.8711 final project codes written by Alina Guo and Olivia Tang.

The repository includes the following scripts:

form_triplets.py: Constructs training triplets (wild type, more thermostable mutation, less thermostable mutation) using data from FireProtDB and HoTMuSiC.

find_constrained_regions.py: Identifies evolutionarily conserved positions for each sequence using multiple sequence alignment and entropy-based conservation scoring.

train_thermodpo.py: Implements ThermoDPO training, including train/validation/test splitting and computation of the DPO loss.

train_mlp_head: Implements finetuning with MLP head on top of a frozen ESM3. Includes train/validation/test splitting and computation of DPO loss.

eval.py: Evaluates the trained model on 11 test proteins by randomly masking a non-conserved residue and prompting the model to predict a thermostabilizing amino acid. 60 positions are sampled per protein sequence.

check_HoTMuSiC_pred.py: Parses melting temperature change predictions from the HoTMuSiC server and retrieves the predicted values for the model-selected mutations generated in eval.py.

finetuned_LoRA_weights.pt: Weights for LoRA weight matrices inserted during ThermoDPO Training

