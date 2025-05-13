This is the repository for the Spring 2025 6.8711 final project codes:

The repository includes the following scripts:
form_triplets.py: Constructs training triplets (wild type, more thermostable mutation, less thermostable mutation) using data from FireProtDB and HoTMuSiC.
find_constrained_regions.py: Identifies evolutionarily conserved positions for each sequence using multiple sequence alignment and entropy-based conservation scoring.
training.py: Implements model training, including train/validation/test splitting and computation of the DPO loss.
eval.py: Evaluates the trained model on 11 test proteins by randomly masking a non-conserved residue and prompting the model to predict a thermostabilizing amino acid. 60 positions are sampled per protein sequence.
check_HoTMuSiC_pred.py: Parses melting temperature change predictions from the HoTMuSiC server and retrieves the predicted values for the model-selected mutations generated in eval.py.

