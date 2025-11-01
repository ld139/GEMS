# GEMS
<div align="center">
  <img src="./fig/fig1_R1.svg" width="600" height="300" alt="GEMS Overview">
</div>

GEMS is a multimodal framework for enzyme engineering that ensembles evolutionary (MSA-based), structure-informed, and sequence-based models (SaProt, ESM-IF1, MSA Transformer, and GEMME) to prioritize functional variants in a zero-shot setting.

---

## 1) Installation

Prerequisites
- Conda (Mamba recommended)
- Docker (for GEMME)
- Python 3.9+ (installed via the provided conda envs)

Clone the repository
```
git clone https://github.com/ld139/GEMS.git
cd GEMS
```

Create environments
```
conda env create -f SaProt_env.yaml
conda env create -f esmfold.yaml
```

Install GEMME (Docker, recommended)
```
docker pull elodielaine/gemme:gemme
```

Model weights and checkpoints
- SaProt (default): ./SaProt/weights/PLMs/SaProt_650M_AF2
  - Put the SaProt checkpoint in this folder or pass a custom path via --ckpt_path_saprot
- ESM-IF1 (default): ~/.cache/torch/hub/checkpoints/esm_if1_gvp4_t16_142M_UR50.pt
- MSA Transformer (default): ~/.cache/torch/hub/checkpoints/esm_msa1b_t12_100M_UR50S.pt

You can pre-download ESM weights (optional) by running in the esmfold env:
```
python -c "import esm; esm.pretrained.esm_if1_gvp4_t16_142M_UR50()"
python -c "import esm; esm.pretrained.esm_msa1b_t12_100M_UR50S()"
```
Then point --ckpt_path_esmif1 and --ckpt_path_msatransformer to the downloaded .pt files if they differ from the defaults.

---

## 2) Inputs (what you must prepare)

You provide one dataset ID (DMS_id) and the required input files/folders:

Required
- Wild-type FASTA
  - Path: DATASET/<DMS_id>/<DMS_id>.fasta
- Multiple sequence alignment (MSA)
  - Path: MSA/<DMS_id>/
  - Contents: your MSA file(s), e.g., <DMS_id>.a3m (or other format supported by your ESM/GEMME setup)
  - Note: The pipeline will compute MSA weights automatically into MSA_weights/
- Wild-type structure (PDB)
  - Path: DATASET/<DMS_id>/<DMS_id>.pdb

Indexing and site selection
- --offset controls the indexing convention for mutation positions (default 1 = 1-based)
- Choose one mutation scope:
  - --all_sites (saturation at every position), or
  - --sites "10,20,30" (specific positions), or
  - --combinatorial "10,20,30" (includes single and multi-site combos on provided positions)

Environment executables
- --saprot_env: Python in the SaProt env (default path in run.py)
- --esmfold_env: Python in the esmfold env (default path in run.py)
Update these to your local paths when running.

---

## 3) How to run

Basic example (all-site saturation, 1-based indexing)
```
python run.py \
  --DMS_id <YOUR_DMS_ID> \
  --data_dir DATASET \
  --all_sites \
  --offset 1 \
  --saprot_env /path/to/miniconda3/envs/SaProt/bin/python \
  --esmfold_env /path/to/miniconda3/envs/esmfold/bin/python
```

Target specific sites
```
python run.py \
  --DMS_id <YOUR_DMS_ID> \
  --data_dir DATASET \
  --sites "10,20,30" \
  --offset 1
```

Combinatorial saturation on specified sites
```
python run.py \
  --DMS_id <YOUR_DMS_ID> \
  --data_dir DATASET \
  --combinatorial "10,20,30" \
  --offset 1
```

Override model checkpoints (if not using defaults)
```
python run.py \
  --DMS_id <YOUR_DMS_ID> \
  --ckpt_path_saprot /models/SaProt_650M_AF2 \
  --ckpt_path_esmif1 /models/esm_if1_gvp4_t16_142M_UR50.pt \
  --ckpt_path_msatransformer /models/esm_msa1b_t12_100M_UR50S.pt
```

Notes
- The script will automatically:
  - Generate mutations (based on --all_sites/--sites/--combinatorial)
  - Compute MSA weights (into ./MSA_weights)
  - Run SaProt, ESM-IF1, MSA Transformer, GEMME
  - Ensemble and rank predictions
- Steps are skipped if the corresponding output folder already exists.

---

## 4) Outputs (what you get)

After a successful run, you should see (non-exhaustive):
- DATASET/<DMS_id>/SaProt/ — SaProt scores
- DATASET/<DMS_id>/ESM-IF1/ — ESM-IF1 scores
- DATASET/<DMS_id>/MSA_Transformer/ — MSA Transformer scores
- DATASET/<DMS_id>/GEMME/ — GEMME outputs (invoked via run_gemme.sh)
- MSA_weights/ — MSA sequence weights (auto-generated)
- Final ranking/ensemble results saved by ranking.py under DATASET/<DMS_id>/ (see that script for filenames)


## 5) Command-line reference (run.py)

Required
- --DMS_id: Dataset ID used to locate inputs and write outputs

Common options
- --data_dir: Root of datasets (default: DATASET)
- --offset: Position indexing offset (default: 1)
- Mutation scope (choose one):
  - --all_sites
  - --sites "p1,p2,..."
  - --combinatorial "p1,p2,..."
- Environments:
  - --saprot_env: Python in SaProt env
  - --esmfold_env: Python in esmfold env
- Checkpoints:
  - --ckpt_path_saprot: SaProt checkpoint dir
  - --ckpt_path_esmif1: ESM-IF1 checkpoint file
  - --ckpt_path_msatransformer: MSA Transformer checkpoint file

---

## 6) Troubleshooting

- Docker/GEMME not found
  - Install Docker; pull the image: docker pull elodielaine/gemme:gemme
  - Ensure run_gemme.sh exists, is executable, and volume paths are correct
- Missing inputs
  - DATASET/<DMS_id>/<DMS_id>.fasta must exist
  - MSA/<DMS_id>/ must contain your MSA file(s)
- Checkpoints not found
  - Supply explicit --ckpt_path_* arguments
  - Or pre-download weights into ~/.cache/torch/hub/checkpoints
- Environment executables invalid
  - Point --saprot_env and --esmfold_env to your local conda env Python paths
- CWD-sensitive paths
  - MSA and MSA_weights are expected at repo root; weights are computed with cwd=./ESM internally

---

## 7) Citation

If you use GEMS in your research, please cite this repository and the upstream models:
- [SaProt](https://github.com/westlake-repl/SaProt)
- [ESM-IF1](https://github.com/facebookresearch/esm)
- [MSA Transformer](https://github.com/facebookresearch/esm)
- [GEMME](https://www.lcqb.upmc.fr/GEMME/Home.html)
