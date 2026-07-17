# Molecular Bayesian Optimization with PolyChord Nested Sampling

## Overview

This repository is a practice-oriented project for molecular Bayesian Optimization using a SELFIES-based Variational Autoencoder (VAE), a Gaussian Process surrogate model, and PolyChord nested sampling.

The main notebook is:

- `bo_selfies_vae_polychord copy.ipynb`

The workflow uses the QM9 molecular dataset and focuses on optimizing the `gap` property, the HOMO-LUMO gap. The purpose of the project is not to produce production-ready molecular candidates, but to practice the full process of molecular Bayesian Optimization with PolyChord Nested Sampling.

## Objective

The objective is to build an end-to-end molecular optimization pipeline that can:

- represent molecules in a machine-learning-friendly format;
- learn a continuous latent space for molecular structures;
- train a probabilistic surrogate model for a molecular property;
- use PolyChord nested sampling to explore promising latent-space regions;
- decode optimized latent vectors back into molecules;
- validate generated molecular suggestions.

In the notebook, the final target is to suggest five valid molecular candidates after optimization.

## Process

The notebook follows these main stages.

### 1. Load QM9 Data

The QM9 dataset is loaded from `qm9.csv`. The workflow keeps the molecular SMILES strings and the selected target property, `gap`, then removes rows with missing values.

### 2. Convert SMILES to SELFIES

SMILES strings are converted into SELFIES strings. SELFIES is used because it is more robust for molecular generation: decoded SELFIES strings are more likely to produce syntactically valid molecular structures than raw SMILES generation.

### 3. Tokenize and Encode SELFIES

The SELFIES alphabet is built from the dataset. Each molecule is converted into a fixed-length integer sequence using token indices and padding. These encoded sequences are used as input to the neural network.

### 4. Train a SELFIES-VAE

A Variational Autoencoder is trained to compress molecules into a continuous latent representation and reconstruct the original SELFIES sequence.

The VAE uses:

- token embeddings for SELFIES input;
- a GRU encoder;
- latent mean and log-variance layers;
- the reparameterization trick;
- a GRU decoder;
- reconstruction loss plus KL-divergence loss.

The notebook uses a 32-dimensional latent space. After training, the model is saved to `models/selfies_vae.pt` so it can be reused without retraining.

### 5. Encode Molecules into Latent Space

After training, all QM9 molecules are passed through the VAE encoder. The encoder mean is used as a deterministic latent representation for each molecule.

This produces the latent-space training data used by the Gaussian Process model.

### 6. Prepare the Target Property

The HOMO-LUMO gap values are standardized using `StandardScaler`. This gives the surrogate model a normalized regression target.

### 7. Train a Gaussian Process Surrogate

A BoTorch `SingleTaskGP` model is trained to predict the standardized `gap` value from the VAE latent vectors.

Because exact Gaussian Process training is computationally expensive, the notebook randomly selects up to 3,000 latent vectors for GP training.

### 8. Optimize with PolyChord Nested Sampling

Instead of using a standard BoTorch acquisition optimizer, this notebook uses PolyChord nested sampling to explore the VAE latent space.

PolyChord samples from a 32-dimensional latent-space search region. A prior function maps PolyChord's unit hypercube into the observed bounds of the VAE latent space.

The GP posterior is used to compute a reference-aware acquisition score. The reference is the best target value in the complete current dataset: the highest observed value for maximization or the lowest for minimization.

- UCB: `mean + sqrt(beta) * std`, used for maximization;
- LCB: `mean - sqrt(beta) * std`, used for minimization.

The confidence-bound improvement is `UCB - reference` for maximization and `reference - LCB` for minimization. The GP also estimates the probability that a candidate improves on the reference by at least `MINIMUM_IMPROVEMENT`. These are combined as:

```text
score = confidence-bound improvement + PI_WEIGHT * log(probability of improvement)
```

`MINIMUM_IMPROVEMENT` is expressed in the original target units, while `PI_WEIGHT` controls how strongly candidates with a low probability of improvement are penalized. PolyChord uses this combined score as its log-likelihood, so it progressively explores latent-space regions with higher reference-aware acquisition scores. The highest-scoring samples are retained for molecular decoding.

For the current notebook configuration:

- `OPTIMIZATION_TARGET = "minimize"`;
- `ACQUISITION_KIND = "auto"`;
- `auto` selects LCB for minimization;
- `BETA = 2.0`;
- `MINIMUM_IMPROVEMENT = 0.0`;
- `PI_WEIGHT = 1.0`;
- the top 20 latent candidates are retained before decoding.

PolyChord writes sampling outputs to `chains/polychord_bo/`.

### 9. Decode Candidate Molecules

The best latent vectors found by PolyChord are passed through the VAE decoder. The decoded token sequence is converted from SELFIES back into SMILES.

### 10. Validate Outputs

Generated SMILES are checked with RDKit. Invalid molecules are discarded. Duplicate molecules are also removed. The workflow stops once five unique valid molecular suggestions are collected.

## Methodologies Applied

This project combines several methods commonly used in molecular machine learning and Bayesian Optimization:

- SELFIES molecular representation;
- Variational Autoencoder latent-space learning;
- Gaussian Process regression;
- uncertainty-aware acquisition functions;
- UCB/LCB Bayesian Optimization logic;
- PolyChord nested sampling;
- RDKit molecular validity checking.

The key methodological idea is to perform optimization in a learned continuous molecular latent space rather than directly over discrete molecular strings.

## Output Validation

Outputs are validated through several checks:

- SELFIES decoding converts generated token sequences back into molecular strings;
- SELFIES outputs are converted into SMILES;
- RDKit parses each generated SMILES string using `Chem.MolFromSmiles`;
- molecules that fail RDKit parsing are rejected;
- duplicate valid molecules are removed;
- only unique valid molecules are kept as final suggestions.

The notebook successfully demonstrates the full optimization, decoding, and validation cycle by producing five valid SMILES suggestions.

Example generated outputs from the notebook include:

```text
C#CCN=C1[NH1]C(=O)O1
COOC1CC1C(C)=O
C1CCC=C=CC([NH1])[NH1]1
N1CC1(CCCC)CCC#N
C=C=NOC(=O)CC=O
```

These outputs should be interpreted as practice outputs from an unconstrained optimization workflow. They are formally valid according to the implemented RDKit check, but they are not guaranteed to be chemically feasible, stable, synthesizable, or useful.

## Environment Note

PolyChordLite includes compiled Fortran/C components and is easiest to install in a Unix-style environment. In this project, PolyChord was installed and run through an Ubuntu Python virtual environment using Windows Subsystem for Linux (WSL).

Before running the PolyChord section of the notebook, activate the appropriate WSL/Ubuntu environment with `pypolychord` installed.

## Repository Contents

```text
.
|-- bo_selfies_vae.ipynb
|-- bo_selfies_vae_polychord.ipynb
|-- bo_selfies_vae_polychord copy.ipynb
|-- qm9.csv
|-- models/
|   `-- selfies_vae.pt
`-- chains/
    `-- polychord_bo/
```

## Scope and Limitations

This repository is for practice and learning purposes.

Current limitations include:

- no synthesizability constraints;
- no chemical stability filtering beyond RDKit validity;
- no multi-objective optimization;
- exact GP training uses a subset of the dataset for memory reasons;
- generated molecules may be unusual even when formally valid;
- the workflow is intended to demonstrate methodology rather than produce deployable molecular designs.

## Key Takeaway

This project demonstrates how molecular Bayesian Optimization can be practiced by combining a SELFIES-VAE latent space, a Gaussian Process surrogate model, and PolyChord Nested Sampling. It provides a complete learning workflow from molecular representation to optimized latent-space sampling, decoding, and molecular validity checking.
