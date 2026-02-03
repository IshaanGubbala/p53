# p53CAD: A Generative CAD Platform for Protein Functional Rescue

**Computational Biology & Bioinformatics (CBIO) | ISEF 2026**

---

## 🔬 Overview
**p53**, famously known as the "Guardian of the Genome," is the most frequently mutated protein in human cancer. Mutations like **R175H** and **Y220C** cause the protein to unfold and lose its tumor-suppressive function.

**p53CAD** is an end-to-end computational engineering platform designed to restore p53 function. It employs a novel **Functional Manifold Rescue (FMR)** algorithm to navigate the latent representation of proteins and discover distal mutation pathways that restabilize the structure and recover biochemical activity.

---

## 🚀 Key Innovations

### 1. Latent Manifold Rescue (FMR)
Traditional protein design operates in "sequence space" ($20^{393}$ possible p53 variants), which is geographically sparse and computationally intractable. p53CAD encodes sequences into a **320-dimensional semantic latent space** using the **ESM-2** Transformer model. This space acts as a "biophysical manifold" where protein function is a continuous, differentiable property.

### 2. Functional Gradient Steering
We treat protein design as a **constrained optimization problem** on the manifold. By training a deep neural **Functional Oracle** ($f: \mathcal{Z} \to \mathcal{S}$) on over 10,000 p53 variants from Deep Mutational Scanning (DMS) datasets (Giacomelli et al., 2018), we can compute the gradient of function with respect to the protein's latent state:
$$ \Delta z = \eta \cdot \nabla_z f(z) $$
The platform "steers" the cancer-mutated latent vector ($z_{mut}$) towards regions of high functional recovery, generating novel design candidates in real-time.

### 3. Grassmannian Subspace Validation
To ensure that generated rescues preserve the "functional geometry" of the Wild Type protein, we employ **Grassmannian Analysis**. We treat the self-attention mechanisms of the transformer as $k$-dimensional subspaces in $\mathbb{R}^D$. By computing the geodesic distance on the Grassmann manifold $Gr(k, D)$ using Principal Angles, we can mathematically prove that a design candidate maintains the structural integrity required for DNA binding.

---

## 🛠️ The Platform in Action

### Install & Setup
```bash
pip install -e .
```

### 1. Functional Landscape Learning
Train the Oracle to recognize the "shape" of p53 activity.
```bash
p53cad train --dms data/raw/p53_DMS_Giacomelli_2018.csv
```

### 2. De Novo Design (Live)
Launch the generative engine. Watch the **Functional Score** and **Mutation Detection** evolve live as the manifold walker navigates the latent space.
```bash
p53cad design R175H --samples 50
```

### 3. Structural Validation
Quantify the divergence between the candidate's attention subspace and the Wild-Type reference.
```bash
p53cad analyze
```

### 4. 3D Inspection
Auto-generate high-fidelity PyMol sessions to inspect the spatial configuration of rescue mutations.
```bash
p53cad visualize
# Generates p53_rescue_session.pml
```

---

## 📂 Project Architecture
- **`p53cad/engine/`**: The Core Generative Engine (Latent Walkers, Gradients).
- **`p53cad/analysis/`**: Mathematical Validation (Grassmannian Metrics).
- **`p53cad/viz/`**: Visual Engineering (PyMol Automation).
- **`p53cad/data/`**: Standardized Proteomic Data Pipelines.

---

## 🎓 Science Fair Pitch
p53CAD moves protein engineering from "random search" to "precision manifold navigation." By treating the protein as a point on a differentiable manifold, we can use the same backpropagation techniques that power modern AI to engineer real biological cures.

**Contact**: [Your Name/School]
**Category**: Computational Biology (CBIO)