# FHV Global Capsid Deformation Analysis

## Purpose
This code package contains Python scripts used to quantify global geometrical
differences between the native FHV capsid structure and an early-disassembly
structure and to generate the associated deformation figure.

## Input structures
Place the following files in the working directory:
- `4ftb_flat2.pdb` — native FHV structure
- `int_flat2.pdb` — early-disassembly structure

The structures are assumed to have already been structurally aligned.

## Analysis
Run:

```bash
python scripts/FHV_global_capsid_deformation_analysis.py
```

The analysis calculates:
- global C-alpha displacement
- radial expansion/contraction
- tangential deformation
- angular reorientation
- global radius distributions
- chain/subunit centroid movement
- chain-chain interface changes
- native contact loss
- residues with the largest deformation

The analysis is geometric and does not calculate mechanical strain energy.

## Figure generation
After the analysis has produced
`FHV_Global_Deformation_per_residue.csv`, run:

```bash
python scripts/FHV_global_deformation_figure.py
```

This generates PNG and PDF versions of the radial-expansion and
angular-deformation figure.

## Dependencies
Recommended Python packages:

```text
numpy
scipy
biopython
pandas
matplotlib
plotly
```

Use the package versions recorded in the computational environment used to
generate the manuscript results.

