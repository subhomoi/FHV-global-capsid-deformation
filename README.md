# FHV Global Capsid Deformation

Python scripts for structural, electrostatic, and inter-subunit interface analysis of the FHV capsid during disassembly.

## Overview

This repository contains the custom Python scripts used for computational analysis of the FHV capsid during early disassembly, including:

- global geometrical deformation analysis
- electrostatic/protonation analysis
- inter-subunit interface analysis
- publication-style figure generation

The scripts are provided to support reproducibility of the computational analyses and associated figures described in the manuscript.

## Repository contents

```text
FHV-global-capsid-deformation/
├── README.md
├── LICENSE
├── requirements.txt
└── scripts/
    ├── FHV_global_capsid_deformation_analysis.py
    ├── FHV_global_deformation_figure.py
    ├── electrostatic_analysis.py
    └── plot_inter_subunit_contact_summary.py
```

### `FHV_global_capsid_deformation_analysis.py`

Performs global geometrical deformation analysis comparing the native and early-disassembly FHV structures.

The analysis calculates:

- global Cα displacement
- radial expansion/contraction
- tangential deformation
- angular reorientation
- global radius distributions
- chain/subunit centroid movement
- chain-chain interface changes
- native contact loss
- residues with the largest deformation

The analysis is geometric and does not calculate mechanical strain energy.

### `FHV_global_deformation_figure.py`

Generates publication-style figures from the output of the global deformation analysis.

The figure-generation script uses the per-residue deformation data to generate radial-expansion and angular-deformation plots in PNG and PDF formats.

### `electrostatic_analysis.py`

Performs PROPKA-based analysis of acidic-residue protonation changes in the FHV structure 4FTB between pH 7.0 and pH 5.5.

The script identifies acidic residues showing a substantial increase in predicted protonation at lower pH.

### `plot_inter_subunit_contact_summary.py`

Generates a publication-style summary figure comparing the BC–HI inter-subunit interface between the native and early-disassembly structures.

**Important:** this is a figure-generation script. It does not independently calculate interface area, salt bridges, hydrogen bonds, or non-bonded contacts from coordinate files. The precomputed values used for the manuscript figure are explicitly defined in the script.

---

## Requirements

Python 3.9+ and the packages listed in `requirements.txt`.

The electrostatic analysis additionally requires PROPKA 3.x with the `propka3` executable available on `PATH`.

Install the Python dependencies using:

```bash
pip install -r requirements.txt
```

---

# Global deformation analysis

## Input structures

Place the following files in the working directory:

- `4ftb_flat2.pdb` — native FHV structure
- `int_flat2.pdb` — early-disassembly structure

The structures are assumed to have already been structurally aligned before analysis.

## Analysis

Run:

```bash
python scripts/FHV_global_capsid_deformation_analysis.py
```

The analysis calculates:

- global Cα displacement
- radial expansion/contraction
- tangential deformation
- angular reorientation
- global radius distributions
- chain/subunit centroid movement
- chain-chain interface changes
- native contact loss
- residues with the largest deformation

The analysis is geometric and does not calculate mechanical strain energy.

The main per-residue output used for figure generation is:

```text
FHV_Global_Deformation_per_residue.csv
```

The analysis also generates summary tables describing the calculated deformation and interface properties.

## Figure generation

After the analysis has produced `FHV_Global_Deformation_per_residue.csv`, run:

```bash
python scripts/FHV_global_deformation_figure.py
```

This generates PNG and PDF versions of the radial-expansion and angular-deformation figure.

---

# Electrostatic / protonation analysis

## Input structure

The electrostatic analysis uses FHV structure:

```text
4FTB
```

The script downloads the structure directly from the RCSB PDB.

## Analysis

Run:

```bash
python scripts/electrostatic_analysis.py
```

The default analysis uses:

- PDB: `4FTB`
- high pH: 7.0
- low pH: 5.5
- minimum net protonation gain: 40 percentage points

The script prepares the structure, runs PROPKA, and calculates the predicted protonation fraction of acidic residues at the two pH values.

The resulting table is:

```text
4FTB_major_neutralization_sites.csv
```

The script reports acidic residues meeting the specified protonation-gain threshold.

---

# BC–HI inter-subunit interface analysis

## Interface comparison

The figure compares the BC–HI inter-subunit interface between:

- native FHV: `4FTB / ugf2`
- early-disassembly intermediate: `9LZL / uge6`

## Figure generation

Run:

```bash
python scripts/plot_inter_subunit_contact_summary.py
```

The script generates:

```text
Nature_Style_BCHI_Summed_BW.png
```

The figure summarizes precomputed measurements for:

- total interface area
- salt bridges
- hydrogen bonds
- non-bonded contacts
- relative reduction between the native and intermediate interfaces

**Important:** the script is a figure-generation script and does not independently calculate these interface properties from structural coordinates. The precomputed values used for the manuscript figure are explicitly defined within the script.

---

# Reproducibility

The global deformation analysis assumes that the input structures have been structurally aligned before geometrical comparison, as described in the manuscript Methods.

Structural coordinates should be obtained from the corresponding public structural database entries cited in the manuscript.

The scripts are provided as the custom computational analysis and figure-generation code used in the study. Numerical parameters, thresholds, and analysis settings are retained in the corresponding scripts.

The computational workflow can be summarized as:

```text
Native and early-disassembly PDB structures
                    │
                    ▼
       Structural alignment
                    │
                    ▼
     Global deformation analysis
                    │
                    ├── Cα displacement
                    ├── Radial deformation
                    ├── Tangential deformation
                    ├── Angular deformation
                    ├── Centroid movement
                    ├── Interface changes
                    └── Contact loss
                    │
                    ▼
FHV_Global_Deformation_per_residue.csv
                    │
                    ▼
        Publication-style figures
```

The electrostatic analysis follows a separate workflow:

```text
4FTB structure
      │
      ▼
Structure preparation
      │
      ▼
     PROPKA
      │
      ▼
Predicted pKa values
      │
      ▼
Protonation at pH 7.0 and pH 5.5
      │
      ▼
Major neutralization sites
```

---

# Code availability

The custom Python scripts used for structural, electrostatic, and inter-subunit interface analyses of the FHV capsid during disassembly are available at:

https://github.com/subhomoi/FHV-global-capsid-deformation

---

# Citation

Please cite the associated FHV capsid-disassembly publication when using these scripts.