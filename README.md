# FHV Global Capsid Deformation

Python scripts for structural, electrostatic, and inter-subunit interface analysis of the FHV capsid during disassembly.

## Repository contents

- `scripts/FHV_global_capsid_deformation_analysis.py` — global geometrical deformation analysis of the native and early-disassembly FHV structures.
- `scripts/FHV_global_deformation_figure.py` — publication-style figure generation from the deformation-analysis output.
- `scripts/electrostatic_analysis.py` — PROPKA-based analysis of acidic-residue protonation changes in 4FTB between pH 7.0 and pH 5.5.
- `scripts/plot_inter_subunit_contact_summary.py` — publication-style BC–HI inter-subunit interface summary figure using precomputed interface measurements.

## Requirements

Python 3.9+ and the packages listed in `requirements.txt`. The electrostatic analysis additionally requires PROPKA 3.x with the `propka3` executable available on `PATH`.

## Electrostatic analysis

Run:

```bash
python scripts/electrostatic_analysis.py
```

Defaults reproduce the supplied analysis: PDB `4FTB`, pH 7.0 versus pH 5.5, and a minimum net protonation gain of 40 percentage points. The output is `4FTB_major_neutralization_sites.csv`.

## BC–HI interface figure

Run:

```bash
python scripts/plot_inter_subunit_contact_summary.py
```

The figure compares the native 4FTB/ugf2 interface with the early-disassembly 9LZL/uge6 interface.

**Important:** this script is a figure-generation script. It does not independently calculate interface area, salt bridges, hydrogen bonds, or non-bonded contacts from coordinate files; the precomputed values used for the manuscript figure are explicitly defined in the script.

## Reproducibility

The global deformation analysis assumes that the input structures have been aligned before geometrical comparison, as described in the manuscript methods. Structural coordinates should be obtained from the corresponding public structural database entries cited in the manuscript.

## Code availability

The custom Python scripts used for structural, electrostatic, and inter-subunit interface analyses of the FHV capsid during disassembly are available at:

https://github.com/subhomoi/FHV-global-capsid-deformation

## Citation

Please cite the associated FHV capsid-disassembly publication when using these scripts.
