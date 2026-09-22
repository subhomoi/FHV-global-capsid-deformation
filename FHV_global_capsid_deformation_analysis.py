#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================
FHV GLOBAL CAPS ID DEFORMATION ANALYSIS
===============================================================

Native:
    4ftb_flat2.pdb

Early disassembly:
    int_flat2.pdb

Purpose
-------
Global geometric comparison of the native FHV capsid and the
early-disassembly structure.

This analysis is designed to address the reviewer request for
the geometrical shift between the native and disassembly states
and its consequence for the overall capsid architecture.

It quantifies:

1. Global C-alpha displacement
2. Radial expansion/contraction
3. Tangential deformation
4. Angular reorientation
5. Global radius distribution
6. Chain/subunit centroid movement
7. Changes in all detected chain-chain interfaces
8. Native contact loss across the capsid
9. Residues showing the largest deformation
10. Spatial distribution of deformation

IMPORTANT
---------
The structures are assumed to have already been structurally
aligned in Chimera.

The script removes only the residual difference between the
global C-alpha centroids before calculating displacement. This
prevents an arbitrary whole-structure translation from appearing
as capsid deformation.

This is a GEOMETRIC analysis. It does not calculate mechanical
strain energy.

Recommended terminology:
    global capsid deformation
    radial expansion
    tangential deformation
    angular reorientation
    interface opening
    contact loss
    deformation consistent with torsional strain

Avoid claiming that this script directly measures mechanical
torsional strain energy.

===============================================================
DEPENDENCIES
===============================================================

    pip install numpy scipy biopython pandas plotly

===============================================================
RUN
===============================================================

Place this script with:

    4ftb_flat2.pdb
    int_flat2.pdb

Then:

    python global_capsid_deformation_analysis.py

===============================================================
"""

import os
import sys
import numpy as np
import pandas as pd

from Bio.PDB import PDBParser
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation


# =============================================================
# USER SETTINGS
# =============================================================

NATIVE_PDB = "4ftb_flat2.pdb"
DISASSEMBLY_PDB = "int_flat2.pdb"

OUTPUT_PREFIX = "FHV_Global_Deformation"

# Heavy atom contact definition
CONTACT_CUTOFF = 4.5

# Minimum contacts for a chain-chain interface
MIN_INTERFACE_CONTACTS = 20

# Top residues in deformation ranking
N_TOP_RESIDUES = 30

# Number of chain interfaces to report
N_TOP_INTERFACES = 20

# HTML visualization
MAKE_HTML = True

# Use all common chains
USE_ALL_COMMON_CHAINS = True


# =============================================================
# STRUCTURE LOADING
# =============================================================

def load_structure(filename):
    if not os.path.isfile(filename):
        raise FileNotFoundError(
            f"\nERROR: file not found:\n    {filename}\n"
        )

    parser = PDBParser(QUIET=True)

    return parser.get_structure(
        os.path.basename(filename),
        filename
    )


def get_first_model(structure):
    return next(structure.get_models())


def get_chains(structure):
    model = get_first_model(structure)

    return {
        chain.id: chain
        for chain in model
    }


# =============================================================
# ATOM EXTRACTION
# =============================================================

def get_ca_atoms(chain):
    """
    Return C-alpha atoms indexed by residue number,
    insertion code and residue name.
    """

    result = {}

    for residue in chain:

        if residue.id[0] != " ":
            continue

        if "CA" not in residue:
            continue

        key = (
            int(residue.id[1]),
            residue.id[2],
            residue.resname.strip()
        )

        result[key] = residue["CA"]

    return result


def get_all_ca_atoms(chains):
    """
    Return all C-alpha atoms with chain/residue identifiers.
    """

    records = []

    for chain_id, chain in chains.items():

        ca_atoms = get_ca_atoms(chain)

        for key, atom in ca_atoms.items():

            resnum, insertion, resname = key

            records.append(
                {
                    "chain": chain_id,
                    "resnum": resnum,
                    "icode": insertion,
                    "resname": resname,
                    "atom": atom
                }
            )

    return records


def get_heavy_atoms(chain):
    atoms = []

    for residue in chain:

        if residue.id[0] != " ":
            continue

        for atom in residue:

            element = atom.element

            if element:

                if element.upper() == "H":
                    continue

            else:

                name = atom.get_name().strip().upper()

                if name.startswith("H"):
                    continue

            atoms.append(atom)

    return atoms


# =============================================================
# MATCH COMMON C-alpha ATOMS
# =============================================================

def match_global_ca(native_chains, target_chains):
    """
    Match C-alpha atoms using:

        chain
        residue number
        insertion code
        residue name
    """

    native_records = get_all_ca_atoms(
        native_chains
    )

    target_records = get_all_ca_atoms(
        target_chains
    )

    target_lookup = {}

    for record in target_records:

        key = (
            record["chain"],
            record["resnum"],
            record["icode"],
            record["resname"]
        )

        target_lookup[key] = record["atom"]

    matched = []

    for record in native_records:

        key = (
            record["chain"],
            record["resnum"],
            record["icode"],
            record["resname"]
        )

        if key not in target_lookup:
            continue

        matched.append(
            {
                "chain":
                    record["chain"],

                "resnum":
                    record["resnum"],

                "icode":
                    record["icode"],

                "resname":
                    record["resname"],

                "native_atom":
                    record["atom"],

                "target_atom":
                    target_lookup[key]
            }
        )

    return matched


# =============================================================
# CENTER STRUCTURES
# =============================================================

def calculate_global_center(matched):
    """
    Mean C-alpha position.
    """

    coords = np.array(
        [
            x["native_atom"].coord
            for x in matched
        ],
        dtype=float
    )

    return coords.mean(axis=0)


def recenter_coordinates(matched):
    """
    Recenter both structures independently on their global
    matched-C-alpha centroid.

    Because the structures were already aligned, this removes
    residual rigid-body translation without performing another
    structural fit.
    """

    native_coords = np.array(
        [
            x["native_atom"].coord
            for x in matched
        ],
        dtype=float
    )

    target_coords = np.array(
        [
            x["target_atom"].coord
            for x in matched
        ],
        dtype=float
    )

    native_center = native_coords.mean(axis=0)
    target_center = target_coords.mean(axis=0)

    native_centered = (
        native_coords -
        native_center
    )

    target_centered = (
        target_coords -
        target_center
    )

    return (
        native_coords,
        target_coords,
        native_center,
        target_center,
        native_centered,
        target_centered
    )


# =============================================================
# GLOBAL DEFORMATION METRICS
# =============================================================

def calculate_global_deformation(matched):
    """
    Calculate displacement, radial and tangential components,
    and angular reorientation for every matched C-alpha atom.
    """

    (
        native_coords,
        target_coords,
        native_center,
        target_center,
        native_centered,
        target_centered
    ) = recenter_coordinates(matched)

    displacement = (
        target_centered -
        native_centered
    )

    rows = []

    for i, record in enumerate(matched):

        r_native = native_centered[i]
        r_target = target_centered[i]
        d = displacement[i]

        r_native_mag = np.linalg.norm(
            r_native
        )

        r_target_mag = np.linalg.norm(
            r_target
        )

        if r_native_mag > 1e-12:

            radial_unit = (
                r_native /
                r_native_mag
            )

            radial_component = float(
                np.dot(
                    d,
                    radial_unit
                )
            )

        else:

            radial_component = 0.0

        displacement_magnitude = float(
            np.linalg.norm(d)
        )

        tangential_component = float(
            np.sqrt(
                max(
                    displacement_magnitude ** 2 -
                    radial_component ** 2,
                    0.0
                )
            )
        )

        radial_change = (
            r_target_mag -
            r_native_mag
        )

        if (
            r_native_mag > 1e-12 and
            r_target_mag > 1e-12
        ):

            cosine = np.dot(
                r_native,
                r_target
            ) / (
                r_native_mag *
                r_target_mag
            )

            cosine = np.clip(
                cosine,
                -1.0,
                1.0
            )

            angular_change = float(
                np.degrees(
                    np.arccos(cosine)
                )
            )

        else:

            angular_change = np.nan

        rows.append(
            {
                "chain":
                    record["chain"],

                "resnum":
                    record["resnum"],

                "icode":
                    record["icode"],

                "resname":
                    record["resname"],

                "native_x":
                    native_centered[i, 0],

                "native_y":
                    native_centered[i, 1],

                "native_z":
                    native_centered[i, 2],

                "target_x":
                    target_centered[i, 0],

                "target_y":
                    target_centered[i, 1],

                "target_z":
                    target_centered[i, 2],

                "dx":
                    d[0],

                "dy":
                    d[1],

                "dz":
                    d[2],

                "displacement_A":
                    displacement_magnitude,

                "radial_change_A":
                    radial_change,

                "radial_component_A":
                    radial_component,

                "tangential_component_A":
                    tangential_component,

                "angular_change_deg":
                    angular_change,

                "native_radius_A":
                    r_native_mag,

                "target_radius_A":
                    r_target_mag
            }
        )

    return (
        pd.DataFrame(rows),
        native_center,
        target_center
    )


# =============================================================
# GLOBAL SUMMARY
# =============================================================

def calculate_summary(df):
    """
    Calculate global deformation statistics.
    """

    displacement = df[
        "displacement_A"
    ].to_numpy()

    radial_change = df[
        "radial_change_A"
    ].to_numpy()

    radial_component = df[
        "radial_component_A"
    ].to_numpy()

    tangential = df[
        "tangential_component_A"
    ].to_numpy()

    angular = df[
        "angular_change_deg"
    ].dropna().to_numpy()

    native_radius = df[
        "native_radius_A"
    ].to_numpy()

    target_radius = df[
        "target_radius_A"
    ].to_numpy()

    summary = {
        "matched_CA":
            len(df),

        "mean_displacement_A":
            float(np.mean(displacement)),

        "median_displacement_A":
            float(np.median(displacement)),

        "std_displacement_A":
            float(np.std(displacement)),

        "max_displacement_A":
            float(np.max(displacement)),

        "mean_radial_change_A":
            float(np.mean(radial_change)),

        "median_radial_change_A":
            float(np.median(radial_change)),

        "max_radial_expansion_A":
            float(np.max(radial_change)),

        "min_radial_change_A":
            float(np.min(radial_change)),

        "mean_radial_component_A":
            float(np.mean(radial_component)),

        "mean_tangential_component_A":
            float(np.mean(tangential)),

        "max_tangential_component_A":
            float(np.max(tangential)),

        "mean_angular_change_deg":
            float(np.mean(angular)),

        "median_angular_change_deg":
            float(np.median(angular)),

        "max_angular_change_deg":
            float(np.max(angular)),

        "native_mean_radius_A":
            float(np.mean(native_radius)),

        "target_mean_radius_A":
            float(np.mean(target_radius)),

        "mean_radius_change_A":
            float(
                np.mean(target_radius) -
                np.mean(native_radius)
            ),

        "native_median_radius_A":
            float(np.median(native_radius)),

        "target_median_radius_A":
            float(np.median(target_radius)),

        "median_radius_change_A":
            float(
                np.median(target_radius) -
                np.median(native_radius)
            )
    }

    # Fraction of residues moving outward
    summary["fraction_outward_radial"] = float(
        np.mean(radial_change > 0.0) * 100.0
    )

    summary["fraction_inward_radial"] = float(
        np.mean(radial_change < 0.0) * 100.0
    )

    # Fraction with substantial deformation
    summary["fraction_displacement_ge_2A"] = float(
        np.mean(displacement >= 2.0) * 100.0
    )

    summary["fraction_displacement_ge_3A"] = float(
        np.mean(displacement >= 3.0) * 100.0
    )

    return summary


# =============================================================
# CHAIN-LEVEL DEFORMATION
# =============================================================

def calculate_chain_statistics(df):
    """
    Summarize deformation for each capsid protein chain.
    """

    records = []

    for chain_id, group in df.groupby(
        "chain",
        sort=True
    ):

        records.append(
            {
                "chain":
                    chain_id,

                "n_CA":
                    len(group),

                "mean_displacement_A":
                    group[
                        "displacement_A"
                    ].mean(),

                "max_displacement_A":
                    group[
                        "displacement_A"
                    ].max(),

                "mean_radial_change_A":
                    group[
                        "radial_change_A"
                    ].mean(),

                "max_radial_change_A":
                    group[
                        "radial_change_A"
                    ].max(),

                "mean_tangential_A":
                    group[
                        "tangential_component_A"
                    ].mean(),

                "mean_angular_change_deg":
                    group[
                        "angular_change_deg"
                    ].mean(),

                "fraction_outward_percent":
                    np.mean(
                        group[
                            "radial_change_A"
                        ] > 0
                    ) * 100.0
            }
        )

    return pd.DataFrame(records)


# =============================================================
# CHAIN CENTROIDS
# =============================================================

def calculate_chain_centroids(
    native_chains,
    target_chains
):
    """
    Calculate C-alpha centroid movement for each common chain.
    """

    rows = []

    common = sorted(
        set(native_chains.keys()) &
        set(target_chains.keys())
    )

    for chain_id in common:

        P, Q, keys = match_chain_ca(
            native_chains[chain_id],
            target_chains[chain_id]
        )

        if P is None:
            continue

        native_centroid = P.mean(axis=0)
        target_centroid = Q.mean(axis=0)

        shift = target_centroid - native_centroid

        rows.append(
            {
                "chain":
                    chain_id,

                "native_x":
                    native_centroid[0],

                "native_y":
                    native_centroid[1],

                "native_z":
                    native_centroid[2],

                "target_x":
                    target_centroid[0],

                "target_y":
                    target_centroid[1],

                "target_z":
                    target_centroid[2],

                "shift_x":
                    shift[0],

                "shift_y":
                    shift[1],

                "shift_z":
                    shift[2],

                "centroid_shift_A":
                    np.linalg.norm(shift),

                "matched_CA":
                    len(keys)
            }
        )

    return pd.DataFrame(rows)


def match_chain_ca(chain_a, chain_b):
    """
    Match C-alpha atoms between two chains.
    """

    ca_a = get_ca_atoms(chain_a)
    ca_b = get_ca_atoms(chain_b)

    common = sorted(
        set(ca_a.keys()) &
        set(ca_b.keys())
    )

    if len(common) < 3:
        return None, None, []

    P = np.array(
        [ca_a[k].coord for k in common],
        dtype=float
    )

    Q = np.array(
        [ca_b[k].coord for k in common],
        dtype=float
    )

    return P, Q, common


# =============================================================
# HEAVY-ATOM INTERFACE ANALYSIS
# =============================================================

def detect_interfaces(chains):
    """
    Detect every chain-chain interface in the native structure.
    """

    chain_ids = sorted(chains.keys())

    results = []

    for i in range(len(chain_ids)):

        for j in range(i + 1, len(chain_ids)):

            chain_a_id = chain_ids[i]
            chain_b_id = chain_ids[j]

            chain_a = chains[chain_a_id]
            chain_b = chains[chain_b_id]

            atoms_a = get_heavy_atoms(chain_a)
            atoms_b = get_heavy_atoms(chain_b)

            if not atoms_a or not atoms_b:
                continue

            coords_a = np.array(
                [x.coord for x in atoms_a],
                dtype=float
            )

            coords_b = np.array(
                [x.coord for x in atoms_b],
                dtype=float
            )

            tree = cKDTree(coords_b)

            neighbors = tree.query_ball_point(
                coords_a,
                r=CONTACT_CUTOFF
            )

            contacts = sum(
                len(x)
                for x in neighbors
            )

            if contacts < MIN_INTERFACE_CONTACTS:
                continue

            P, Q, keys = match_chain_ca(
                chain_a,
                chain_b
            )

            if P is None:
                continue

            centroid_a = coords_a.mean(axis=0)
            centroid_b = coords_b.mean(axis=0)

            centroid_distance = np.linalg.norm(
                centroid_a -
                centroid_b
            )

            results.append(
                {
                    "chain_A":
                        chain_a_id,

                    "chain_B":
                        chain_b_id,

                    "native_contacts":
                        contacts,

                    "native_centroid_distance_A":
                        centroid_distance,

                    "matched_CA":
                        len(keys)
                }
            )

    return results


# =============================================================
# INTERFACE CONTACT LOSS
# =============================================================

def build_target_atom_lookup(chain):
    """
    Lookup for target heavy atoms.
    """

    lookup = {}

    for residue in chain:

        if residue.id[0] != " ":
            continue

        for atom in residue:

            element = atom.element

            if element:

                if element.upper() == "H":
                    continue

            else:

                if atom.get_name().strip().upper().startswith("H"):
                    continue

            key = (
                int(residue.id[1]),
                residue.id[2],
                residue.resname.strip(),
                atom.get_name().strip()
            )

            lookup[key] = atom

    return lookup


def interface_contact_change(
    native_A,
    native_B,
    target_A,
    target_B
):
    """
    Compare native contacts with their corresponding distances
    in the target structure.
    """

    atoms_A = get_heavy_atoms(native_A)
    atoms_B = get_heavy_atoms(native_B)

    coords_B = np.array(
        [x.coord for x in atoms_B],
        dtype=float
    )

    tree = cKDTree(coords_B)

    lookup_A = build_target_atom_lookup(
        target_A
    )

    lookup_B = build_target_atom_lookup(
        target_B
    )

    native_contacts = 0
    retained = 0
    lost = 0

    distance_changes = []

    for atom_A in atoms_A:

        neighbors = tree.query_ball_point(
            atom_A.coord,
            r=CONTACT_CUTOFF
        )

        for index in neighbors:

            atom_B = atoms_B[index]

            residue_A = atom_A.get_parent()
            residue_B = atom_B.get_parent()

            key_A = (
                int(residue_A.id[1]),
                residue_A.id[2],
                residue_A.resname.strip(),
                atom_A.get_name().strip()
            )

            key_B = (
                int(residue_B.id[1]),
                residue_B.id[2],
                residue_B.resname.strip(),
                atom_B.get_name().strip()
            )

            if (
                key_A not in lookup_A or
                key_B not in lookup_B
            ):
                continue

            target_atom_A = lookup_A[key_A]
            target_atom_B = lookup_B[key_B]

            d_native = np.linalg.norm(
                atom_A.coord -
                atom_B.coord
            )

            d_target = np.linalg.norm(
                target_atom_A.coord -
                target_atom_B.coord
            )

            native_contacts += 1

            if d_target <= CONTACT_CUTOFF:
                retained += 1
            else:
                lost += 1

            distance_changes.append(
                d_target - d_native
            )

    if native_contacts == 0:

        return {
            "native_contacts": 0,
            "retained_contacts": 0,
            "lost_contacts": 0,
            "contact_loss_percent": np.nan,
            "mean_distance_change_A": np.nan,
            "median_distance_change_A": np.nan,
            "max_distance_change_A": np.nan
        }

    distance_changes = np.array(
        distance_changes,
        dtype=float
    )

    return {
        "native_contacts":
            native_contacts,

        "retained_contacts":
            retained,

        "lost_contacts":
            lost,

        "contact_loss_percent":
            100.0 * lost / native_contacts,

        "mean_distance_change_A":
            float(np.mean(distance_changes)),

        "median_distance_change_A":
            float(np.median(distance_changes)),

        "max_distance_change_A":
            float(np.max(distance_changes))
    }


def analyze_all_interfaces(
    native_chains,
    target_chains
):
    """
    Analyze every detected native interface.
    """

    detected = detect_interfaces(
        native_chains
    )

    rows = []

    for item in detected:

        A = item["chain_A"]
        B = item["chain_B"]

        if (
            A not in target_chains or
            B not in target_chains
        ):
            continue

        changes = interface_contact_change(
            native_chains[A],
            native_chains[B],
            target_chains[A],
            target_chains[B]
        )

        P_A, Q_A, _ = match_chain_ca(
            native_chains[A],
            target_chains[A]
        )

        P_B, Q_B, _ = match_chain_ca(
            native_chains[B],
            target_chains[B]
        )

        if P_A is None or P_B is None:
            continue

        native_centroid_A = P_A.mean(axis=0)
        native_centroid_B = P_B.mean(axis=0)

        target_centroid_A = Q_A.mean(axis=0)
        target_centroid_B = Q_B.mean(axis=0)

        native_sep = np.linalg.norm(
            native_centroid_A -
            native_centroid_B
        )

        target_sep = np.linalg.norm(
            target_centroid_A -
            target_centroid_B
        )

        rows.append(
            {
                "chain_A":
                    A,

                "chain_B":
                    B,

                "native_contacts":
                    changes["native_contacts"],

                "retained_contacts":
                    changes["retained_contacts"],

                "lost_contacts":
                    changes["lost_contacts"],

                "contact_loss_percent":
                    changes["contact_loss_percent"],

                "native_centroid_separation_A":
                    native_sep,

                "target_centroid_separation_A":
                    target_sep,

                "centroid_separation_change_A":
                    target_sep - native_sep,

                "mean_contact_distance_change_A":
                    changes["mean_distance_change_A"],

                "median_contact_distance_change_A":
                    changes["median_distance_change_A"],

                "maximum_contact_distance_change_A":
                    changes["max_distance_change_A"]
            }
        )

    return pd.DataFrame(rows)


# =============================================================
# TOP DEFORMATION RESIDUES
# =============================================================

def save_top_residues(df):
    """
    Save residues ranked by global displacement.
    """

    ranked = df.sort_values(
        "displacement_A",
        ascending=False
    ).copy()

    ranked["residue"] = (
        ranked["resname"] +
        ranked["resnum"].astype(str)
    )

    columns = [
        "chain",
        "residue",
        "resnum",
        "displacement_A",
        "radial_change_A",
        "radial_component_A",
        "tangential_component_A",
        "angular_change_deg",
        "native_radius_A",
        "target_radius_A"
    ]

    filename = (
        OUTPUT_PREFIX +
        "_top_deformed_residues.csv"
    )

    ranked[
        columns
    ].head(
        N_TOP_RESIDUES
    ).to_csv(
        filename,
        index=False
    )

    return ranked


# =============================================================
# SAVE GLOBAL SUMMARY
# =============================================================

def save_summary(
    summary,
    native_center,
    target_center
):

    summary = dict(summary)

    summary["native_center_x_A"] = native_center[0]
    summary["native_center_y_A"] = native_center[1]
    summary["native_center_z_A"] = native_center[2]

    summary["target_center_x_A"] = target_center[0]
    summary["target_center_y_A"] = target_center[1]
    summary["target_center_z_A"] = target_center[2]

    filename = (
        OUTPUT_PREFIX +
        "_summary.csv"
    )

    pd.DataFrame(
        [summary]
    ).to_csv(
        filename,
        index=False
    )

    return filename


# =============================================================
# PRINT GLOBAL RESULTS
# =============================================================

def print_global_results(summary):

    print("\n")
    print("=" * 80)
    print("GLOBAL CAPS ID DEFORMATION")
    print("=" * 80)

    print(
        f"\nMatched Cα atoms               : "
        f"{summary['matched_CA']}"
    )

    print("\nGLOBAL DISPLACEMENT")

    print(
        f"Mean displacement              : "
        f"{summary['mean_displacement_A']:.3f} Å"
    )

    print(
        f"Median displacement            : "
        f"{summary['median_displacement_A']:.3f} Å"
    )

    print(
        f"Maximum displacement            : "
        f"{summary['max_displacement_A']:.3f} Å"
    )

    print(
        f"Residues ≥2 Å displacement     : "
        f"{summary['fraction_displacement_ge_2A']:.2f}%"
    )

    print(
        f"Residues ≥3 Å displacement     : "
        f"{summary['fraction_displacement_ge_3A']:.2f}%"
    )

    print("\nGLOBAL RADIAL GEOMETRY")

    print(
        f"Native mean radius             : "
        f"{summary['native_mean_radius_A']:.3f} Å"
    )

    print(
        f"Disassembly mean radius        : "
        f"{summary['target_mean_radius_A']:.3f} Å"
    )

    print(
        f"Mean radius change             : "
        f"{summary['mean_radius_change_A']:+.3f} Å"
    )

    print(
        f"Median radius change           : "
        f"{summary['median_radius_change_A']:+.3f} Å"
    )

    print(
        f"Maximum radial expansion       : "
        f"{summary['max_radial_expansion_A']:+.3f} Å"
    )

    print(
        f"Residues moving outward        : "
        f"{summary['fraction_outward_radial']:.2f}%"
    )

    print(
        f"Residues moving inward         : "
        f"{summary['fraction_inward_radial']:.2f}%"
    )

    print("\nTANGENTIAL / ANGULAR DEFORMATION")

    print(
        f"Mean tangential displacement   : "
        f"{summary['mean_tangential_component_A']:.3f} Å"
    )

    print(
        f"Maximum tangential displacement: "
        f"{summary['max_tangential_component_A']:.3f} Å"
    )

    print(
        f"Mean angular change             : "
        f"{summary['mean_angular_change_deg']:.3f}°"
    )

    print(
        f"Median angular change           : "
        f"{summary['median_angular_change_deg']:.3f}°"
    )

    print(
        f"Maximum angular change          : "
        f"{summary['max_angular_change_deg']:.3f}°"
    )


# =============================================================
# PRINT CHAIN RESULTS
# =============================================================

def print_chain_results(chain_df):

    print("\n")
    print("=" * 100)
    print("CHAIN-LEVEL DEFORMATION")
    print("=" * 100)

    print(
        f"{'Chain':<10}"
        f"{'N Cα':>10}"
        f"{'Mean Δ':>14}"
        f"{'Max Δ':>14}"
        f"{'Mean radial':>16}"
        f"{'Mean tang.':>16}"
        f"{'Mean angular':>16}"
    )

    print("-" * 100)

    for _, row in chain_df.iterrows():

        print(
            f"{row['chain']:<10}"
            f"{int(row['n_CA']):>10}"
            f"{row['mean_displacement_A']:>14.3f}"
            f"{row['max_displacement_A']:>14.3f}"
            f"{row['mean_radial_change_A']:>16.3f}"
            f"{row['mean_tangential_A']:>16.3f}"
            f"{row['mean_angular_change_deg']:>16.3f}"
        )


# =============================================================
# PRINT INTERFACE RESULTS
# =============================================================

def print_interface_results(interface_df):

    print("\n")
    print("=" * 120)
    print("GLOBAL CHAIN-CHAIN INTERFACE DEFORMATION")
    print("=" * 120)

    if interface_df.empty:

        print(
            "No interfaces were detected."
        )

        return

    ranked = interface_df.sort_values(
        [
            "contact_loss_percent",
            "centroid_separation_change_A"
        ],
        ascending=False
    )

    print(
        f"{'Pair':<10}"
        f"{'Native':>10}"
        f"{'Lost':>10}"
        f"{'Loss %':>12}"
        f"{'Sep Δ':>14}"
        f"{'Mean d Δ':>14}"
        f"{'Max d Δ':>14}"
    )

    print("-" * 120)

    for _, row in ranked.head(
        N_TOP_INTERFACES
    ).iterrows():

        pair = (
            f"{row['chain_A']}-"
            f"{row['chain_B']}"
        )

        print(
            f"{pair:<10}"
            f"{int(row['native_contacts']):>10}"
            f"{int(row['lost_contacts']):>10}"
            f"{row['contact_loss_percent']:>12.2f}"
            f"{row['centroid_separation_change_A']:>14.3f}"
            f"{row['mean_contact_distance_change_A']:>14.3f}"
            f"{row['maximum_contact_distance_change_A']:>14.3f}"
        )


# =============================================================
# INTERACTIVE HTML
# =============================================================

def make_html(
    df,
    native_center,
    target_center,
    filename
):
    """
    Interactive global deformation map.

    Native structure:
        solid backbone traces

    Disassembly:
        dashed backbone traces

    Each C-alpha is represented by a marker whose size is
    proportional to displacement.

    No arrows are used.
    """

    if not MAKE_HTML:
        return

    try:
        import plotly.graph_objects as go

    except ImportError:

        print(
            "\nPlotly is not installed; HTML skipped."
        )

        return

    fig = go.Figure()

    # ---------------------------------------------------------
    # Group by chain
    # ---------------------------------------------------------

    for chain_id, group in df.groupby(
        "chain",
        sort=True
    ):

        # Native
        fig.add_trace(
            go.Scatter3d(
                x=group["native_x"],
                y=group["native_y"],
                z=group["native_z"],
                mode="lines",
                name=f"Native {chain_id}",
                line=dict(
                    width=5
                ),
                hoverinfo="skip"
            )
        )

        # Target
        fig.add_trace(
            go.Scatter3d(
                x=group["target_x"],
                y=group["target_y"],
                z=group["target_z"],
                mode="lines",
                name=f"Disassembly {chain_id}",
                line=dict(
                    width=4,
                    dash="dash"
                ),
                hoverinfo="skip"
            )
        )

    # ---------------------------------------------------------
    # Deformation points
    # ---------------------------------------------------------

    marker_sizes = (
        3.0 +
        3.0 *
        df["displacement_A"].to_numpy()
    )

    hover_text = []

    for _, row in df.iterrows():

        hover_text.append(
            f"Chain {row['chain']} "
            f"{row['resname']}{row['resnum']}<br>"
            f"Displacement: "
            f"{row['displacement_A']:.2f} Å<br>"
            f"Radial change: "
            f"{row['radial_change_A']:+.2f} Å<br>"
            f"Tangential: "
            f"{row['tangential_component_A']:.2f} Å<br>"
            f"Angular change: "
            f"{row['angular_change_deg']:.2f}°"
        )

    fig.add_trace(
        go.Scatter3d(
            x=df["target_x"],
            y=df["target_y"],
            z=df["target_z"],
            mode="markers",
            name="Disassembly deformation",
            marker=dict(
                size=marker_sizes,
                opacity=0.75
            ),
            text=hover_text,
            hovertemplate="%{text}<extra></extra>"
        )
    )

    # ---------------------------------------------------------
    # Centers
    # ---------------------------------------------------------

    fig.add_trace(
        go.Scatter3d(
            x=[native_center[0]],
            y=[native_center[1]],
            z=[native_center[2]],
            mode="markers",
            name="Native center",
            marker=dict(
                size=8,
                symbol="circle"
            )
        )
    )

    fig.add_trace(
        go.Scatter3d(
            x=[target_center[0]],
            y=[target_center[1]],
            z=[target_center[2]],
            mode="markers",
            name="Disassembly center",
            marker=dict(
                size=8,
                symbol="diamond"
            )
        )
    )

    fig.update_layout(
        title=(
            "FHV Global Capsid Deformation<br>"
            "Native vs Early Disassembly"
        ),

        scene=dict(
            xaxis_title="X (Å)",
            yaxis_title="Y (Å)",
            zaxis_title="Z (Å)",
            aspectmode="data"
        ),

        template="simple_white",

        legend=dict(
            x=0.01,
            y=0.99
        )
    )

    fig.write_html(
        filename,
        include_plotlyjs="cdn"
    )


# =============================================================
# MAIN
# =============================================================

def main():

    print("\n")
    print("=" * 80)
    print("FHV GLOBAL CAPS ID DEFORMATION ANALYSIS")
    print("=" * 80)

    print(
        f"\nNative      : {NATIVE_PDB}"
    )

    print(
        f"Disassembly : {DISASSEMBLY_PDB}"
    )

    # ---------------------------------------------------------
    # Load
    # ---------------------------------------------------------

    native = load_structure(
        NATIVE_PDB
    )

    target = load_structure(
        DISASSEMBLY_PDB
    )

    native_chains = get_chains(
        native
    )

    target_chains = get_chains(
        target
    )

    print("\nNative chains:")
    print(
        "    " +
        ", ".join(
            sorted(native_chains.keys())
        )
    )

    print("\nDisassembly chains:")
    print(
        "    " +
        ", ".join(
            sorted(target_chains.keys())
        )
    )

    # ---------------------------------------------------------
    # Common chains
    # ---------------------------------------------------------

    common_chains = sorted(
        set(native_chains.keys()) &
        set(target_chains.keys())
    )

    if len(common_chains) == 0:

        raise RuntimeError(
            "No common chain IDs found."
        )

    print("\nCommon chains:")
    print(
        "    " +
        ", ".join(common_chains)
    )

    analysis_native_chains = {
        cid: native_chains[cid]
        for cid in common_chains
    }

    analysis_target_chains = {
        cid: target_chains[cid]
        for cid in common_chains
    }

    # ---------------------------------------------------------
    # Global C-alpha matching
    # ---------------------------------------------------------

    print("\nMatching C-alpha atoms...")

    matched = match_global_ca(
        analysis_native_chains,
        analysis_target_chains
    )

    if len(matched) < 10:

        raise RuntimeError(
            "Too few matched C-alpha atoms."
        )

    print(
        f"Matched Cα atoms: {len(matched)}"
    )

    # ---------------------------------------------------------
    # Global deformation
    # ---------------------------------------------------------

    (
        deformation_df,
        native_center,
        target_center
    ) = calculate_global_deformation(
        matched
    )

    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------

    summary = calculate_summary(
        deformation_df
    )

    print_global_results(
        summary
    )

    # ---------------------------------------------------------
    # Chain statistics
    # ---------------------------------------------------------

    chain_df = calculate_chain_statistics(
        deformation_df
    )

    print_chain_results(
        chain_df
    )

    # ---------------------------------------------------------
    # Chain centroid statistics
    # ---------------------------------------------------------

    centroid_df = calculate_chain_centroids(
        analysis_native_chains,
        analysis_target_chains
    )

    # ---------------------------------------------------------
    # Interface analysis
    # ---------------------------------------------------------

    print("\n")
    print(
        "Analyzing all detected chain-chain interfaces..."
    )

    interface_df = analyze_all_interfaces(
        analysis_native_chains,
        analysis_target_chains
    )

    print_interface_results(
        interface_df
    )

    # ---------------------------------------------------------
    # Save residue deformation table
    # ---------------------------------------------------------

    deformation_file = (
        OUTPUT_PREFIX +
        "_per_residue.csv"
    )

    deformation_df.to_csv(
        deformation_file,
        index=False
    )

    print(
        f"\nSaved:\n    {deformation_file}"
    )

    # ---------------------------------------------------------
    # Save top residues
    # ---------------------------------------------------------

    ranked_df = save_top_residues(
        deformation_df
    )

    print(
        f"Saved:\n    "
        f"{OUTPUT_PREFIX}_top_deformed_residues.csv"
    )

    # ---------------------------------------------------------
    # Save chain statistics
    # ---------------------------------------------------------

    chain_file = (
        OUTPUT_PREFIX +
        "_chain_statistics.csv"
    )

    chain_df.to_csv(
        chain_file,
        index=False
    )

    print(
        f"Saved:\n    {chain_file}"
    )

    # ---------------------------------------------------------
    # Save chain centroid table
    # ---------------------------------------------------------

    centroid_file = (
        OUTPUT_PREFIX +
        "_chain_centroids.csv"
    )

    centroid_df.to_csv(
        centroid_file,
        index=False
    )

    print(
        f"Saved:\n    {centroid_file}"
    )

    # ---------------------------------------------------------
    # Save interface table
    # ---------------------------------------------------------

    interface_file = (
        OUTPUT_PREFIX +
        "_interfaces.csv"
    )

    interface_df.to_csv(
        interface_file,
        index=False
    )

    print(
        f"Saved:\n    {interface_file}"
    )

    # ---------------------------------------------------------
    # Save global summary
    # ---------------------------------------------------------

    summary_file = save_summary(
        summary,
        native_center,
        target_center
    )

    print(
        f"Saved:\n    {summary_file}"
    )

    # ---------------------------------------------------------
    # HTML
    # ---------------------------------------------------------

    if MAKE_HTML:

        html_file = (
            OUTPUT_PREFIX +
            "_3D.html"
        )

        make_html(
            deformation_df,
            native_center,
            target_center,
            html_file
        )

        print(
            f"Saved:\n    {html_file}"
        )

    # ---------------------------------------------------------
    # Final interpretation
    # ---------------------------------------------------------

    print("\n")
    print("=" * 80)
    print("GLOBAL INTERPRETATION")
    print("=" * 80)

    print(
        "\nThe global analysis separates deformation into:"
    )

    print(
        "    • radial deformation"
    )

    print(
        "    • tangential deformation"
    )

    print(
        "    • angular reorientation"
    )

    print(
        "    • inter-subunit interface opening/contact loss"
    )

    print(
        "\nThese measurements can be used to determine whether "
        "the early-disassembly structure represents a globally "
        "expanded and geometrically distorted capsid."
    )

    print(
        "\nFor the manuscript, distinguish geometric deformation "
        "from direct mechanical torsional strain energy."
    )

    print("\n")
    print("=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
