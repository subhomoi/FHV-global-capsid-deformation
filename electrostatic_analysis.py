#!/usr/bin/env python3
"""PROPKA-based protonation analysis of acidic residues in FHV 4FTB."""
from __future__ import annotations
import argparse
import subprocess
import urllib.request
from pathlib import Path
import pandas as pd

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
ACIDIC_RESIDUES = {"ASP", "GLU"}


def download_pdb(pdb_id: str, output: Path) -> None:
    """Download a PDB coordinate file from RCSB."""
    url = RCSB_URL.format(pdb_id=pdb_id.upper())
    print(f"Downloading {pdb_id.upper()} from RCSB PDB...")
    urllib.request.urlretrieve(url, output)


def clean_pdb(input_path: Path, output_path: Path) -> None:
    """Remove HOH and EPE records, preserving all other PDB records."""
    with input_path.open() as src, output_path.open("w") as dst:
        for line in src:
            if line.startswith("HETATM") and line[17:20] in {"HOH", "EPE"}:
                continue
            dst.write(line)


def run_propka(pdb_path: Path) -> Path:
    """Run PROPKA and return its .pka output file."""
    subprocess.run(["propka3", str(pdb_path)], check=True)
    pka_path = pdb_path.with_suffix(".pka")
    if not pka_path.exists():
        raise FileNotFoundError(f"PROPKA output not found: {pka_path}")
    return pka_path


def protonated_fraction(pka: float, ph: float) -> float:
    """Calculate the protonated fraction of an acidic group."""
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def parse_propka(pka_path: Path, threshold: float = 40.0) -> pd.DataFrame:
    """Extract ASP/GLU residues with >= threshold percentage-point gain."""
    results = []
    reading = False

    with pka_path.open() as handle:
        for line in handle:
            if "SUMMARY OF THIS PREDICTION" in line:
                reading = True
                continue
            if not reading:
                continue
            if "Free energy of" in line:
                break
            if "----------------" in line or line.startswith("   Residue") or not line.strip():
                continue

            parts = line.split()
            if len(parts) < 4:
                continue
            res_name, res_num, chain = parts[0], parts[1], parts[2]
            try:
                pka = float(parts[3])
                res_num = int(res_num)
            except ValueError:
                continue
            if res_name not in ACIDIC_RESIDUES:
                continue

            pct_70 = 100.0 * protonated_fraction(pka, 7.0)
            pct_55 = 100.0 * protonated_fraction(pka, 5.5)
            net_gain = pct_55 - pct_70
            if net_gain >= threshold:
                results.append({
                    "Residue": res_name,
                    "ResNum": res_num,
                    "Chain": chain,
                    "Predicted_pKa": round(pka, 2),
                    "Protonated_%_pH7.0": round(pct_70, 1),
                    "Protonated_%_pH5.5": round(pct_55, 1),
                    "Net_Proton_Gain_%": round(net_gain, 1),
                })

    return pd.DataFrame(results)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb-id", default="4FTB")
    parser.add_argument("--threshold", type=float, default=40.0)
    parser.add_argument("--output", default="4FTB_major_neutralization_sites.csv")
    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    pdb_path = Path(f"{pdb_id}.pdb")
    clean_path = Path(f"{pdb_id}_clean.pdb")

    download_pdb(pdb_id, pdb_path)
    clean_pdb(pdb_path, clean_path)
    pka_path = run_propka(clean_path)
    df = parse_propka(pka_path, args.threshold)

    if df.empty:
        print(f"No residues met the >= {args.threshold:g}% net gain threshold.")
        return

    df = df.sort_values(["ResNum", "Chain"]).reset_index(drop=True)
    df["Residue"] = df["Residue"] + df["ResNum"].astype(str)
    df = df.drop(columns="ResNum")
    df.to_csv(args.output, index=False)
    print(f"Found {len(df)} major neutralization sites.")
    print(f"Data saved to {args.output}")


if __name__ == "__main__":
    main()
