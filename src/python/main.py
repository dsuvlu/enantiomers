"""Provide the primary functions."""
import pandas as pd
from scipy.constants import Avogadro
from openff.toolkit.topology import Molecule, Topology

def example_molecule_smiles():
    # Define the data as a list of dictionaries
    data = [
        {'Molecule': 'Alanine', 'Enantiomer': 'L-Alanine (S)', 'SMILES': '[H:2][C:1]([H:3])([H:4])[C:5](=[O:6])[N:7]([H:8])[C@:9]([H:10])([C:11]([H:12])([H:13])[H:14])[C:15](=[O:16])[N:17]([H:18])[C:19]([H:20])([H:21])[H:22]'},
        {'Molecule': 'Alanine', 'Enantiomer': 'D-Alanine (R)', 'SMILES': '[H:2][C:1]([H:3])([H:4])[C:5](=[O:6])[N:7]([H:8])[C@@:9]([H:10])([C:11]([H:12])([H:13])[H:14])[C:15](=[O:16])[N:17]([H:18])[C:19]([H:20])([H:21])[H:22]'},
        {'Molecule': 'Alanine', 'Enantiomer': 'L-Alanine (S)', 'SMILES': 'CC(=O)N[C@H](C)C(=O)N(C)'},
        {'Molecule': 'Alanine', 'Enantiomer': 'D-Alanine (R)', 'SMILES': 'CC(=O)N[C@@H](C)C(=O)N(C)'},
        {'Molecule': 'Lactic Acid', 'Enantiomer': '(S)-Lactic Acid', 'SMILES': 'C[C@H](O)C(=O)O'},
        {'Molecule': 'Lactic Acid', 'Enantiomer': '(R)-Lactic Acid', 'SMILES': 'C[C@@H](O)C(=O)O'},
        {'Molecule': 'Glyceraldehyde', 'Enantiomer': 'D-Glyceraldehyde (R)', 'SMILES': 'O=CC([C@H](O))CO'},
        {'Molecule': 'Glyceraldehyde', 'Enantiomer': 'L-Glyceraldehyde (S)', 'SMILES': 'O=CC([C@@H](O))CO'},
        {'Molecule': '2-Butanol', 'Enantiomer': '(R)-2-Butanol', 'SMILES': 'CC[C@H](O)C'},
        {'Molecule': '2-Butanol', 'Enantiomer': '(S)-2-Butanol', 'SMILES': 'CC[C@@H](O)C'},
        {'Molecule': '2-Chlorobutane', 'Enantiomer': '(R)-2-Chlorobutane', 'SMILES': 'CC[C@H](Cl)C'},
        {'Molecule': '2-Chlorobutane', 'Enantiomer': '(S)-2-Chlorobutane', 'SMILES': 'CC[C@@H](Cl)C'},
        {'Molecule': 'Mandelic Acid', 'Enantiomer': '(R)-Mandelic Acid', 'SMILES': 'O=C(O)[C@H](O)c1ccccc1'},
        {'Molecule': 'Mandelic Acid', 'Enantiomer': '(S)-Mandelic Acid', 'SMILES': 'O=C(O)[C@@H](O)c1ccccc1'},
        {'Molecule': 'Epichlorohydrin', 'Enantiomer': '(R)-Epichlorohydrin', 'SMILES': 'Cl[C@H]1CO1'},
        {'Molecule': 'Epichlorohydrin', 'Enantiomer': '(S)-Epichlorohydrin', 'SMILES': 'Cl[C@@H]1CO1'},
        {'Molecule': 'Propylene Oxide', 'Enantiomer': '(R)-Propylene Oxide', 'SMILES': 'C[C@H]1CO1'},
        {'Molecule': 'Propylene Oxide', 'Enantiomer': '(S)-Propylene Oxide', 'SMILES': 'C[C@@H]1CO1'},
        {'Molecule': 'Carvone', 'Enantiomer': '(R)-Carvone', 'SMILES': 'CC(=C)[C@H]1CC=C(C)C(=O)C1'},
        {'Molecule': 'Carvone', 'Enantiomer': '(S)-Carvone', 'SMILES': 'CC(=C)[C@@H]1CC=C(C)C(=O)C1'},
        {'Molecule': 'Limonene', 'Enantiomer': '(R)-Limonene', 'SMILES': 'CC(=C)[C@H]1CC=C(C)C=C1'},
        {'Molecule': 'Limonene', 'Enantiomer': '(S)-Limonene', 'SMILES': 'CC(=C)[C@@H]1CC=C(C)C=C1'},
    ]

    # Create the DataFrame
    df = pd.DataFrame(data)

    # Return the DataFrame
    return df

def alanine_pdb():
    smiles = example_molecule_smiles()

    s_enantiomer = Molecule.from_mapped_smiles(smiles.iloc[0][2])
    s_enantiomer.generate_conformers(n_conformers=1)
    s_enantiomer.visualize()

    top_s = Topology.from_molecules(s_enantiomer)
    top_s.to_file('alanine_s.pdb', file_format='pdb')

    r_enantiomer = Molecule.from_smiles(smiles.iloc[1][2])
    r_enantiomer.generate_conformers(n_conformers=1)
    r_enantiomer.visualize()

    top_r = Topology.from_molecules(r_enantiomer)
    top_r.to_file('alanine_r.pdb', file_format='pdb')


def water_box(L):
    """Calculate number of water molecules needed for a box with
    side length L."""

    # Create a box of water molecules
    lx = L
    ly = L
    lz = L
    volume = lx * ly * lz * 1e-24 # cm^3
    
    density = 0.993 # g/cm^3
    conc_wat = 55.5 # mol/L
    nwat = density*volume*Avogadro/18.0
    nwat = round(nwat)
    nwat = int(nwat)
    
    return nwat

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(example_molecule_smiles())
    print(water_box(30))
