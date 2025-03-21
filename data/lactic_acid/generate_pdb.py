# %%
import sys
from openff.toolkit import ForceField, Molecule, unit
from openmm.app import ForceField, PDBFile

# %%
sys.path.insert(0, '/home/dsuvlu/git/openmm_chiral_water/openmm_chiral_water')
import setup

# %%
smiles = setup.example_molecule_smiles()
smiles

# %%
s_enantiomer = Molecule.from_smiles(smiles.iloc[4][2])
s_enantiomer.generate_conformers(n_conformers=1)
s_top = Molecule.to_topology(s_enantiomer)
s_top.to_file('s_lactic_acid.pdb', file_format='pdb')
s_enantiomer.visualize()


# %%
r_enantiomer = Molecule.from_smiles(smiles.iloc[5][2])
r_enantiomer.generate_conformers(n_conformers=1)
r_top = Molecule.to_topology(r_enantiomer)
r_top.to_file('r_lactic_acid.pdb', file_format='pdb')
r_enantiomer.visualize()

# %%
