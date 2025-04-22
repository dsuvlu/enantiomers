# %%
from openmm import *
from openmm.app import *
from openmm import unit
import sys
from openff.toolkit import Molecule

sys.path.insert(0, '/home/dsuvlu/git/openmm_chiral_water/openmm_chiral_water')
import setup

# %%
smiles = setup.example_molecule_smiles()
smiles

# %%
s_enantiomer = Molecule.from_smiles(smiles.iloc[0][2])
s_enantiomer.generate_conformers(n_conformers=1)
s_top = Molecule.to_topology(s_enantiomer)
s_enantiomer.visualize()

# %%
r_enantiomer = Molecule.from_smiles(smiles.iloc[1][2])
r_enantiomer.generate_conformers(n_conformers=1)
r_top = Molecule.to_topology(r_enantiomer)
r_enantiomer.visualize()


# %%
