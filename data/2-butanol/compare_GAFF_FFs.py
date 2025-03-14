# %%
from openmm import *
from openmm.app import *
from openmm import unit
import sys
from openff.toolkit import Molecule

sys.path.insert(0, '/home/dsuvlu/git/enantiomers/src/python/')
import main

# %%
smiles = main.example_molecule_smiles()

# %%
s_enantiomer = Molecule.from_smiles(smiles.iloc[9][2])
r_enantiomer = Molecule.from_smiles(smiles.iloc[8][2])

# %%
s_enantiomer.assign_partial_charges('am1bcc')
r_enantiomer.assign_partial_charges('am1bcc')
r_enantiomer.partial_charges = s_enantiomer.partial_charges

# %%
from openmmforcefields.generators import (
    GAFFTemplateGenerator,
)

gaff = GAFFTemplateGenerator(molecules=r_enantiomer)


# %%
# Input Files
#os.chdir('/home/dsuvlu/git/openmm_chiral_water/openmm_chiral_water/data/2-butanol')
pdb = PDBFile('r_2-butanol.pdb')
forcefield = ForceField('amber14/opc.xml')
forcefield.registerTemplateGenerator(gaff.generator)

# %%
# System Configuration

nonbondedMethod = PME
nonbondedCutoff = 1.0*unit.nanometers
ewaldErrorTolerance = 0.0005
constraints = HBonds
rigidWater = True
constraintTolerance = 0.000001

# Integration Options

dt = 0.002*unit.picoseconds
temperature = 300*unit.kelvin
friction = 1.0/unit.picosecond
pressure = 1.0*unit.atmospheres
barostatInterval = 25

# Simulation Options

steps = 500000000
equilibrationSteps = 100000
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'DeviceIndex':'2', 'Precision': 'mixed'}
dcdReporter = DCDReporter('tmp.dcd', 500)
#pdbReporter = PDBReporter('s_trajectory.pdb', 500)
dataReporter = StateDataReporter('tmp.txt', 1000, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
checkpointReporter = CheckpointReporter('tmp.chk', 10000)

# %%
# Prepare the Simulation

print('Building system...')
topology = pdb.topology
positions = pdb.positions
system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))
integrator = LangevinMiddleIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)
simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)

# %%
with open("r_enantiomer-s-charges.xml", mode="w") as file:
    file.write(XmlSerializer.serialize(system))

# %%
