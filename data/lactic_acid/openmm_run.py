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
smiles

# %%
s_enantiomer = Molecule.from_smiles(smiles.iloc[4][2])
r_enantiomer = Molecule.from_smiles(smiles.iloc[5][2])

# %%
from openmmforcefields.generators import (
    GAFFTemplateGenerator,
)

gaff = GAFFTemplateGenerator(molecules=s_enantiomer)


# %%
# Input Files

pdb = PDBFile('s_lactic_acid_water.pdb')
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

steps = 5000000
equilibrationSteps = 100000
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'DeviceIndex':'2', 'Precision': 'mixed'}
dcdReporter = DCDReporter('s_trajectory.dcd', 500)
#pdbReporter = PDBReporter('s_trajectory.pdb', 500)
dataReporter = StateDataReporter('s_log.txt', 1000, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
checkpointReporter = CheckpointReporter('s_checkpoint.chk', 10000)

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
# Minimize and Equilibrate

print('Performing energy minimization...')
simulation.minimizeEnergy()
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

# Simulate

print('Simulating...')
simulation.reporters.append(dcdReporter)
#simulation.reporters.append(pdbReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0
simulation.step(steps)

# %%
