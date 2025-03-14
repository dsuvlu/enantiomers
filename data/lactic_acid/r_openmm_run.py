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
s_enantiomer = Molecule.from_smiles(smiles.iloc[4][2])
r_enantiomer = Molecule.from_smiles(smiles.iloc[5][2])

# %%
r_enantiomer.visualize()
s_enantiomer.visualize()

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

pdb = PDBFile('r_lactic_acid_water.pdb')
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
dcdReporter = DCDReporter('r_trajectory_1000ns.dcd', 500)
#pdbReporter = PDBReporter('s_trajectory.pdb', 500)
dataReporter = StateDataReporter('r_log.txt', 1000, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
checkpointReporter = CheckpointReporter('r_checkpoint.chk', 10000)

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
