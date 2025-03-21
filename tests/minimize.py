# %%
from openmm import *
from openmm.app import *
from openmm import unit

# %%
# Input Files

pdb = PDBFile('rotated.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/opc.xml')

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
platformProperties = {'DeviceIndex':'1', 'Precision': 'mixed'}
dcdReporter = DCDReporter('trajectory.dcd', 500)
#pdbReporter = PDBReporter('r_trajectory.pdb', 500)
dataReporter = StateDataReporter('log.txt', 1000, totalSteps=steps,
    step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
checkpointReporter = CheckpointReporter('checkpoint.chk', 10000)

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
# set oxygen mass to zero
for atom in topology.atoms():
    if atom.name == 'O':
        mass = 0.0 * unit.dalton
        system.setParticleMass(atom.index, mass)


# %%
# Minimize and Equilibrate

print('Performing energy minimization...')
simulation.minimizeEnergy()
PDBFile.writeFile(topology, simulation.context.getState(getPositions=True).getPositions(), 'minimized.pdb')


# %%
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
#simulation.step(equilibrationSteps)

# Simulate

print('Simulating...')
simulation.reporters.append(dcdReporter)
#simulation.reporters.append(pdbReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0
#simulation.step(steps)

# %%
