# %%
from openmm import *
from openmm.app import *
from openmm import unit

# %%
# Input Files

pdb = PDBFile('s_alanine.pdb')
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

steps = 500000000
equilibrationSteps = 100000
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'DeviceIndex':'2', 'Precision': 'mixed'}
dcdReporter = DCDReporter('r_trajectory_1000ns.dcd', 500)
#pdbReporter = PDBReporter('r_trajectory.pdb', 500)
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
with open("s_system.xml", mode="w") as file:
    file.write(XmlSerializer.serialize(system))

# %%
