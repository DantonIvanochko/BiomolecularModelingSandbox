#!/programs/x86_64-linux/system/sbgrid_bin/python.openmm

### simple simulation with openmm from sbgrid 

from openmm import app, unit, Platform, MonteCarloBarostat
import openmm as mm

from sys import stdout
from os.path import basename
from argparse import ArgumentParser



# testInput = '../test001.pdb'


def main(args): 

        inputFileBasename = basename(args.inputFile)
        
        print(f"Working on {inputFileBasename}")

        # load in input PDB file 
        pdb = app.PDBFile(args.inputFile)

        # force field XML files - explicit solvent
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml') # amber14-all.xml == amber14/protein.ff14SB.xml

        # use app.Modeller to add hydrogens and solvent
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield, pH=7.4)
        modeller.addSolvent(forcefield, padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)
        modeller.addExtraParticles(forcefield)
        app.PDBFile.writeFile(modeller.topology, modeller.positions, open(f'{inputFileBasename}-coordinates.pdb', 'w')) 

        # prepare system and integrator
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME,
                nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
                ewaldErrorTolerance=0.0005) # possibly use contrains=app.None
        integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1.0/unit.picoseconds,
                2.0*unit.femtoseconds)
    
        integrator.setConstraintTolerance(0.00001)

        # prepare simulation
        simulation = app.Simulation(modeller.topology, system, integrator, Platform.getPlatformByName('CUDA'))
        simulation.context.setPositions(modeller.positions)

        # for GPU calculations
        platformProperties = {'Precision','mixed'}

        ### equilibration

        # minimize
        simulation.minimizeEnergy() # default tolerance=10*kilojoule/mole

        # NVT equillibration for 1000 steps 
        simulation.context.setVelocitiesToTemperature(310*unit.kelvin)
        simulation.step(1000)

        # setup NPT equillibration for 1000 steps
        system.addForce(MonteCarloBarostat(1*unit.bar, 310*unit.kelvin))
        simulation.context.reinitialize(preserveState=True)

        # NPT simulated annealing for 200 ps
        for i in range(100):
                integrator.setTemperature((3*(1+i)+10)*unit.kelvin)
                simulation.step(100)

        # last NPT equillibration for 1000 steps
        simulation.step(1000)

        ### Production
 
        #STEPS
        #MD_steps = 5000000 # 5000000 x 2 fs = 10 ns
        #MD_steps = 25000000 # 50 ns
        MD_steps = 250000 # .5 ns

        # append reporters, take snapshot every 10 ps
        simulation.reporters.append(app.DCDReporter(f'{inputFileBasename}-trajectory.dcd', 5000))
        simulation.reporters.append(app.StateDataReporter(stdout, 5000, step=True,
                potentialEnergy=True, temperature=True, progress=True, remainingTime=True,
                speed=True, totalSteps=MD_steps, separator='\t'))

        # NPT production MD
        system.addForce(MonteCarloBarostat(1*unit.bar, 310*unit.kelvin))
        simulation.context.reinitialize(preserveState=True)

        # run production simulation
        simulation.step(MD_steps) 




if __name__ == "__main__":
        
        #Parse arguments
        parser = ArgumentParser()
        parser.add_argument("--inputFile", type=str, required=True, help="Input PDB file")
        args = parser.parse_args()

        main(args)






