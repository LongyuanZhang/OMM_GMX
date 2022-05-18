import numpy as np
import csv
import networkx as nx
import random
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from sys import argv
from math import sqrt
import parmed
from parmed import gromacs

import time

# Gromacs force field library

gromacs.GROMACS_TOPDIR = "/home/lzhang657/anaconda3/envs/CLIPS2/share/gromacs/top"
output =  open('ti.txt', 'w')

parm = parmed.load_file("ghost.top",xyz="solv_ions.gro")
positions = parm.positions
#print(type(parm))
system = parm.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometers, constraints=HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)

forces = system.getForces()
for force in forces:
    if type(force) is NonbondedForce:
        force.setUseDispersionCorrection(True)
        nonbondedforce = force
        break

# comb_rule = 3 leads to a ParmEd defined customNB == nonbonded LJ
for force in forces:
    if type(force) is CustomNonbondedForce:
        #force.setUseLongRangeCorrection(True)
        CustomLJbyParmed = force
        break

index = 1020
alchemical_particles = {index,index+1}
chemical_particles = set(range(system.getNumParticles())) - alchemical_particles

customnonbond = CustomNonbondedForce('lambda*epsilon*x*(x-1.0);'+
                                     'x = sigma^6/(0.5*(1.0-lambda)*sigma^6 + r^6);'+
                                     'sigma = sigma1*sigma2; epsilon = epsilon1*epsilon2; lambda = sqrt(lambda1*lambda2)')
# set LJ Scaling to 0.0

customnonbond.addPerParticleParameter('epsilon')
customnonbond.addPerParticleParameter('sigma')
customnonbond.addPerParticleParameter('lambda')
#for i in range(system.getNumParticles()):
#    [charge, sigma, epsilon] = nonbondedforce.getParticleParameters(i)
#    customnonbond.addParticle([sigma, epsilon])
#    if i in alchemical_particles:
#        nonbondedforce.setParticleParameters(i, 0.0, sigma, 0.0)

for i in range(system.getNumParticles()):
    [epsilon, sigma] = CustomLJbyParmed.getParticleParameters(i)
    customnonbond.addParticle([epsilon, sigma, 0.0])
    if i in alchemical_particles:
        CustomLJbyParmed.setParticleParameters(i, (0.0, sigma))
        nonbondedforce.setParticleParameters(i, 0.0, sigma, 0.0)

customnonbond.addInteractionGroup(alchemical_particles, chemical_particles)
customnonbond.addInteractionGroup({index},{index+1})
#customnonbond.setUseLongRangeCorrection(True)
customnonbond.setCutoffDistance(nonbondedforce.getCutoffDistance())
customnonbond.setNonbondedMethod(2)
for i in range(nonbondedforce.getNumExceptions()):
    (p1, p2, q, sig, eps) = nonbondedforce.getExceptionParameters(i)
    # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
    # FORCE
    #print("p1 = %.2f, p2 = %.2f"%(p1, p2))
    customnonbond.addExclusion(p1, p2)

system.addForce(customnonbond)

integrator = LangevinIntegrator(298.15*kelvin, 1.0/picoseconds,
    2.0*femtoseconds)

# No Barostat bec. no solvent
#system.addForce(MonteCarloBarostat(1*atmospheres, 298.15*kelvin, 25))

platform = Platform.getPlatformByName('CPU')

simulation = Simulation(parm.topology, system, integrator, platform)
simulation.context.setPositions(positions)

def setParameters(LJScaling, chargeScaling):
    nonbondedforce.setParticleParameters(index,1.0*chargeScaling,0.2584,0.0)
    nonbondedforce.setParticleParameters(index+1,-1.0*chargeScaling,0.4036,0.0)
    #simulation.context.setParameter('lambda', LJScaling)
    nonbondedforce.updateParametersInContext(simulation.context)
    for i in range(system.getNumParticles()):
        [epsilon, sigma, lambda0] = customnonbond.getParticleParameters(i)
        customnonbond.setParticleParameters(i, [epsilon, sigma, LJScaling])
    customnonbond.updateParametersInContext(simulation.context)

setParameters(0.0, 0.0)
simulation.context.setVelocitiesToTemperature(298.15*kelvin)
simulation.reporters.append(app.DCDReporter('traj.dcd', 5000))

nequlibrium = 20000  #run short (10–100 ps) simulations to equlibrate at each lambda state
#nequlibrium = 10000  #run short (10–100 ps) simulations to equlibrate at each lambda state
#print(parmed.openmm.energy_decomposition_system(parm, system,nrg=kilojoules_per_mole))

simulation.step(nequlibrium)
nsteps =1000
niterations = [125, 250]

x1,y1 = np.polynomial.legendre.leggauss(10)
t1 = 0.5*(x1+1)*(1.0-0.0)+0.0
w1 = y1*0.5*(1.0-0.0)

x2,y2 = np.polynomial.legendre.leggauss(10)
t2 = 0.5*(x2+1)*(1.0-0.0)+0.0
w2 = y2*0.5*(1.0-0.0)

LJScalings = np.concatenate([t1,np.ones(10)])
chargeScalings = np.concatenate([np.zeros(10),t2])

nstates = len(LJScalings) 
derivatives = [0]*nstates
    
chargeScalings = np.concatenate([np.zeros(10),t2]) 

nstates = len(LJScalings)
derivatives = [0]*nstates

kT = AVOGADRO_CONSTANT_NA*BOLTZMANN_CONSTANT_kB*integrator.getTemperature()
for k in range(nstates):
    print(k)
    setParameters(LJScalings[k], chargeScalings[k])
    simulation.step(nequlibrium)
    
    for iteration in range(niterations[k//10]):
        setParameters(LJScalings[k], chargeScalings[k])
        output.write('state %5d iteration %5d / %5d\n' % (k, iteration, niterations[k//10]))
        simulation.step(nsteps)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()/kT
        if (k <= 10):
            setParameters(LJScalings[k]-0.0005, chargeScalings[k])
            a = simulation.context.getState(getEnergy=True).getPotentialEnergy()/kT
            setParameters(LJScalings[k]+0.0005, chargeScalings[k])
            b = simulation.context.getState(getEnergy=True).getPotentialEnergy()/kT
            output.write("energy derivative%f\n" %((b-a)/0.001))
            derivatives[k] += (b-a)/0.001
        else:
            setParameters(LJScalings[k], chargeScalings[k]-0.0005)
            a = simulation.context.getState(getEnergy=True).getPotentialEnergy()/kT
            setParameters(LJScalings[k], chargeScalings[k]+0.0005)
            b = simulation.context.getState(getEnergy=True).getPotentialEnergy()/kT
            output.write("energy derivative%f\n" %((b-a)/0.001))
            derivatives[k] += (b-a)/0.001
lj = 0
for i in range(10):
    lj += derivatives[i]*w1[i]
lj = lj/niterations[0]
elec = 0
for i in range(10):
    elec += derivatives[i+10]*w2[i]
elec = elec/niterations[1]
output.write("LJ %f elec %f total %f\n" %(lj, elec, lj+elec))

#output.write("%f\n" %lj+elec)
            
with open('nacl_out.pdb', 'w') as config:
    boxVec = simulation.context.getState(getPositions=True).getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxVec)
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), config)

#output.write("LJ + elec + free_energy %f\n" %lj+elec)
#output.write("LJ + elec %f\n" %lj+elec)
output.write("time: %f s" %(time.time()-start_t))
output.close()
