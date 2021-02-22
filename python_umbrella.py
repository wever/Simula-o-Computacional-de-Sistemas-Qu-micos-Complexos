# Example illustrating the application of MBAR to compute a 1D PMF from an umbrella sampling simulation.
#
# The data represents an umbrella sampling simulation for the dist  of a valine sidechain in lysozyme L99A with benzene bound in the cavity.
#
# REFERENCE
#
# D. L. Mobley, A. P. Graves, J. D. Chodera, A. C. McReynolds, B. K. Shoichet and K. A. Dill, "Predicting absolute ligand binding free energies to a simple model site," Journal of Molecular Biology 371(4):1118-1134 (2007).
# http://dx.doi.org/10.1016/j.jmb.2007.06.002

from __future__ import print_function
import numpy # numerical array library
import pymbar # multistate Bennett acceptance ratio
from pymbar import timeseries # timeseries analysis

import matplotlib.pyplot as plt
import numpy as np
# Constants.
kB = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K

temperature = 298.15 # assume a single temperature -- can be overridden with data from center.dat
# Parameters
K = 13 # number of umbrellas
N_max = 50001 # maximum number of snapshots/simulation
T_k = numpy.ones(K,float)*temperature # inital temperatures are all equal
beta = 1.0 / (kB * temperature) # inverse temperature of simulations (in 1/(kJ/mol))
dist_min = 1 # min for PMF
dist_max = 2.2 # max for PMF
nbins = 36 # number of bins for 1D PMF

# Allocate storage for simulation data
N_k = numpy.zeros([K], dtype = int) # N_k[k] is the number of snapshots from umbrella simulation k
K_k = numpy.zeros([K]) # K_k[k] is the spring constant (in kJ/mol/nm**2) for umbrella simulation k
dist0_k = numpy.zeros([K]) # dist0_k[k] is the spring center location (in nm) for umbrella simulation k
dist_kn = numpy.zeros([K,N_max]) # dist_kn[k,n] is the torsion angle (in nm) for snapshot n from umbrella simulation k
u_kn = numpy.zeros([K,N_max]) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
g_k = numpy.zeros([K]);

for k in range(K):
    # Read enthalpy data.
    filename = 'data/%d.Enthalpy.xvg' % k
    print("Reading %s..." % filename)
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    # Parse data.
    n = 0
    for line in lines:
        if line[0] != '#' and line[0] != '@':
            tokens = line.split()
            u_kn[k,n] = float(tokens[1]) # torsion angle
    
            n+= 1



# Read in umbrella spring constants and centers.
infile = open('data/centers.dat', 'r')
lines = infile.readlines()
infile.close()
for k in range(K):
    # Parse line k.
    line = lines[k]
    tokens = line.split()
    dist0_k[k] = float(tokens[0]) # dist center locatiomn (in nm)
    K_k[k] = float(tokens[1])
    if len(tokens) > 2:
        T_k[k] = float(tokens[2])  # temperature the kth simulation was run at.

beta_k = 1.0/(kB*T_k)   # beta factor for the different temperatures

for k in range(K):
    # Read distance data.
    filename = 'data/%d.dist.xvg' % k
    print("Reading %s..." % filename)
    infile = open(filename, 'r')
    lines = infile.readlines()
    infile.close()
    # Parse data.
    n = 0
    for line in lines:
        if line[0] != '#' and line[0] != '@':
            tokens = line.split()
            dist = float(tokens[1]) # torsion angle
            # wrap dist_kn to be within [dist_min,dist_max)
            dist_kn[k,n] = dist
           
            n += 1
    N_k[k] = n
            
    binPlot, edges = np.histogram(dist_kn[k], bins=50)
    plt.plot(edges[:-1], binPlot, label=k)


    dist_ = dist_kn[k,0:N_k[k]]
    g_dist = timeseries.statisticalInefficiency(dist_)
    print("g_dist = %.1f " % (g_dist))
    g_k[k] = g_dist
    print("Correlation time for set %5d is %10.3f" % (k,g_k[k]))
    indices = timeseries.subsampleCorrelatedData(dist_, g=g_k[k])
    # Subsample data.
    N_k[k] = len(indices)
    u_kn[k,0:N_k[k]] = u_kn[k,indices]
    dist_kn[k,0:N_k[k]] = dist_kn[k,indices]

plt.xlabel(r"$nm$")
plt.ylabel(r"$Frequency$")
plt.legend()        
plt.show()
plt.close('all')

N_max = numpy.max(N_k) # shorten the array size
u_kln = numpy.zeros([K,K,N_max], numpy.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l

# Set zero of u_kn -- this is arbitrary.
u_kn -= u_kn.min()

# Construct torsion bins
print("Binning data...")
delta = (dist_max - dist_min) / float(nbins)
# compute bin centers
bin_center_i = numpy.zeros([nbins], numpy.float64)
for i in range(nbins):
    bin_center_i[i] = dist_min + delta/2 + delta * i
# Bin data
bin_kn = numpy.zeros([K,N_max], numpy.int32)
for k in range(K):

    for n in range(N_k[k]):
        # Compute bin assignment.
        bin_kn[k,n] = int((dist_kn[k,n] - dist_min) / delta)

# Evaluate reduced energies in all umbrellas
print("Evaluating reduced potential energies...")
for k in range(K):
    for n in range(N_k[k]):
        # Compute minimum-image torsion deviation from umbrella center l
        ddist = dist_kn[k,n] - dist0_k 
       
        # Compute energy of snapshot n from simulation k in umbrella potential l
        u_kln[k,:,n] = u_kn[k,n] + beta_k[k] * (K_k/2.0) * ddist**2

# Initialize MBAR.
print("Running MBAR...")
mbar = pymbar.MBAR(u_kln, N_k, verbose = True)

# Compute PMF in unbiased potential (in units of kT).
results = mbar.computePMF(u_kn, bin_kn, nbins, return_dict=True)
f_i = results['f_i']
df_i = results['df_i']

# Write out PMF
print("PMF (in units of kT)")
print("%8s %8s %8s" % ('bin', 'f', 'df'))
for i in range(nbins):
    print("%8.1f %8.3f %8.3f" % (bin_center_i[i], f_i[i], df_i[i]))

target_temperature = 300
target_beta = 1.0 / (8.314e-3 * target_temperature)

f_i /= target_beta
df_i /= target_beta

dih = np.array(bin_center_i)

plt.errorbar(dih, f_i, yerr=df_i)

plt.xlabel(r"$nm$")
plt.ylabel(r"$\Delta F (kJ/mol)$")
plt.show()
