import msprime
import pyslim
import numpy as np
import matplotlib
import glob
matplotlib.use('Agg')
from matplotlib import pyplot as plt

outdir = "."
treefiles = glob.glob(outdir + "pop_*.trees")

for treefile in treefiles:
    outfile = treefile + ".locs.png"
    ts = pyslim.load(treefile, slim_format=True)
    locations = np.array([i.location for i in ts.individuals()])
    fig = plt.figure(figsize=(max(locations[:,0]),max(locations[:,1])))
    plt.scatter(locations[:,0], locations[:,1], marker='.')
    plt.savefig(outfile, dpi=288)

