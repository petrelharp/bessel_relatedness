#!/usr/bin/env python3
description = '''Compute mean pairwise divergences across a discritization of the landscape.'''

import msprime
import pyslim
import numpy as np
import os, glob

### command-line parsing: execute this with '-h' for help
import argparse

parser = argparse.ArgumentParser(description=description)
parser.add_argument("--basedir", "-o", type=str, dest="basedir", required=True,
                    help="name of directory to save output files to.")
parser.add_argument("--logfile", "-g", type=str, dest="logfile", 
                    help="Name of log file [default: outdir/divergences.log]")
parser.add_argument("--num_samples", "-s", type=int, dest="num_samples", 
                    help="Number of samples.")
parser.add_argument("--mutation_rate", "-u", type=float, dest="mutation_rate", 
                    help="Mutation rate.")

args = parser.parse_args()

if not os.path.isdir(args.basedir):
    os.mkdir(args.basedir)

if args.logfile is None:
    args.logfile = os.path.join(args.basedir, "divergences.log")

logfile = open(args.logfile, "w")

print(args, file=logfile)

outdir = args.basedir
num_samples = args.num_samples
mutation_rate = args.mutation_rate

for treefile in glob.glob(os.path.join(outdir, "*.trees")):
    logfile.write("Reading {}\n".format(treefile))
    logfile.flush()
    decap = pyslim.load(treefile)

    # individuals, in random order
    sample_inds = [i for i in decap.individuals() if i.flags & pyslim.INDIVIDUAL_ALIVE]
    np.random.shuffle(sample_inds)
    # genomes corresponding to the first num_samples individuals
    samples = []
    for k in range(num_samples):
         samples.extend(sample_inds[k].nodes)

    # record info about the samples
    samplefile = open(treefile + ".samples.tsv", 'w')
    samplefile.write("\t".join(['id', 'x', 'y', 'individual']) + "\n")
    for a in samples:
        ind = decap.individual(decap.node(a).individual)
        samplefile.write("{}\t{}\t{}\t{}\n".format(a, ind.location[0],
    		 ind.location[1], decap.node(a).individual))
    samplefile.close()

    # simplify first for speed
    recap = decap.recapitate(recombination_rate=1e-9, Ne=1e3)
    sub_ts = pyslim.SlimTreeSequence(recap.simplify(samples))
    mut_ts = msprime.mutate(sub_ts, rate=mutation_rate)
    logfile.write("Recapitated; added {} mutations to {} trees. Now computing statistics.\n".format(mut_ts.num_mutations, mut_ts.num_trees))
    logfile.flush()
    
    # calculate divergences:
    # writes out the upper triangle, row-by-row
    windows = np.linspace(0.0, mut_ts.sequence_length, 2)
    bs = msprime.SiteStatCalculator(mut_ts)
    divs = np.array(bs.divergence([[x] for x in range(len(samples))], 
    			      windows=windows))
    outfile = treefile + ".divs.tsv"
    np.savetxt(outfile, divs, delimiter='\t')

logfile.flush()
