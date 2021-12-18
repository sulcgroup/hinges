import numpy as np
from sys import path, argv, stderr
from os import getenv
import matplotlib.pyplot as plt
from multiprocessing import Pool

#set the path to include oxDNA analysis tools
path.append('/home/epopplet/software/')

from oxdna_analysis_tools.UTILS.readers import LorenzoReader2, get_input_parameter
from oxdna_analysis_tools.output_bonds import output_bonds

#get id file
with open('/scratch/epopplet/hinge/structures/ss_ids.txt', 'r') as f:
    ids = [int(i) for i in f.read().split(' ')]

id_dict = {n : i for i, n in enumerate(ids)}

#get files
prefix = argv[1]
inputfile = './'+prefix+'1/input'
topology = './'+prefix+'1/'+get_input_parameter(inputfile, "topology")

def get_hinge_bonds(file):
    interaction_matrix = np.zeros((len(ids), len(ids)))
    r = LorenzoReader2(file, topology)
    s = r._get_system()
    bond_count = []
    count = 0

    while s != False:
        s.map_nucleotides_to_strands()
        out = output_bonds(inputfile, s)
        s.read_H_bonds_output_bonds(out)
        bonds_at_t = 0

        for ii, i in enumerate(ids):
            for j in s._nucleotides[i].interactions:
                #print(i, j)
                bonds_at_t += 1
                try:
                    interaction_matrix[ii][id_dict[j]] += 1
                except:
                    print("Unexpected interaction:", i, j)

        bond_count.append(bonds_at_t / 2)
        s = r._get_system()
        count += 1


    return (interaction_matrix / count, bond_count)

# NEEDS TO BE FIXED SO IT CAN RUN ON PULL/RELEASE
results = Pool().map(get_hinge_bonds, ["{n}{i}/aligned.dat".format(n=prefix, i=i) for i in range(1, int(getenv('SLURM_CPUS_PER_TASK'))+1)])
counts = [r[1] for r in results]
bonds = [r[0] for r in results]

fig, ax = plt.subplots()
for i, d in enumerate(counts):
    ax.plot(np.arange(len(d)), d, label='run{}'.format(i))
plt.xlabel('Configuration')
plt.ylabel('Bonds in flexure')
plt.legend()
plt.savefig('{}_bonds.png'.format(prefix))

with open("{}_bonds.txt".format(prefix), 'w+') as f:
    for d in counts:
        [f.write('{} '.format(v)) for v in d]
        f.write('\n')
print('wrote data to {}_bonds.txt'.format(prefix))

for i, r in enumerate(bonds):
    #print(results, flush=True)
    fig, ax = plt.subplots()
    a = ax.imshow(r, cmap='viridis', origin='lower')
    ax.set(ylabel="Nucleotide", xlabel="Nucleotide")
    b = fig.colorbar(a, ax=ax)
    b.set_label("Interaction_freq", rotation=270)

    plt.axvspan(0, 52, ymin=-0.03, ymax=0, color='#1f77b4', clip_on=False)
    plt.axvspan(53, 105, ymin=-0.03, ymax=0, color='#ff7f0e', clip_on=False)
    plt.axvspan(106, 159, ymin=-0.03, ymax=0, color='#2ca02c', clip_on=False)
    plt.axvspan(160, 212, ymin=-0.03, ymax=0, color='#d62728', clip_on=False)
    plt.axvspan(213, 265, ymin=-0.03, ymax=0, color='#9467bd', clip_on=False)
    plt.axvspan(266, 318, ymin=-0.03, ymax=0, color='#8c564b', clip_on=False)
    plt.axhspan(0, 52, xmin=-0.03, xmax=0, color='#1f77b4', clip_on=False)
    plt.axhspan(53, 105, xmin=-0.03, xmax=0, color='#ff7f0e', clip_on=False)
    plt.axhspan(106, 159, xmin=-0.03, xmax=0, color='#2ca02c', clip_on=False)
    plt.axhspan(160, 212, xmin=-0.03, xmax=0, color='#d62728', clip_on=False)
    plt.axhspan(213, 265, xmin=-0.03, xmax=0, color='#9467bd', clip_on=False)
    plt.axhspan(266, 318, xmin=-0.03, xmax=0, color='#8c564b', clip_on=False)
    plt.savefig(prefix+str(i)+".png")
