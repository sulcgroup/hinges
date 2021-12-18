import numpy as np
from sys import path, argv, stderr
import matplotlib.pyplot as plt
from multiprocessing import Pool
from os import getenv

#set the path to include oxDNA analysis tools
path.append('/home/epopplet/software/oxdna_analysis_tools')

from UTILS.readers import ErikReader

if len(argv) != 4:
    print("Usage is {} directory_prefix arm1_index_file arm2_index_file".format(argv[0]))
prefix = argv[1]
arm1_file = argv[2]
arm2_file = argv[3]

#Get arm IDs
with open(arm1_file, 'r') as f:
    id1 = f.read().split()
    
with open(arm2_file, 'r') as f:
    id2 = f.read().split()
    
id1 = [int(i) for i in id1]
id2 = [int(i) for i in id2]

#the commented out ones are for the old design only found in hinge_correct and hinge_no
#arm1_ref_id = [2609, 6336] #for hinge
arm1_ref_id = [2034, 6657] #for all others
#arm2_ref_id = [2687, 2864]
arm2_ref_id = [4593, 4411]
def analyze_data(file):
    r = ErikReader(file)
    s = r.read()
    angles = []
    distances = []
    
    while s:
        #s.inbox()
        print(s.time, file=stderr, flush=True)
    
        arm1 = s.positions[id1]
        arm2 = s.positions[id2]
    
        arm1_ref = s.positions[arm1_ref_id[1]] - s.positions[arm1_ref_id[0]]
        arm2_ref = s.positions[arm2_ref_id[1]] - s.positions[arm2_ref_id[0]]

        mean1 = arm1.mean(axis=0)
        mean2 = arm2.mean(axis=0)

        #SVD on mean-centered data
        uu1, dd1, vv1 = np.linalg.svd(arm1 - mean1)
        uu2, dd2, vv2 = np.linalg.svd(arm2 - mean2)
    
        #flip vectors that point the wrong direction
        if np.dot(vv1[0], arm1_ref) < 0:
            vv1[0] = -vv1[0]
        if np.dot(vv2[0], arm2_ref) < 0:
            vv2[0] = -vv2[0]
    
        angles.append(np.arccos(np.dot(vv1[0], vv2[0]))*180/np.pi)
        distances.append(np.linalg.norm(mean1-mean2)*0.85)
    
        s = r.read()
        
    return(angles, distances)

#bring threads back together
results = Pool().map(analyze_data, ["{n}{i}/aligned.dat".format(n=prefix, i=i) for i in range(1, int(getenv('SLURM_CPUS_PER_TASK'))+1)])
#["{}1/aligned.dat".format(prefix), "{}2/aligned.dat".format(prefix), "{}3/aligned.dat".format(prefix)])#, "{}4/aligned.dat".format(prefix)])

angles = [a[0] for a in results]
distances = [d[1] for d in results]

#angles1 = results[0][0]
#distances1 = results[0][1]
#angles2 = results[1][0]
#distances2 = results[1][1]
#angles3 = results[2][0]
#distances3 = results[2][1]
#angles4 = results[3][0]
#distances4 = results[3][1]

#trajectory plots
for i, a in enumerate(angles):
    plt.plot(np.arange(len(a)), a, label='run{}'.format(i+1))
#plt.plot(np.arange(len(angles2)), angles2, label='run2')
#plt.plot(np.arange(len(angles3)), angles3, label='run3')
#plt.plot(np.arange(len(angles4)), angles4, label='run4')
plt.xlabel('Configuration')
plt.ylabel('Angle(deg)')
plt.legend()
plt.savefig('{}_angles.png'.format(prefix))

plt.clf()
for i, d in enumerate(distances):
    plt.plot(np.arange(len(d)), d, label='run{}'.format(i+1))
#plt.plot(np.arange(len(distances2)), distances2, label='run2')
#plt.plot(np.arange(len(distances3)), distances3, label='run3')
#plt.plot(np.arange(len(distances4)), distances4, label='run4')
plt.xlabel('Configuration')
plt.ylabel('Distance(nm)')
plt.legend()
plt.savefig('{}_distances.png'.format(prefix))


#histograms
lower_angles = min([min(a) for a in angles])
#lower_angles = min(min(angles1), min(angles2), min(angles3), min(angles4))
upper_angles = max([max(a) for a in angles])
#upper_angles = max(max(angles1), max(angles2), max(angles2), max(angles4))

lower_dists = min([min(d) for d in distances])
#lower_dists = min(min(distances1), min(distances2), min(distances3))#, min(angles4))
upper_dists = max([max(d) for d in distances])
#upper_dists = max(max(distances1), max(distances2), max(distances3))#, max(angles4))

plt.clf()
for a in angles:
    print(a, sep=' ')
#print(angles1, angles2, angles3)#, angles4)
bins = np.linspace(np.floor(lower_angles-(lower_angles*0.1)), np.ceil(upper_angles+(upper_angles*0.1)), 60)
for i, a in enumerate(angles):
    plt.hist(a, bins, weights=np.ones(len(a)) / len(a), alpha=0.5, histtype=u'stepfilled', edgecolor='k', label='run{}'.format(i))
#plt.hist(angles2, bins, weights=np.ones(len(angles2)) / len(angles2), alpha=0.5, histtype=u'stepfilled', edgecolor='k', label='run2')
#plt.hist(angles3, bins, weights=np.ones(len(angles3)) / len(angles3), alpha=0.5, histtype=u'stepfilled', edgecolor='k', label='run3')
#plt.hist(angles4, bins, weights=np.ones(len(angles4)) / len(angles4), alpha=0.5, histtype=u'stepfilled', edgecolor='k', label='run4')
plt.xlabel('Angle(deg)')
plt.ylabel('Normalized Frequency')
plt.legend()
plt.savefig('{}_angles_hist.png'.format(prefix))

plt.clf()
bins = np.linspace(np.floor(lower_dists-(lower_dists*0.1)), np.ceil(upper_dists+(upper_dists*0.1)), 60)
for i, d in enumerate(distances):
    plt.hist(d, bins, weights=np.ones(len(d)) / len(d), alpha=0.5, histtype=u'stepfilled', edgecolor='k', label='run{}'.format(i))
#plt.hist(distances2, bins, weights=np.ones(len(distances2)) / len(distances2), alpha=0.5, histtype=u'stepfilled', edgecolor='k', label='run2')
#plt.hist(distances3, bins, weights=np.ones(len(distances3)) / len(distances3), alpha=0.5, histtype=u'stepfilled', edgecolor='k', label='run3')
#plt.hist(distances4, bins, weights=np.ones(len(distances4)) / len(distances3), alpha=0.5, histtype=u'stepfilled', edgecolor='k', label='run4')
plt.xlabel('distance(deg)')
plt.ylabel('Normalized Frequency')
plt.legend()
plt.savefig('{}_distances_hist.png'.format(prefix))

#stats
angle_means = [np.mean(a) for a in angles]
#angle_mean1 = np.mean(angles1)
#angle_mean2 = np.mean(angles2)
#angle_mean3 = np.mean(angles3)
#angle_mean4 = np.mean(angles4)

distance_means = [np.mean(d) for d in distances]
#distance_mean1 = np.mean(distances1)
#distance_mean2 = np.mean(distances2)
#distance_mean3 = np.mean(distances3)
#distance_mean4 = np.mean(distances4)

angle_vars = [np.var(a) for a in angles]
#angle_var1 = np.var(angles1)
#angle_var2 = np.var(angles2)
#angle_var3 = np.var(angles3)
#angle_var4 = np.var(angles4)

distance_vars = [np.var(d) for d in distances]
#distance_var1 = np.var(distances1)
#distance_var2 = np.var(distances2)
#distance_var3 = np.var(distances3)
#distance_var4 = np.var(distances4)

print('run:\t', end='')
for i in range(1, len(angles)+1):
    print('{}\t'.format(i), end='')
print()
#print('run:\t1\t2\t3\t')
print('ang:\t', end='')
for am in angle_means:
    print('{:.4}\t'.format(am), end='')
print()
#print('ang:\t{:.4}\t{:.4}\t{:.4}\t{:.4}'.format(angle_mean1, angle_mean2, angle_mean3, angle_mean4))
print('dist:\t', end='')
for dm in distance_means:
    print('{:.4}\t'.format(dm), end='')
print()
#print('dist:\t{:.4}\t{:.4}\t{:.4}'.format(distance_mean1, distance_mean2, distance_mean3))
#print('dist:\t{:.4}\t{:.4}\t{:.4}\t{:.4}'.format(distance_mean1, distance_mean2, distance_mean3, distance_mean4))
print()
print('ang_v:\t', end='')
for av in angle_vars:
    print('{:.4}\t'.format(av), end='')
print()
#print('ang_v\t{:.4}\t{:.4}\t{:.4}'.format(angle_var1, angle_var2, angle_var3))
#print('ang_v\t{:.4}\t{:.4}\t{:.4}\t{:.4}'.format(angle_var1, angle_var2, angle_var3, distance_var4))
print('dist_v:\t', end='')
for dv in distance_vars:
    print('{:.4}\t'.format(dv), end='')
print()
#print('dist_v\t{:.4}\t{:.4}\t{:.4}'.format(distance_var1, distance_var2, distance_var3))
#print('dist_v\t{:.4}\t{:.4}\t{:.4}\t{:.4}'.format(distance_var1, distance_var2, distance_var3, distance_var4))

#equipartition theorem: K(x-x0)^2 = 1/2KbT
#variance = sum((x-x0)^2)/(n-1)
Kas = [0.5/av for av in angle_vars]
Kds = [0.5/dv for dv in distance_vars]
#K1a = 0.5/angle_var1
#K1d = 0.5/distance_var1
#K2a = 0.5/angle_var2
#K2d = 0.5/distance_var2
#K3a = 0.5/angle_var3
#K3d = 0.5/distance_var3
#K4d = 0.5/distance_var4
#K4a = 0.5/angle_var4

#spring constants (units of [F]/[L]*B)
print()
print('ang_k:\t', end='')
for ak in Kas:
    print('{:.4}\t'.format(ak), end='')
print()
print('dist_k:\t', end='')
for dk in Kds:
    print('{:.4}\t'.format(dk), end='')
print()
#print('ang_k:\t{:.4}\t{:.4}\t{:.4}'.format(K1a, K2a, K3a))#K4a))
#print('dist_k:\t{:.4}\t{:.4}\t{:.4}'.format(K1d, K2d, K3d))#K4d))
