import numpy as np
import sys
import time
import math
import matplotlib.pyplot as plt
import os

dirs = ['default', 'a_icpc1', 'a_icpc1.5', 'b_ml6', 'b_ml7', 'c_sc20', 'c_sc5']

def relative_entropy(motif, predictedmotif, key):    
    relative_entropy = 0

    for l in range(len(motif[key])):
        for x in range(4):
            P = float(motif[key][l][x])
            Q = float(predictedmotif[key][l][x])
            relative_entropy += P * math.log(((P+0.00001) / (Q+0.00001)), 2)

    return relative_entropy



def num_overlapping(sites, predictedsites, motif_lengths, key):
    overlap_sites = 0
    
    site = sites[key]
    predictedsite = predictedsites[key]
    ml = motif_lengths[key]

    site.sort()
    predictedsite.sort()

    # print 'site', len(site), site
    # print 'predictedsite', len(predictedsite), predictedsite
    # print()

    pred_idx = 0
    site_idx = 0

    while pred_idx < len(predictedsite) or site_idx < len(site):
        # print 'pred_idx',pred_idx, 'site_idx',site_idx

        if predictedsite[pred_idx]+ml > site[site_idx] and predictedsite[pred_idx] <= site[site_idx]:
            overlap_sites += 1
            # print 'predictedsite',predictedsite[pred_idx], 'site',site[site_idx]
            
        elif site[site_idx]+ml > predictedsite[pred_idx] and site[site_idx] <= predictedsite[pred_idx]:
            overlap_sites += 1
            # print 'predictedsite',predictedsite[pred_idx], 'site',site[site_idx]

        # at least one reaches last site
        if pred_idx == len(predictedsite)-1 and site_idx < len(site)-1:
            site_idx += 1
        elif pred_idx < len(predictedsite)-1 and site_idx == len(site)-1:
            pred_idx += 1
        elif pred_idx == len(predictedsite)-1 and site_idx == len(site)-1:
            break

        # none reaches last site
        else:
            if predictedsite[pred_idx] < site[site_idx]:
                pred_idx += 1
            elif predictedsite[pred_idx] > site[site_idx]:
                site_idx += 1
            else:
                if predictedsite[pred_idx+1] < site[site_idx+1]:
                    pred_idx += 1
                else:
                    site_idx += 1

        
      
    for x in site:
        plt.hlines(0.5, x, x+ml, 'r', alpha = 0.7, lw=4)
    for x in predictedsite:
        plt.hlines(0.48, x, x+ml, 'b', alpha = 0.7, lw=4)    
    
    plt.xlim(0, 600)
    # plt.xlim(180, 230)
    plt.ylim(0, 1)
    # plt.show()
    plt.savefig('./evaluations/overlaps/'+key+'_overlap.png', bbox_inches='tight')

    
    return overlap_sites

def main():

    motif = {}
    predictedmotif = {}

    for data_folder in dirs:
        with open('./data_sets/%s/motif.txt' % data_folder, 'r') as f:
            motif[data_folder] = []

            for i, line in enumerate(f):
                if i == 0:
                    continue
                if line == '<':
                    continue
                
                data = line.split()
                l = []
                l.append(float(data[0]))
                l.append(float(data[1]))
                l.append(float(data[2]))
                l.append(float(data[3]))

                motif[data_folder].append(l)

    for data_folder in dirs:
        with open('./outcomes/%s/predictedmotif.txt' % data_folder, 'r') as f:
            predictedmotif[data_folder] = []

            for i, line in enumerate(f):
                if i == 0:
                    continue
                if line == '<':
                    continue
                
                data = line.split()
                l = []
                l.append(float(data[0]))
                l.append(float(data[1]))
                l.append(float(data[2]))
                l.append(float(data[3]))

                predictedmotif[data_folder].append(l)


    # for key in dirs:
    #     print key
    #     for l in motif[key]:
    #         print l
    # print()        
    # for key in dirs:
    #     print key
    #     for l in predictedmotif[key]:
    #         print l

    # for key in dirs:
    #     print key
    #     print relative_entropy(motif, predictedmotif, key)

    if not os.path.exists('./evaluations'):
        os.makedirs('./evaluations')

    rel_entropy_default = relative_entropy(motif, predictedmotif, 'default')

    rel_entropy_icpc1 = relative_entropy(motif, predictedmotif, 'a_icpc1')
    rel_entropy_icpc1_5 = relative_entropy(motif, predictedmotif, 'a_icpc1.5')
    plt.plot(1, rel_entropy_icpc1, 'ro')
    plt.plot(1.5,rel_entropy_icpc1_5, 'ro')
    plt.plot(2, rel_entropy_default, 'ro')
    plt.axis([0, 3.0, 0.0, 100.0])
    plt.xlabel("icpc")
    plt.ylabel("relative entropy")
    # plt.legend(numpoints=1)
    plt.title('icpc vs relative entropy')
    # plt.show()
    plt.savefig('./evaluations/icpc_vs_relative_entropy.png', bbox_inches='tight')

    plt.figure()
    rel_entropy_b_ml6 = relative_entropy(motif, predictedmotif, 'b_ml6')
    rel_entropy_b_ml7 = relative_entropy(motif, predictedmotif, 'b_ml7')
    plt.plot(6, rel_entropy_b_ml6, 'ro')
    plt.plot(7,rel_entropy_b_ml7, 'ro')
    plt.plot(8, rel_entropy_default, 'ro')
    plt.axis([0, 15.0, 0.0, 100.0])
    plt.xlabel("ml")
    plt.ylabel("relative entropy")
    plt.title('ml vs relative entropy')
    plt.savefig('./evaluations/ml_vs_relative_entropy.png', bbox_inches='tight')

    plt.figure()
    rel_entropy_c_sc5 = relative_entropy(motif, predictedmotif, 'c_sc5')
    rel_entropy_c_sc20 = relative_entropy(motif, predictedmotif, 'c_sc20')
    plt.plot(5, rel_entropy_c_sc5, 'ro')
    plt.plot(20,rel_entropy_c_sc20, 'ro')
    plt.plot(10, rel_entropy_default, 'ro')
    plt.axis([0, 25.0, 0.0, 150.0])
    plt.xlabel("sc")
    plt.ylabel("relative entropy")
    plt.title('sc vs relative entropy')
    plt.savefig('./evaluations/sc_vs_relative_entropy.png', bbox_inches='tight')
#############################################################################################

    sites = {}
    predictedsites = {}

    for data_folder in dirs:
        with open('./data_sets/%s/sites.txt' % data_folder, 'r') as f:
            sites[data_folder] = []
            for l in f:
                sites[data_folder].append(int(l))

    for data_folder in dirs:
        with open('./outcomes/%s/predictedsites.txt' % data_folder, 'r') as f:
            predictedsites[data_folder] = []
            for l in f:
                predictedsites[data_folder].append(int(l))

    motif_lengths = {}
    for data_folder in dirs:
        with open('./data_sets/%s/motiflength.txt' % data_folder, 'r') as f:
            motif_lengths[data_folder] = int(f.read())

    # for key in dirs:
    #     print key
    #     for l in sites[key]:
    #         print l
    # print()        
    # for key in dirs:
    #     print key
    #     for l in predictedsites[key]:
    #         print l

    # for k in motif_lengths:
    #     print k, motif_lengths[k]
    
    overlapping = {}

    if not os.path.exists('./evaluations/overlaps'):
        os.makedirs('./evaluations/overlaps')

    for key in dirs:
        plt.figure()
        overlapping[key] = num_overlapping(sites, predictedsites, motif_lengths, key)

    # for key in dirs:
    #     print 'num_overlapping', key, overlapping[key]

    plt.figure()
    plt.plot(1, overlapping['a_icpc1'], 'ro')
    plt.plot(1.5,overlapping['a_icpc1.5'], 'ro')
    plt.plot(2, overlapping['default'], 'ro')
    plt.axis([0, 3.0, 0.0, 10.0])
    plt.xlabel("icpc")
    plt.ylabel("overlapping")
    plt.title('icpc vs overlapping')
    plt.savefig('./evaluations/icpc_vs_overlapping.png', bbox_inches='tight')

    plt.figure()
    plt.plot(6, overlapping['b_ml6'], 'ro')
    plt.plot(7,overlapping['b_ml7'], 'ro')
    plt.plot(8, overlapping['default'], 'ro')
    plt.axis([0, 15.0, 0.0, 10.0])
    plt.xlabel("ml")
    plt.ylabel("overlapping")
    plt.title('ml vs overlapping')
    plt.savefig('./evaluations/ml_vs_overlapping.png', bbox_inches='tight')

    plt.figure()
    plt.plot(5, overlapping['c_sc5'], 'ro')
    plt.plot(20,overlapping['c_sc20'], 'ro')
    plt.plot(10, overlapping['default'], 'ro')
    plt.axis([0, 25.0, 0.0, 25.0])
    plt.xlabel("sc")
    plt.ylabel("overlapping")
    plt.title('sc vs overlapping')
    plt.savefig('./evaluations/sc_vs_overlapping.png', bbox_inches='tight')
#############################################################################################

    runtime = {}

    for data_folder in dirs:
        with open('./outcomes/%s/runtime.txt' % data_folder, 'r') as f:
            runtime[data_folder] = float(f.read())

    # for key in dirs:
    #     print key, runtime[key]

    plt.figure()
    plt.plot(1, runtime['a_icpc1'], 'ro')
    plt.plot(1.5,runtime['a_icpc1.5'], 'ro')
    plt.plot(2, runtime['default'], 'ro')
    plt.axis([0, 3.0, 0.0, 1.0])
    plt.xlabel("icpc")
    plt.ylabel("runtime")
    plt.title('icpc vs runtime')
    plt.savefig('./evaluations/icpc_vs_runtime.png', bbox_inches='tight')

    plt.figure()
    plt.plot(6, runtime['b_ml6'], 'ro')
    plt.plot(7,runtime['b_ml7'], 'ro')
    plt.plot(8, runtime['default'], 'ro')
    plt.axis([0, 15.0, 0.0, 1.0])
    plt.xlabel("ml")
    plt.ylabel("runtime")
    plt.title('ml vs runtime')
    plt.savefig('./evaluations/ml_vs_runtime.png', bbox_inches='tight')

    plt.figure()
    plt.plot(5, runtime['c_sc5'], 'ro')
    plt.plot(20,runtime['c_sc20'], 'ro')
    plt.plot(10, runtime['default'], 'ro')
    plt.axis([0, 25.0, 0.0, 2.0])
    plt.xlabel("sc")
    plt.ylabel("runtime")
    plt.title('sc vs runtime')
    plt.savefig('./evaluations/sc_vs_runtime.png', bbox_inches='tight')


if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))