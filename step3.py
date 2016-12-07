from __future__ import print_function
import numpy as np
import sys
import time
import math
import matplotlib.pyplot as plt
import os

dirs = ['default', 'a_icpc1', 'a_icpc1.5', 'b_ml6', 'b_ml7', 'c_sc20', 'c_sc5']

def relative_entropy(motif, predictedmotif, key):    
    relative_entropy = 0
    whitening_factor = 0.00001

    for row in range(len(motif[key])):
        for col in range(4):
            # print(motif)
            # print(key, motif[key])
            # print(key, row, motif[key][row])
            # print()
            P = float(motif[key][row][col])
            Q = float(predictedmotif[key][row][col])
            relative_entropy += P * math.log((P+whitening_factor)/(Q+whitening_factor), 2)

    return relative_entropy


def num_overlaps(sites, predictedsites, motif_lengths, key, num):
    '''
    return number of overlaps between sites & predictedsites
    generate graphs visualize the overlaps
    '''
    overlap_sites = 0
    
    site = sites[key]
    predictedsite = predictedsites[key]
    ml = motif_lengths[key]

    site = sites[key]
    predictedsite = predictedsites[key]
    motif_length = motif_lengths[key]

    for i, j in zip(site, predictedsite):
        if i <= j and j < i + motif_length:
            overlap_sites+=1
        elif j <= i and i < j + motif_length:
            overlap_sites+=1

    # site.sort()
    # predictedsite.sort()

    # # print 'site', len(site), site
    # # print 'predictedsite', len(predictedsite), predictedsite
    # # print()

    # pred_idx = 0
    # site_idx = 0

    # while pred_idx < len(predictedsite) or site_idx < len(site):
    #     # print 'pred_idx',pred_idx, 'site_idx',site_idx

    #     if predictedsite[pred_idx]+ml > site[site_idx] and predictedsite[pred_idx] <= site[site_idx]:
    #         overlap_sites += 1
    #         # print 'predictedsite',predictedsite[pred_idx], 'site',site[site_idx]
            
    #     elif site[site_idx]+ml > predictedsite[pred_idx] and site[site_idx] <= predictedsite[pred_idx]:
    #         overlap_sites += 1
    #         # print 'predictedsite',predictedsite[pred_idx], 'site',site[site_idx]

    #     # at least one reaches last site
    #     if pred_idx == len(predictedsite)-1 and site_idx < len(site)-1:
    #         site_idx += 1
    #     elif pred_idx < len(predictedsite)-1 and site_idx == len(site)-1:
    #         pred_idx += 1
    #     elif pred_idx == len(predictedsite)-1 and site_idx == len(site)-1:
    #         break

    #     # none reaches last site
    #     else:
    #         if predictedsite[pred_idx] < site[site_idx]:
    #             pred_idx += 1
    #         elif predictedsite[pred_idx] > site[site_idx]:
    #             site_idx += 1
    #         else:
    #             if predictedsite[pred_idx+1] < site[site_idx+1]:
    #                 pred_idx += 1
    #             else:
    #                 site_idx += 1    
      
    # # generate graph
    # plt.figure()
    # for x in site:
    #     plt.hlines(0.5, x, x+ml, 'r', alpha = 0.7, lw=4)
    # for x in predictedsite:
    #     plt.hlines(0.48, x, x+ml, 'b', alpha = 0.7, lw=4)    
    
    # plt.xlim(0, 600)
    # # plt.xlim(180, 230)
    # plt.ylim(0, 1)
    # # plt.show()
    
    # if not os.path.exists('./evaluations/{}/overlaps'.format(str(num))):
    #     os.makedirs('./evaluations/{}/overlaps'.format(str(num)))
    # plt.savefig('./evaluations/{}/overlaps/{}_overlap.png'.format(str(num), key), bbox_inches='tight')

    return overlap_sites


def read_motifs(num):
    '''
    read the 7 motifs of the num'th dataset
    '''
    motifs = {}
    for data_folder in dirs:
        with open('./data_sets/{}/{}/motif.txt'.format(data_folder, str(num)), 'r') as f:
            motifs[data_folder] = []

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

                motifs[data_folder].append(l)
    return motifs

def read_predictedmotifs(num):
    '''
    read the 7 predicted motifs of the num'th dataset
    '''
    predictedmotifs = {}
    for data_folder in dirs:
        with open('./outcomes/{}/{}/predictedmotif.txt'.format(data_folder, str(num)), 'r') as f:
            predictedmotifs[data_folder] = []

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

                predictedmotifs[data_folder].append(l)
    return predictedmotifs





def evaluate_relative_entropy(motifs, predictedmotifs, num):
    '''
    evaluate relative_entropy
    '''
    global tot_rel_entropy

    rel_entropy = {}
    for key in dirs:
        rel_entropy[key] = relative_entropy(motifs, predictedmotifs, key)
    
    # save sums for evaluating average values
    for data_folder in dirs:
        tot_rel_entropy[data_folder] += rel_entropy[data_folder]

    plt.plot(1, rel_entropy['a_icpc1'], 'ro')
    plt.plot(1.5,rel_entropy['a_icpc1.5'], 'ro')
    plt.plot(2, rel_entropy['default'], 'ro')
    plt.axis([0, 3.0, 0.0, 150.0])
    plt.ylabel("relative entropy")
    plt.xlabel("icpc")
    plt.title('relative entropy vs icpc')
    plt.savefig('./evaluations/{}/relative_entropy_vs_icpc.png'.format(str(num)), bbox_inches='tight')

    plt.figure()
    plt.plot(6, rel_entropy['b_ml6'], 'ro')
    plt.plot(7,rel_entropy['b_ml7'], 'ro')
    plt.plot(8, rel_entropy['default'], 'ro')
    plt.axis([0, 15.0, 0.0, 150.0])
    plt.ylabel("relative entropy")
    plt.xlabel("ml")
    plt.title('relative entropy vs ml')
    plt.savefig('./evaluations/{}/relative_entropy_vs_ml.png'.format(str(num)), bbox_inches='tight')

    plt.figure()
    plt.plot(5, rel_entropy['c_sc5'], 'ro')
    plt.plot(20,rel_entropy['c_sc20'], 'ro')
    plt.plot(10, rel_entropy['default'], 'ro')
    plt.axis([0, 25.0, 0.0, 150.0])
    plt.ylabel("relative entropy")
    plt.xlabel("sc")
    plt.title('relative entropy vs sc')
    plt.savefig('./evaluations/{}/relative_entropy_vs_sc.png'.format(str(num)), bbox_inches='tight')
    plt.close('all')
    # for data_folder in dirs:
    #     print (rel_entropy[data_folder], end=' ')
    # print()

def evaluate_avg_relative_entropy(num):
    '''
    evaluate average relative_entropy of the 10 datasets
    '''
    global tot_rel_entropy

    fig = plt.figure()
    plt.plot(1, tot_rel_entropy['a_icpc1']/num, 'ro')
    plt.plot(1.5,tot_rel_entropy['a_icpc1.5']/num, 'ro')
    plt.plot(2, tot_rel_entropy['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([1, 1.5, 2],[tot_rel_entropy['a_icpc1']/num,tot_rel_entropy['a_icpc1.5']/num,tot_rel_entropy['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))    
    plt.axis([0, 3.0, 0.0, max(tot_rel_entropy['a_icpc1']/num, tot_rel_entropy['a_icpc1.5']/num, tot_rel_entropy['default']/num)+20.0])
    plt.ylabel("relative entropy")
    plt.xlabel("icpc")
    plt.title('average relative entropy vs icpc')
    plt.savefig('./evaluations/avg_relative_entropy_vs_icpc.png', bbox_inches='tight')

    fig = plt.figure()
    plt.plot(6, tot_rel_entropy['b_ml6']/num, 'ro')
    plt.plot(7,tot_rel_entropy['b_ml7']/num, 'ro')
    plt.plot(8, tot_rel_entropy['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([6, 7, 8],[tot_rel_entropy['b_ml6']/num,tot_rel_entropy['b_ml7']/num,tot_rel_entropy['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 15.0, 0.0, max(tot_rel_entropy['b_ml6']/num, tot_rel_entropy['b_ml7']/num, tot_rel_entropy['default']/num)+20.0])
    plt.ylabel("relative entropy")
    plt.xlabel("ml")
    plt.title('average relative entropy vs ml')
    plt.savefig('./evaluations/avg_relative_entropy_vs_ml.png', bbox_inches='tight')

    fig = plt.figure()
    plt.plot(5, tot_rel_entropy['c_sc5']/num, 'ro')
    plt.plot(20,tot_rel_entropy['c_sc20']/num, 'ro')
    plt.plot(10, tot_rel_entropy['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([5, 20, 10],[tot_rel_entropy['c_sc5']/num,tot_rel_entropy['c_sc20']/num,tot_rel_entropy['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 25.0, 0.0, max(tot_rel_entropy['c_sc5']/num, tot_rel_entropy['c_sc20']/num, tot_rel_entropy['default']/num)+20.0])
    plt.ylabel("relative entropy")
    plt.xlabel("sc")
    plt.title('average relative entropy vs sc')
    plt.savefig('./evaluations/avg_relative_entropy_vs_sc.png', bbox_inches='tight')
    plt.close('all')




def evaluate_overlaps(sites, predictedsites, motif_lengths, num):
    '''
    evaluate overlaps
    '''
    global tot_overlaps
    overlaps = {}
    # calculate number of overlaps & visualize the overlaps
    for key in dirs:
        overlaps[key] = num_overlaps(sites, predictedsites, motif_lengths, key, num)
    # for key in dirs:
    #     print 'num_overlaps', key, overlaps[key]
    
    # save sums for evaluating average values
    for data_folder in dirs:
        tot_overlaps[data_folder] += overlaps[data_folder]

    plt.figure()
    plt.plot(1, overlaps['a_icpc1'], 'ro')
    plt.plot(1.5,overlaps['a_icpc1.5'], 'ro')
    plt.plot(2, overlaps['default'], 'ro')
    plt.axis([0, 3.0, 0.0, 15.0])
    plt.ylabel("overlaps")
    plt.xlabel("icpc")
    plt.title('overlaps vs icpc')
    plt.savefig('./evaluations/{}/overlaps_vs_icpc.png'.format(str(num)), bbox_inches='tight')

    plt.figure()
    plt.plot(6, overlaps['b_ml6'], 'ro')
    plt.plot(7,overlaps['b_ml7'], 'ro')
    plt.plot(8, overlaps['default'], 'ro')
    plt.axis([0, 15.0, 0.0, 15.0])
    plt.ylabel("overlaps")
    plt.xlabel("ml")
    plt.title('overlaps vs ml')
    plt.savefig('./evaluations/{}/overlaps_vs_ml.png'.format(str(num)), bbox_inches='tight')

    plt.figure()
    plt.plot(5, overlaps['c_sc5'], 'ro')
    plt.plot(20,overlaps['c_sc20'], 'ro')
    plt.plot(10, overlaps['default'], 'ro')
    plt.axis([0, 25.0, 0.0, 25.0])
    plt.ylabel("overlaps")
    plt.xlabel("sc")
    plt.title('overlaps vs sc')
    plt.savefig('./evaluations/{}/overlaps_vs_sc.png'.format(str(num)), bbox_inches='tight')
    plt.close('all')

    # for data_folder in dirs:
    #     print (overlaps[data_folder], end=' ')
    # print()

def evaluate_avg_overlaps(num):
    '''
    evaluate average overlaps
    '''
    global tot_overlaps

    fig = plt.figure()
    plt.plot(1, tot_overlaps['a_icpc1']/num, 'ro')
    plt.plot(1.5,tot_overlaps['a_icpc1.5']/num, 'ro')
    plt.plot(2, tot_overlaps['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([1, 1.5, 2],[tot_overlaps['a_icpc1']/num,tot_overlaps['a_icpc1.5']/num,tot_overlaps['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 3.0, 0.0, max(tot_overlaps['a_icpc1']/num, tot_overlaps['a_icpc1.5']/num, tot_overlaps['default']/num)+3.0])
    plt.ylabel("overlaps")
    plt.xlabel("icpc")
    plt.title('average overlaps vs icpc')
    plt.savefig('./evaluations/avg_overlaps_vs_icpc.png', bbox_inches='tight')

    fig = plt.figure()
    plt.plot(6, tot_overlaps['b_ml6']/num, 'ro')
    plt.plot(7,tot_overlaps['b_ml7']/num, 'ro')
    plt.plot(8, tot_overlaps['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([6, 7, 8],[tot_overlaps['b_ml6']/num,tot_overlaps['b_ml7']/num,tot_overlaps['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 15.0, 0.0, max(tot_overlaps['b_ml6']/num, tot_overlaps['b_ml7']/num, tot_overlaps['default']/num)+3.0])
    plt.ylabel("overlaps")
    plt.xlabel("ml")
    plt.title('average overlaps vs ml')
    plt.savefig('./evaluations/avg_overlaps_vs_ml.png', bbox_inches='tight')

    fig = plt.figure()
    plt.plot(5, tot_overlaps['c_sc5']/num, 'ro')
    plt.plot(20,tot_overlaps['c_sc20']/num, 'ro')
    plt.plot(10, tot_overlaps['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([5, 20, 10],[tot_overlaps['c_sc5']/num,tot_overlaps['c_sc20']/num,tot_overlaps['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 25.0, 0.0, max(tot_overlaps['c_sc5']/num, tot_overlaps['c_sc20']/num, tot_overlaps['default']/num)+3.0])
    plt.ylabel("overlaps")
    plt.xlabel("sc")
    plt.title('average overlaps vs sc')
    plt.savefig('./evaluations/avg_overlaps_vs_sc.png', bbox_inches='tight')
    plt.close('all')


def evaluate_runtime(runtime, num):
    '''
    evaluate runtime
    '''
    global tot_runtime
    for data_folder in dirs:
        tot_runtime[data_folder] += runtime[data_folder]

    fig = plt.figure()
    plt.plot(1, runtime['a_icpc1'], 'ro')        
    plt.plot(1.5,runtime['a_icpc1.5'], 'ro')
    plt.plot(2, runtime['default'], 'ro')
    # ax = fig.add_subplot(111)
    # for i,j in zip([1, 1.5, 2],[runtime['a_icpc1'],runtime['a_icpc1.5'],runtime['default']]):
    #     ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 3.0, 0.0, max(runtime['a_icpc1'], runtime['a_icpc1.5'], runtime['default'])+1.0])
    plt.ylabel("running time (secs)")
    plt.xlabel("icpc")
    plt.title('running time vs icpc')
    plt.savefig('./evaluations/{}/runtime_vs_icpc.png'.format(str(num)), bbox_inches='tight')

    fig = plt.figure()
    plt.plot(6, runtime['b_ml6'], 'ro')
    plt.plot(7,runtime['b_ml7'], 'ro')
    plt.plot(8, runtime['default'], 'ro')
    # ax = fig.add_subplot(111)
    # for i,j in zip([6, 7, 8],[runtime['b_ml6'],runtime['b_ml7'],runtime['default']]):
    #     ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 15.0, 0.0, max(runtime['b_ml6'], runtime['b_ml7'], runtime['default'])+1.0])
    plt.ylabel("running time (secs)")
    plt.xlabel("ml")
    plt.title('running time vs ml')
    plt.savefig('./evaluations/{}/runtime_vs_ml.png'.format(str(num)), bbox_inches='tight')

    fig = plt.figure()
    plt.plot(5, runtime['c_sc5'], 'ro')
    plt.plot(20,runtime['c_sc20'], 'ro')
    plt.plot(10, runtime['default'], 'ro')
    # ax = fig.add_subplot(111)
    # for i,j in zip([5, 20, 10],[runtime['c_sc5'],runtime['c_sc20'],runtime['default']]):
    #     ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 25.0, 0.0, max(runtime['c_sc5'], runtime['c_sc20'], runtime['default'])+1.0])
    plt.ylabel("running time (secs)")
    plt.xlabel("sc")
    plt.title('running time vs sc')
    plt.savefig('./evaluations/{}/runtime_vs_sc.png'.format(str(num)), bbox_inches='tight')
    plt.close('all')

    # for data_folder in dirs:
    #     print (runtime[data_folder], end=' ')
    # print()


def evaluate_avg_runtime(num):
    '''
    evaluate average runtime
    '''
    global tot_runtime

    fig = plt.figure()
    plt.plot(1, tot_runtime['a_icpc1']/num, 'ro')
    plt.plot(1.5,tot_runtime['a_icpc1.5']/num, 'ro')
    plt.plot(2, tot_runtime['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([1, 1.5, 2],[tot_runtime['a_icpc1']/num,tot_runtime['a_icpc1.5']/num,tot_runtime['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 3.0, 0.0, max(tot_runtime['a_icpc1']/num, tot_runtime['a_icpc1.5']/num, tot_runtime['default']/num)+1.0])
    plt.ylabel("running time (secs)")
    plt.xlabel("icpc")
    plt.title('average running time vs icpc')
    plt.savefig('./evaluations/avg_runtime_vs_icpc.png', bbox_inches='tight')

    fig = plt.figure()
    plt.plot(6, tot_runtime['b_ml6']/num, 'ro')
    plt.plot(7,tot_runtime['b_ml7']/num, 'ro')
    plt.plot(8, tot_runtime['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([6, 7, 8],[tot_runtime['b_ml6']/num,tot_runtime['b_ml7']/num,tot_runtime['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 15.0, 0.0, max(tot_runtime['b_ml6']/num, tot_runtime['b_ml7']/num, tot_runtime['default']/num)+1.0])
    plt.ylabel("running time (secs)")
    plt.xlabel("ml")
    plt.title('average running time vs ml')
    plt.savefig('./evaluations/avg_runtime_vs_ml.png', bbox_inches='tight')

    fig = plt.figure()
    plt.plot(5, tot_runtime['c_sc5']/num, 'ro')
    plt.plot(20,tot_runtime['c_sc20']/num, 'ro')
    plt.plot(10, tot_runtime['default']/num, 'ro')
    ax = fig.add_subplot(111)
    for i,j in zip([5, 20, 10],[tot_runtime['c_sc5']/num,tot_runtime['c_sc20']/num,tot_runtime['default']/num]):
        ax.annotate(str("{0:.2f}".format(j)),xy=(i,j+0.02))
    plt.axis([0, 25.0, 0.0, max(tot_runtime['c_sc5']/num, tot_runtime['c_sc20']/num, tot_runtime['default']/num)+1.0])
    plt.ylabel("running time (secs)")
    plt.xlabel("sc")
    plt.title('average running time vs sc')
    plt.savefig('./evaluations/avg_runtime_vs_sc.png', bbox_inches='tight')
    plt.close('all')


def main():

    if not os.path.exists('./evaluations'):
        os.makedirs('./evaluations')

    global tot_rel_entropy
    tot_rel_entropy = {}
    for data_folder in dirs:
        tot_rel_entropy[data_folder] = 0.0

    global tot_overlaps
    tot_overlaps = {}
    for data_folder in dirs:
        tot_overlaps[data_folder] = 0.0

    global tot_runtime
    tot_runtime = {}
    for data_folder in dirs:
        tot_runtime[data_folder] = 0.0

    trials = 10
    for num in range(trials):
        print('num', num)

        if not os.path.exists('./evaluations/{}'.format(str(num))):
            os.makedirs('./evaluations/{}'.format(str(num)))

        # evaluate relative_entropy
        motifs = read_motifs(num)
        predictedmotifs = read_predictedmotifs(num)
        # for key in dirs:
        #     print key
        #     for l in motif[key]:
        #         print l
        # print()        
        # for key in dirs:
        #     print key
        #     for l in predictedmotif[key]:
        #         print l

        evaluate_relative_entropy(motifs, predictedmotifs, num)

    #############################################################################################
        
        # evaluate overlaps
        sites = {}
        predictedsites = {}
        motif_lengths = {}

        for data_folder in dirs:
            with open('./data_sets/{}/{}/sites.txt'.format(data_folder, str(num)), 'r') as f:
                sites[data_folder] = []
                for l in f:
                    sites[data_folder].append(int(l))

        for data_folder in dirs:
            with open('./outcomes/{}/{}/predictedsites.txt'.format(data_folder, str(num)), 'r') as f:
                predictedsites[data_folder] = []
                for l in f:
                    predictedsites[data_folder].append(int(l))
        
        for data_folder in dirs:
            with open('./data_sets/{}/{}/motiflength.txt'.format(data_folder, str(num)), 'r') as f:
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

        evaluate_overlaps(sites, predictedsites, motif_lengths, num)

    #############################################################################################

        # evaluate runtime
        runtime = {}

        for data_folder in dirs:
            with open('./outcomes/{}/{}/runtime.txt'.format(data_folder, str(num)), 'r') as f:
                runtime[data_folder] = float(f.read())
        # for key in dirs:
        #     print key, runtime[key]

        evaluate_runtime(runtime, num)
   

    # evaluate average values
    evaluate_avg_relative_entropy(trials)
    evaluate_avg_overlaps(trials)
    evaluate_avg_runtime(trials)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))