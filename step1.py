# CS466 Mini Project (Bioinformatics)
# Edward, Hank, William
# Step 1

from __future__ import print_function
import numpy as np
import sys
import time
import os

base_list = ['A', 'C', 'G', 'T']
# Key: icpc -> (float). Value: p -> (float). To be used for generate_motif()
icpc_to_p_dct = {1.0:0.8105, 1.5:0.9245, 2.0:1.0}

def generate_random_sequences():
    '''
    Generates SC random sequences with uniform probability. Each random sequence
    has length SL.
    '''
    sequence_list = []
    for i in range(sc):
      sequence_list += [np.random.choice(base_list, sl)]
    return sequence_list

def generate_motif():
    '''
    Generates a random motif with length ML. Total information content = ICPC *
    ML.
    '''
    # Get the value of p corresponding to the ICPC.
    p = icpc_to_p_dct[icpc]
    # Initialize position weight matrix to all (1 - p) / 3.
    pwm = np.ones([ml, 4]) * (1 - p) / 3
    # Assume uniform background distribution, log base 2.
    for i in range(ml):
      # Preferred nucleotide gets probability p. Others stay the same.
      preferred_nucleotide = np.random.choice(4)
      pwm[i][preferred_nucleotide] = p
    return pwm
  
def plant_binding_sites(sequence_list, pwm):
    '''
    Mutates the sequence list in-place.
    '''
    plant_site_list = []
    # print('motifs:')
    for current_sequence in sequence_list:
        # Pick an index to plant.
        plant_index = np.random.choice(sl - ml + 1)
        plant_site_list += [plant_index]

        for prob_list in pwm:
            # Generate a nucleotide from the motif PWM.
            motif_nucleotide = base_list[np.random.choice(4, p=prob_list)]
            # Plant the nucleotide at the current index.
            current_sequence[plant_index] = motif_nucleotide            
            plant_index += 1
            # print(motif_nucleotide, end='')
        # print('\r')
    return plant_site_list

def main():
    if len(sys.argv) != 6:
      print('Usage: python %s ICPC ML SL SC data_folder' % (sys.argv[0]))
      exit()
    global icpc, ml, sl, sc
    icpc = float(sys.argv[1])
    ml, sl, sc = map(int, sys.argv[2:-1])
    data_folder = sys.argv[-1]
    # print (ml, sl, sc,data_folder)
  
    sequence_list = generate_random_sequences()
    pwm = generate_motif()
    # print('original sequence_list:')
    # for l in sequence_list:
    #     print(''.join(l))
    # print('pwm:')    
    # for l in pwm:
    #     print(l)        
    plant_site_list = plant_binding_sites(sequence_list, pwm)
    # print('plant site indices:')
    # for i in plant_site_list:
    #     print(''.join(str(i)))
    # print('planted sequence_list:')
    # for l in sequence_list:
    #     print(''.join(l))

    if not os.path.exists('./data_sets'):
        os.makedirs('./data_sets')

    if not os.path.exists('./data_sets/{}'.format(data_folder)):
        os.makedirs('./data_sets/{}'.format(data_folder))

    # generate 10 data sets for input configuration
    for num in range(10):
        if not os.path.exists('./data_sets/{}/{}'.format(data_folder, str(num))):
            os.makedirs('./data_sets/{}/{}'.format(data_folder, str(num)))
        # generate files
        with open('./data_sets/{}/{}/sequences.fa'.format(data_folder, str(num)), 'w') as f:
            f.write('>SEQUENCES\n')
            for l in sequence_list:
                f.write(''.join(l) + '\n')

        with open('./data_sets/{}/{}/sites.txt'.format(data_folder, str(num)), 'w') as f:
            for x in plant_site_list:
                f.write(str(x) + '\n')
        
        with open('./data_sets/{}/{}/motif.txt'.format(data_folder, str(num)), 'w') as f:
            f.write('>MOTIF1 ' + str(ml) + '\n')
            for l in pwm:
                for s in l:
                    f.write(str(s) + ' ')
                f.write('\n')
            f.write('<')
        
        with open('./data_sets/{}/{}/motiflength.txt'.format(data_folder, str(num)), 'w') as f:
            f.write(str(ml))

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))