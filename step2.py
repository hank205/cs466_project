# CS466 Mini Project (Bioinformatics)
# Edward, Hank, William
# Step 2

import numpy as np
import sys
import time

# Seed
np.random.seed(seed=1)

def read_sequence_file():
    sequence_list = []
    f = open('./data_sets/%s/sequences.fa' % data_folder, 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        sequence_list += [line.strip()]
    f.close()
    return sequence_list
  
def initialize_random_positions():
    '''
    Initializes a random starting position for each sequence.
    '''
    sequence_length = len(sequence_list[0])
    # The highest possible starting index.
    highest_index = sequence_length - ml

    return [np.random.choice(range(highest_index + 1)
        ) for sequence in sequence_list]
    
def gibbs():
    '''
    Performs Gibbs sampling on the list of sequences, given a motif length.
    '''
    curr_position_list = initialize_random_positions()

    # Repeat until position_list doesn't change.
    prev_position_list = []

    while curr_position_list != prev_position_list:
        prev_position_list = curr_position_list[:]
        pwm = []
        for i in curr_position_list:
            # Build a profile Q using sequences in curr_position_list, except i.
            build_profile_matrix(sequence_list, ml)



  
def main():
    if len(sys.argv) != 2:
        print('Usage: python %s data_folder' % (sys.argv[0]))
        exit()
    global data_folder, sequence_list, ml
    data_folder = sys.argv[1]
  
    # set up background matrix
    with open('./data_sets/%s/motiflength.txt' % data_folder, 'r') as f:
        ml = int(f.read())
        
    # background_matrix = []
    # for x in xrange(0, ml):
    #     l = [0.25, 0.25, 0.25, 0.25]
    #     background_matrix.append(l)    
    
    # Read sequence file and run Gibbs sampler.
    sequence_list = read_sequence_file()
    gibbs()
  

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))