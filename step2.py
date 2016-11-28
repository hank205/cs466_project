# CS466 Mini Project (Bioinformatics)
# Edward, Hank, William
# Step 2

import numpy as np
import sys
import time
import os

# Seed
np.random.seed(seed=1)
base_list = ['A', 'C', 'G', 'T']

def read_sequence_file(num):
    sequence_list = []
    f = open('./data_sets/{}/{}/sequences.fa'.format(data_folder, str(num)), 'r')
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
    
def build_profile_matrix(position_list, i):
    '''
    Build a PWM for every sequence except i.
    '''
    # Initialize weight matrix.
    pwm = [np.zeros(4) for j in range(ml)]
    for sequence_i, sequence in enumerate(sequence_list):
        # Skip the ith sequence.
        if sequence_i == i:
            continue
        position = position_list[sequence_i]
        possible_motif = sequence[position:position+ml]
        for base_i, base in enumerate(possible_motif):
            pwm[base_i][base_list.index(base)] += 1
    # Normalize each sequence position's probability.
    pwm = np.array([row / np.sum(row) for row in pwm])
    return pwm

def string_score(subsequence, pwm):
    '''
    Given a subsequence, compute the score according to the profile matrix.
    '''
    score_total = 0.0
    for base_i, base in enumerate(subsequence):
        base_prob = pwm[base_i][base_list.index(base)]
        if base_prob == 0:
            continue
        score_total += np.log(pwm[base_i][base_list.index(base)] / 0.25)
    return score_total

def find_best_index(sequence, pwm):
    '''
    Given a sequence and a PWM, find the index where the PWM best matches the 
    sequence.
    '''
    best_index, best_index_score = -1, -float('inf')
    # Loop over possible starting indices.
    for x_i in range(len(sequence) - ml + 1):
        possible_motif = sequence[x_i:x_i+ml]
        current_score = string_score(possible_motif, pwm)
        if current_score > best_index_score:
            best_index_score = current_score
            best_index = x_i
    return best_index

def gibbs():
    '''
    Performs Gibbs sampling on the list of sequences, given a motif length.
    '''
    curr_position_list = initialize_random_positions()

    # Repeat until position_list doesn't change.
    prev_position_list = []

    while curr_position_list != prev_position_list:
        prev_position_list = curr_position_list[:]
        for sequence_i, position in enumerate(curr_position_list):
            # Build a profile Q using sequences in curr_position_list, except i.
            pwm = build_profile_matrix(curr_position_list, sequence_i)
            # Find where the profile matches best in sequence_i.
            curr_position_list[sequence_i] = find_best_index(sequence_list[
                sequence_i], pwm)
    return curr_position_list
  
def main():
    if len(sys.argv) != 2:
        print('Usage: python %s data_folder' % (sys.argv[0]))
        exit()
    global data_folder, sequence_list, ml
    data_folder = sys.argv[1]
  

    if not os.path.exists('outcomes/{}'.format(data_folder)):
        os.makedirs('outcomes/{}'.format(data_folder))

    for num in range(10):
        run_start = time.time()

        # set up background matrix
        with open('./data_sets/{}/{}/motiflength.txt'.format(data_folder, str(num)), 'r') as f:
            ml = int(f.read())

        # Read sequence file and run Gibbs sampler.
        sequence_list = read_sequence_file(num)
        curr_position_list = gibbs()


        # Make PWM predictedmotif out of the curr_position_list
        predictedmotif = [[0.0, 0.0, 0.0, 0.0] for x in range(ml)]
        for i in range(len(curr_position_list)):        
            start_pos = curr_position_list[i]
            motif = sequence_list[i][start_pos:start_pos+ml]
            for j in range(ml):
                char = motif[j]
                predictedmotif[j][base_list.index(char)] += 1
        
        # Normalize predictedmotif
        sc = len(sequence_list)
        for i in range(ml):
            predictedmotif[i] = map(lambda x: x/sc, predictedmotif[i])

        # Generate files
        if not os.path.exists('outcomes/{}/{}'.format(data_folder, str(num))):
            os.makedirs('outcomes/{}/{}'.format(data_folder, str(num)))

        with open('outcomes/{}/{}/predictedmotif.txt'.format(data_folder, str(num)), 'w') as f:
            f.write('>MOTIF ' + str(ml) + ' ' + data_folder + '\n')
            for l in predictedmotif:
                for s in l:
                    f.write(str(s) + ' ')
                f.write('\n')
            f.write('<')

        with open('outcomes/{}/{}/predictedsites.txt'.format(data_folder, str(num)), 'w') as f:
            for l in curr_position_list:
                f.write(str(l) + '\n')

        runtime = (time.time() - run_start)
        with open('outcomes/{}/{}/runtime.txt'.format(data_folder, str(num)), 'w') as f:
            f.write(str(runtime))


if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
