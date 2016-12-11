# CS466 Mini Project (Bioinformatics)
# Edward, Hank, William
# Step 2

import numpy as np
import operator
import sys
import time
import os

# Seed
np.random.seed(seed=int(time.time()))
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
    # Normalize each sequence position's probability. Added a smoothing factor (0.5/2)
    pwm = np.array([(row+0.5) / (np.sum(row)+2) for row in pwm])
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
        score_total += np.log(base_prob / background[base])
    return score_total

def find_best_index(sequence, pwm):
    '''
    Given a sequence and a PWM, find the index where the PWM best matches the 
    sequence. Probabilistically choose the index.
    '''
    index_prob_list = []
    # Loop over possible starting indices.
    for x_i in range(len(sequence) - ml + 1):
        possible_motif = sequence[x_i:x_i+ml]
        current_score = string_score(possible_motif, pwm)
        index_prob_list += [current_score]
    # assert len(index_prob_list) == ml
    # Normalize the score list with softmax.
    exp_list = np.exp(index_prob_list)
    index_prob_list = exp_list / np.sum(exp_list)
    best_index = np.random.choice(len(index_prob_list), p=index_prob_list)
    return best_index

def compute_information_content(pwm):
    '''
    Given a profile weight matrix, compute its information content.
    '''
    information_content = 0.0
    for row in pwm:
        for base_prob in row:
            if base_prob == 0:
                continue
            information_content += base_prob * np.log(base_prob / 0.25)
    return information_content

def gibbs():
    '''
    Performs Gibbs sampling on the list of sequences, given a motif length.
    '''
    pos_info_dct = {}
<<<<<<< Updated upstream

    for i in range(100):
        curr_position_list = initialize_random_positions()
        # Repeat a hundred times. Get the position list with the best
        # information content.
        for i in range(1000):
            for sequence_i, position in enumerate(curr_position_list):
                # Build a profile Q using sequences in curr_position_list,
                # except i.
                pwm = build_profile_matrix(curr_position_list, sequence_i)
                # Find where the profile matches best in sequence_i.
                curr_position_list[sequence_i] = find_best_index(sequence_list[
                    sequence_i], pwm)

            key = tuple(curr_position_list[:])
            if key in pos_info_dct:
                continue
            s = build_profile_matrix(curr_position_list, 'dummy_arg')
            pos_info_dct[key] = compute_information_content(s)

=======
    # Repeat a hundred times. Get the position list with the best information
    # content.
    for i in range(10):
        for sequence_i, position in enumerate(curr_position_list):
            # Build a profile Q using sequences in curr_position_list, except i.
            pwm = build_profile_matrix(curr_position_list, sequence_i)
            # Find where the profile matches best in sequence_i.
            curr_position_list[sequence_i] = find_best_index(sequence_list[
                sequence_i], pwm)

        key = tuple(curr_position_list[:])
        if key in pos_info_dct:
            continue
        s = build_profile_matrix(curr_position_list, 'dummy_arg')
        pos_info_dct[key] = compute_information_content(s)
>>>>>>> Stashed changes
    pos_info_dct = sorted(pos_info_dct.items(), key=operator.itemgetter(1),
        reverse=True)

    return pos_info_dct[0][0]

def main():
    if len(sys.argv) != 2:
        print('Usage: python %s data_folder' % (sys.argv[0]))
        exit()
    global data_folder, sequence_list, ml, background
    data_folder = sys.argv[1]
  
    if not os.path.exists('outcomes2/{}'.format(data_folder)):
        os.makedirs('outcomes2/{}'.format(data_folder))

    for num in range(10):
        print num
        run_start = time.time()

        # set up background dict
        with open('./background_ACGT.txt', 'r') as f:
            background = {}
            for (l, c) in zip(f, base_list):
                background[c] = float(l)
        
        # load in motif length
        with open('./data_sets/{}/{}/motiflength.txt'.format(data_folder, str(num)), 'r') as f:
            ml = int(f.read())

        # Read sequence file and run Gibbs sampler.
        sequence_list = read_sequence_file(num)
        curr_position_list = gibbs()

        # Make PWM predictedmotif out of the curr_position_list
        predictedmotif = build_profile_matrix(curr_position_list, 'dummy_arg')

        # Generate files
        if not os.path.exists('outcomes2/{}/{}'.format(data_folder, str(num))):
            os.makedirs('outcomes2/{}/{}'.format(data_folder, str(num)))

        with open('outcomes2/{}/{}/predictedmotif.txt'.format(data_folder, str(num)), 'w') as f:
            # f.write('>MOTIF ' + str(ml) + ' ' + data_folder + '\n')
            f.write('>MOTIF ' + str(ml) + '\n')
            for l in predictedmotif:
                for s in l:
                    f.write(str(s) + ' ')
                f.write('\n')
            f.write('<')

        with open('outcomes2/{}/{}/predictedsites.txt'.format(data_folder, str(num)), 'w') as f:
            for l in curr_position_list:
                f.write(str(l) + '\n')

        runtime = (time.time() - run_start)
        with open('outcomes2/{}/{}/runtime.txt'.format(data_folder, str(num)), 'w') as f:
            f.write(str(runtime))

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
