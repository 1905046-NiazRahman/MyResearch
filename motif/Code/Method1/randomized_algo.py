from random import randint
import os
import time

def score(motifs):
    score = 0
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        max_freq = max(column.count(nucleotide) for nucleotide in 'ACGT')
        score += len(motifs) - max_freq
    return score

def GenerateProfileMatWithPesudoCount(dna):
    seqNum = float(len(dna))
    nucleotide = ['A', 'C', 'G', 'T']
    matrix = []
    for i in range(len(dna[0])):
        base_i = [seq[i] for seq in dna]
        colProfile = [float(base_i.count(n) + 1) / float(seqNum + 4) for n in nucleotide]
        matrix.append(colProfile)
    return [list(i) for i in zip(*matrix)]

def FindPortableKmer(text, k, matrix):
    maxPr = 0
    prList = []
    for i in range(len(text)-k+1):
        KmerPr = 1
        pattern = text[i:i+k]
        for j in range(len(pattern)):
            profile = GenerateProfileDict(matrix, j)
            KmerPr *= profile[pattern[j]]
        prList.append(KmerPr)
    max_index = prList.index(max(prList))
    maxKmer = text[max_index:max_index+k]
    return maxKmer, max_index  # Return the motif and its starting position

def GenerateProfileDict(matrix, j):
    return {'A': matrix[0][j], 'C': matrix[1][j], 'G': matrix[2][j], 'T': matrix[3][j]}

def RandomizedMotifSearch(dna, k, max_iterations=1000, score_threshold=0.000000001):
    # Store initial random positions
    kmerIndex = [randint(0, len(dna[0]) - k) for _ in range(len(dna))]
    Motifs = [dna[i][j:j+k] for i, j in enumerate(kmerIndex)]
    BestMotifs = Motifs
    BestPositions = kmerIndex
    iterations = 0

    while iterations < max_iterations:
        ProfileMatrix = GenerateProfileMatWithPesudoCount(Motifs)
        Motifs = []
        Positions = []
        for i in range(len(dna)):
            motif, position = FindPortableKmer(dna[i], k, ProfileMatrix)
            Motifs.append(motif)
            Positions.append(position)
        
        if score(Motifs) + score_threshold < score(BestMotifs):
            BestMotifs = Motifs
            BestPositions = Positions
        else:
            break
        
        iterations += 1

    return BestMotifs, BestPositions, score(BestMotifs)  # Also return the score of BestMotifs

input_files = ['hm03.txt', 'yst04r.txt', 'yst08r.txt']
ks = [10, 15, 20]
max_iterations_values = [i*10000 for i in range(1, 6)]
itr_cnt = 100

for input_file in input_files:
    with open(f'inputs/{input_file}', 'r') as file:
        dnaSet = file.read().splitlines()
    
    for k in ks:
        output_dir = f'outputs/randomized/{input_file.split(".")[0]}/motif_len_{k}'
        os.makedirs(output_dir, exist_ok=True)
        
        with open(os.path.join(output_dir, 'motifs_scores.txt'), 'a') as scores_file, \
             open(os.path.join(output_dir, 'timing_info.txt'), 'a') as timing_file:
            
            for max_iterations in max_iterations_values:
                tot_time = 0
                best_motifs = []
                best_positions = []
                best_score = int(1e9)
                for i in range(itr_cnt):
                    start_time = time.time()
                    motifs, positions, motifs_score = RandomizedMotifSearch(dnaSet, k, max_iterations)
                    end_time = time.time()
                    tot_time += end_time - start_time
                    if motifs_score < best_score:
                        best_motifs = motifs
                        best_positions = positions
                        best_score = motifs_score

                with open(os.path.join(output_dir, 'motifs_output.txt'), 'a') as motif_file, \
                     open(os.path.join(output_dir, 'motifs_indices.txt'), 'a') as index_file:
                    motif_file.write(f"Max Iterations: {max_iterations}\n" + "\n".join(best_motifs) + "\n\n")
                    index_file.write(f"Max Iterations: {max_iterations}\n" + " ".join(map(str, best_positions)) + "\n\n")

                scores_file.write(f"Max Iterations: {max_iterations}\nScore: {best_score}\n\n")
                timing_file.write(f"Max Iterations: {max_iterations}\nExecution Time: {tot_time / itr_cnt} seconds\n\n")