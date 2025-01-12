from random import randint
import time 
import math
import random

def CalHammingDistance(seq1, seq2):
    assert len(seq1) == len(seq2)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist

def GenerateProfileMatWithPesudoCount(dna):
    seqNum = float(len(dna))
    nucleotide = ['A', 'C', 'G', 'T']
    matrix = []
    for i in range(len(dna[0])):
        base_i = [seq[i] for seq in dna]

        colProfile = [float(base_i.count(n) + 1)/float(seqNum + 4) for n in nucleotide]

        matrix.append(colProfile)
    return [list(i) for i in zip(*matrix)]



def FindPortableKmer(text, k, matrix):

    prList = []
    for i in range(len(text)-k+1):
        KmerPr = 1
        pattern = text[i:i+k]

        for j in range(len(pattern)):
            profile = GenerateProfileDict(matrix, j)

            KmerPr *= profile[pattern[j]]

        prList.append(KmerPr)
    i = prList.index(max(prList))
    maxKmer = text[i:i+k]
    return maxKmer, i


def GenerateProfileDict(matrix, j):
    return {'A': matrix[0][j], 'C': matrix[1][j], 'G': matrix[2][j], 'T': matrix[3][j]}


def GenerateProfileDictAlt(matrix, j):
    nucleotides = ['A', 'C', 'G', 'T']
    profile_dict = {}
    for idx, nucleotide in enumerate(nucleotides):
        profile_dict[nucleotide] = matrix[idx][j]
    return profile_dict


def score(motifs):
    '''Returns the score of the dna list motifs.'''
    score = 0
    for i in range(len(motifs[0])):
        motif = ''.join([motifs[j][i] for j in range(len(motifs))])

        score += min([CalHammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
    return score
    


def GibbsSamplerWithSA(dna, k, t, N, initial_temperature, cooling_rate):
    randomIndex = [random.randint(0, len(dna[0])-k) for _ in range(len(dna))]
    Motifs = [(seq[i:i+k], i) for seq, i in zip(dna, randomIndex)]  # Store motif and its starting index

    BestMotifs = Motifs
    best_score = score([motif[0] for motif in BestMotifs])
    temperature = initial_temperature
    time_taken = []
    
    for _ in range(N):
        start_time = time.time()  # Start time for each iteration
        i = random.randint(0, len(dna) - 1)
        SampleMotifs = [dna[index] for index in range(len(dna)) if index != i]
        ProfileMatrix = GenerateProfileMatWithPesudoCount([motif[0] for motif in Motifs])
        NewMotif = FindPortableKmer(dna[i], k, ProfileMatrix)
        Motifs[i] = NewMotif

        new_score = score([motif[0] for motif in Motifs])
        delta_score = new_score - best_score
        
        if delta_score < 0 or random.random() < math.exp(-delta_score / temperature):
            BestMotifs = Motifs
            best_score = new_score

        time_taken.append(time.time() - start_time)

        # Decrease temperature
        temperature *= cooling_rate
    
    return BestMotifs, time_taken


res = []
k_motif = int(input("Enter k: ")) 
initial_temperature = 100
cooling_rate = 0.95
N_values = [500, 1000, 1500, 2000]

with open("hm03.txt", "r") as f:
     dnaSet = [line.strip() for line in f]

for N in N_values:
    run_results = []
    for _ in range(20):
        best_motifs, time_taken = GibbsSamplerWithSA(dnaSet, k_motif, len(dnaSet), N, initial_temperature, cooling_rate)
        run_results.append((best_motifs, time_taken))
    res.append((N, run_results))

# Save the output to a text file
with open("output_with_SA.txt", "w") as output_file:
    for N, run_results in res:
        output_file.write(f"k = {k_motif}, N = {N}\n")
        best_index = min(range(len(run_results)), key=lambda i: score([motif[0] for motif in run_results[i][0]]))
        output_file.write("Minimum Score: {}\n".format(score([motif[0] for motif in run_results[best_index][0]])))
        output_file.write("Index of the best motifs: {}\n".format(best_index))
        for motif, time_taken in zip(run_results[best_index][0], run_results[best_index][1]):
            output_file.write("Motif: {}, Start Index: {}, Time Taken: {}\n".format(motif[0], motif[1], time_taken))
        output_file.write("\n")
