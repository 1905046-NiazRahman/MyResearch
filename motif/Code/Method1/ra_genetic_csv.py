import os
import numpy as np
from random import randint, choice, random
import time
from scipy.stats import binom
import csv
import numpy as np

def calculate_motif_p_value(motif, dna):
    prob_of_motif = (1/4) ** len(motif) 
    occurrences = sum(seq.count(motif) for seq in dna)
    total_positions = sum(len(seq) - len(motif) + 1 for seq in dna)
    p_value = 1 - binom.cdf(occurrences - 1, total_positions, prob_of_motif)
    return p_value

def fitness(motifs, dna):
    '''Calculate fitness using the negative logarithm of p-values.'''
    total_fitness = 0
    for motif in motifs:
        p_value = calculate_motif_p_value(motif, dna)
        if p_value > 0:
            # Apply negative logarithm to transform the p-value
            # Adding a small constant to p_value to avoid log(0)
            total_fitness += -np.log(p_value + 1e-100)
        else:
            # Handle the case where p_value is 0 by assigning a large fitness value
            total_fitness += 1000  # Arbitrary large value for motifs with p_value of 0
    return total_fitness


def generate_initial_population(dna, k, population_size):
    population = []
    positions = []  # Track start positions
    for _ in range(population_size):
        start_positions = [randint(0, len(seq) - k) for seq in dna]
        motifs = [seq[pos:pos+k] for seq, pos in zip(dna, start_positions)]
        population.append(motifs)
        positions.append(start_positions) 
    return population, positions  

def selection_with_positions(population, positions, dna, top_n=2):
    fitness_scores = [(motifs, pos, fitness(motifs, dna)) for motifs, pos in zip(population, positions)]
    sorted_population = sorted(fitness_scores, key=lambda x: x[2], reverse=True)
    selected_motifs = [motifs for motifs, _, _ in sorted_population[:top_n]]
    selected_positions = [pos for _, pos, _ in sorted_population[:top_n]]
    return selected_motifs, selected_positions

def crossover_with_positions(motif_set1, motif_set2, positions1, positions2):
    crossover_point = randint(1, len(motif_set1) - 1)
    new_motif_set1 = motif_set1[:crossover_point] + motif_set2[crossover_point:]
    new_motif_set2 = motif_set2[:crossover_point] + motif_set1[crossover_point:]
    return new_motif_set1, new_motif_set2, positions1, positions2

def mutate_with_positions(motifs, positions, mutation_rate=0.1):
    mutated_motifs = []
    for motif in motifs:
        if random() < mutation_rate:
            mutation_pos = randint(0, len(motif) - 1)
            new_nucleotide = choice(['A', 'C', 'G', 'T'])
            motif = motif[:mutation_pos] + new_nucleotide + motif[mutation_pos+1:]
        mutated_motifs.append(motif)
    return mutated_motifs, positions  

def score(motifs):
    score = 0
    for i in range(len(motifs[0])):
        column = [motif[i] for motif in motifs]
        max_freq = max(column.count(nucleotide) for nucleotide in 'ACGT')
        score += len(motifs) - max_freq
    return score*.4

def genetic_algorithm(dna, k, population_size, generations, score_threshold):
    population, positions = generate_initial_population(dna, k, population_size)  
    best_score = 1e9
    best_motifs, best_positions = None, None

    for _ in range(generations):
        selected, selected_positions = selection_with_positions(population, positions, dna)  
        new_population, new_positions = [], []
        for i in range(len(selected) // 2):
            offspring1, offspring2, pos1, pos2 = crossover_with_positions(selected[i], selected[-i-1], selected_positions[i], selected_positions[-i-1]) 
            offspring1, pos1 = mutate_with_positions(offspring1, pos1)  
            offspring2, pos2 = mutate_with_positions(offspring2, pos2)
            new_population.extend([offspring1, offspring2])
            new_positions.extend([pos1, pos2])
        population, positions = new_population, new_positions

        current_best_motifs, current_best_positions = max(zip(population, positions), key=lambda x: fitness(x[0], dna))
        current_best_score = score(current_best_motifs)
        if current_best_score < best_score:
            best_score = current_best_score
            best_motifs, best_positions = current_best_motifs, current_best_positions

    return best_motifs, best_positions, best_score 

def create_directory(path):
    """Create a directory if it does not exist."""
    if not os.path.exists(path):
        os.makedirs(path)

def run_genetic_algorithm_to_csv(dna, population_sizes, mutation_rates, k, generations, csv_path):
    itr_cnt = 1
    with open(csv_path, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Population Size', 'Mutation Rate', 'Best Score', 'Execution Time', 'Best Motifs', 'Start Positions'])
        
        for size in population_sizes:
            for rate in mutation_rates:
                
                tot_time = 0
                best_motifs = []
                best_positions = []
                best_score = int(1e9)

                for _ in range(itr_cnt):
                    start_time = time.time()
                    motifs, positions, score = genetic_algorithm(dna, k, size, generations, rate)
                    end_time = time.time()
                    tot_time += end_time - start_time
                    if score < best_score:
                        best_motifs = motifs
                        best_positions = positions
                        best_score = score
                
                # Write results to the CSV file
                positions_str = "; ".join(map(str, best_positions))
                writer.writerow([size, rate, int(best_score), tot_time / itr_cnt, "; ".join(best_motifs), positions_str])

population_sizes = range(10, 101, 10)  # 10 to 100 in steps of 10
mutation_rates = np.arange(0.1, 1.1, 0.1)  # 0.1 to 1.0 in steps of 0.1
generations = 100  # Number of generations

input_files = ['hm03.txt', 'yst04r.txt', 'yst08r.txt']
ks = [10, 15, 20]

for input_file in input_files:
    with open(f'inputs/{input_file}', 'r') as file:
        dnaSet = file.read().splitlines()
    
    for k in ks:
        output_csv = f'outputs/genetic/{input_file.split(".")[0]}_motif_len_{k}.csv'
        # Ensure the output directory exists
        os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        
        # Execute the algorithm with varying parameters and save outputs to a CSV
        run_genetic_algorithm_to_csv(dnaSet, population_sizes, mutation_rates, k, generations, output_csv)
