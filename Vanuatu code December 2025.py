import random
import statistics
import winsound
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

class individual:
    
    def __init__(self, chunks_id = [],chunk_mutations = [], fitness = 1):
        self.chunk_ids = chunks_id
        self.chunk_mutations = chunk_mutations
        self.fitness = fitness

    def __str__(self):
        return 'ind with chunks: {self.chunk_ids}, mutations: {self.chunk_mutations}, fitness: {self.fitness}\n'.format(self=self)

    def from_scratch(n_chunks, pop_1_or_2 = 1):
        if pop_1_or_2 == 1:
            ids = 1
        elif pop_1_or_2 == 2:
            ids = -1
        ind = individual()
        ind.chunk_ids = individual.generate_chunks(n_chunks, ids)
        ind.chunk_mutations = np.zeros([n_chunks,2])
        ind.fitness = 1
        return ind

    def generate_chunks(n_chunks, ids):
        chunks_id = []
        for chunk in range(n_chunks):
            chunks_id += [[random.random()*ids, random.random()*ids]]
        return chunks_id

    def from_parents(parent_a, parent_b, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate):
        ind = individual(chunks_id = [],chunk_mutations = [], fitness = 1)
        n_chunks = len(parent_a.chunk_ids)
        a_chunk_choices = recombine(n_chunks, recomb_rate)
        b_chunk_choices = recombine(n_chunks, recomb_rate)
        mutations = np.random.poisson(del_mut_lambda/(2*n_chunks),2*n_chunks)
        mut_a = mutations[:n_chunks]
        mut_b = mutations[n_chunks:]

        for chunk in range(n_chunks):
            ind.chunk_ids.append([parent_a.chunk_ids[chunk][a_chunk_choices[chunk]],parent_b.chunk_ids[chunk][b_chunk_choices[chunk]]])
            ind.chunk_mutations.append([parent_a.chunk_mutations[chunk][a_chunk_choices[chunk]] + mut_a[chunk],
                                        parent_b.chunk_mutations[chunk][b_chunk_choices[chunk]] + mut_b[chunk]])
            ind.fitness *= del_mut_harm**sum(ind.chunk_mutations[-1])
            if ind.chunk_ids[-1][0] == ind.chunk_ids[-1][1]:
                ind.fitness *= homo_harm

        #print("example chunk:",[parent_a.chunk_ids[chunk][a_chunk_choices[chunk]]],"!")
        return ind
            

class population:
    
    def __init__(self, ind_lst):
        self.ind_lst = ind_lst
        self.chunk_ids = []
        self.chunk_mutations = []
        self.fitness = []
        for ind in ind_lst:
            self.chunk_ids.append([ind.chunk_ids])
            self.chunk_mutations.append(ind.chunk_mutations)
            self.fitness.append(ind.fitness)

    def __rpr__(self):
        return self.ind_lst

    def add(pop1,pop2):
        pop1.ind_lst += pop2.ind_lst
        pop1.chunk_ids += pop2.chunk_ids
        pop1.chunk_mutations += pop2.chunk_mutations
        pop1.fitness += pop2.fitness

def new_pop(pop_size,n_chunks,pop_1_or_2 = 1):
    pop = []
    for person in range(pop_size):
        ind = individual.from_scratch(n_chunks, pop_1_or_2)
        pop.append(ind)
    pop = population(pop)
    return pop

def new_homogenous_pop(pop_size,n_chunks, homo_harm,pop_1_or_2 = 1):
    pop = []
    genome = []
    if pop_1_or_2 == 1:
            ids = 1
    elif pop_1_or_2 == 2:
        ids = -1
    allele = random.random()*ids
    genome += [[allele, allele]]
    for person in range(pop_size):
        ind = individual(chunks_id = genome, chunk_mutations = np.zeros([n_chunks,2]), fitness = homo_harm**n_chunks)
        pop.append(ind)
    pop = population(pop)
    return pop

def recombine2(n_chunks, recomb_rate = 1):
    
    chosen_chunks = random.choices([0,1])
    for chunk in range(n_chunks-1):
        chosen_chunks += random.choices([chosen_chunks[-1],0**chosen_chunks[-1]], weights = [1,recomb_rate])

    return chosen_chunks

def recombine(n_chunks, recomb_rate):

    current_chrom = random.choices([0,1])[0]
    recombine_at = random.choices([0,1],[1,recomb_rate],k=n_chunks)
    chosen_chunks = [current_chrom]
    
    for chunk in range(n_chunks-1):
        if recombine_at[chunk] == 1:
            current_chrom = 0**current_chrom
        chosen_chunks += [current_chrom]

    return chosen_chunks

def next_gen2(pop, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity):
    #Mates the existing individuals, based on their fitness, and creates a new generation#
    
    if type(pop) != population:
        pop = population(pop)
    
    n_chunks = len(pop.ind_lst[0].chunk_ids)-1
        
    fitness_lst = pop.fitness
    
    pop_size = 0
    new_gen = []
    ind_ids = list(range(len(pop.ind_lst)))
    
    while pop_size < capacity:
        a_chunk_choices = recombine(n_chunks, recomb_rate)
        b_chunk_choices = recombine(n_chunks, recomb_rate)
        
        parent_a = random.choices(ind_ids, fitness_lst)[0]
        parent_b = random.choices(ind_ids, fitness_lst)[0]
        
        while parent_a == parent_b:            
            parent_b = random.choices(ind_ids, fitness_lst)[0]

        offspring = individual.from_parents(pop.ind_lst[parent_a], pop.ind_lst[parent_b], homo_harm, del_mut_lambda, del_mut_harm, recomb_rate)

        new_gen += [offspring]
        pop_size += 1

    new_gen = population(new_gen)

    return new_gen

def next_gen(pop, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity):
    #Mates the existing individuals, based on their fitness, and creates a new generation#
    
    if type(pop) != population:
        pop = population(pop)
    
    n_chunks = len(pop.ind_lst[0].chunk_ids)-1
        
    fitness_lst = pop.fitness
    
    pop_size = 0
    new_gen = []
    ind_ids = list(range(len(pop.ind_lst)))
    parents = random.choices(ind_ids, fitness_lst, k = 2*capacity)
    
    for i in range(capacity):
        
        parent_a = parents[2*i]
        parent_b = parents[2*i+1]
        
        while parent_a == parent_b:            
            parent_b = random.choices(ind_ids, fitness_lst)[0]

        offspring = individual.from_parents(pop.ind_lst[parent_a], pop.ind_lst[parent_b], homo_harm, del_mut_lambda, del_mut_harm, recomb_rate)

        new_gen += [offspring]
        pop_size += 1

    new_gen = population(new_gen)

    return new_gen


def allele_freq(pop):
    #Calculates the frequencies of alleles originating in each population#

    pop1_alleles = 0
    pop2_alleles = 0 #the total number of alleles originating in population 2
    
    for ind in pop.chunk_ids:
        for chunk in ind[0]:
            for allele in chunk:
                if allele > 0:
                    pop1_alleles += 1
                else:
                    pop2_alleles += 1

    return [pop1_alleles/(pop1_alleles+pop2_alleles), pop2_alleles/(pop1_alleles+pop2_alleles)]

def detailed_allele_freq(pop):

    total_alleles = len(pop.chunk_ids)*len(pop.chunk_ids[0][0])*2 #the total number of alleles in the population
    d = {}
    n_alleles = 0

    for ind in pop.chunk_ids:
        for chunk in ind[0]:
            for allele in chunk:
                try:
                    d[allele] += 1/total_alleles
                except:
                    d[allele] = 1/total_alleles
                    n_alleles += 1

    print(n_alleles)

    return d
    

def graph(freq_matrix, legend = False):
    #Creates a graph of allele frequencies from both sources per generation, starting from the arrival of the second wave
    
    matrix = np.array(freq_matrix)
    from_1 = matrix[:,0]
    from_2 = matrix[:,1]
    x=range(len(freq_matrix))
    plt.stackplot(x,from_1, from_2, labels=['From population 1','From population 2'])
    if legend == True:
        plt.legend(loc='upper left')

    return plt.show()

def pop_test(n_gens):
    
    pop = new_pop(2500,5,2)
    pop = population(pop)
    gen = 1
    
    while gen < n_gens:
        gen += 1
        pop = next_gen(pop,1, 1, 0.999, 1, 10000)
    print ("\n", np.mean(pop.fitness))
    
    return None

def drift_test():
    for i in range(1000):
        pop = new_pop(3,1,2)
        pop = next_gen(pop,1,1,1,1,3)
        detailed_allele_freq(pop)
    return None

def ind_test(n_gens):

    gen = 1
    print ("GEN", gen)
    ind1 = individual.from_scratch(1,2)
    ind2 = individual.from_scratch(1,2)

    print("ind 1:", ind1.chunk_ids, "\nind 2:", ind2.chunk_ids)

    while gen < n_gens:
        gen += 1
        print ("\nGEN", gen)
        ind3 = individual.from_parents(ind1,ind2, 0.9, 1, 0.95, 0.5)
        ind4 = individual.from_parents(ind1,ind2, 0.9, 1, 0.95, 0.5)
        ind1 = ind3
        ind2 = ind4
        print("\nind 1:", ind1.chunk_ids, "\nind 2:", ind2.chunk_ids)

    return None

def exp(t1 = 60, t2 = 40, t3 = 100, n_chunks = 100, wave_1_size = 25, capacity = 25000, homo_harm = 0.99, waves_2_size = 20,
         waves_interval = None, del_mut_lambda = 2.2, del_mut_harm = 0.99, recomb_rate = 0.01, sd = None, min_wave_size = 2):
    #Simulates the entire process#

    now = datetime.now()

    current_time = now.strftime("%H:%M:%S")
    print("Start Time =", current_time)
    if del_mut_harm == 1:
        del_mut_lambda = 0
        print("Deleterious mutations are not modelled to save time, since they don't affect fitness anyway")
    
    #Creating the waves of migration#
    wave1 = new_pop(pop_size = wave_1_size ,n_chunks = n_chunks ,pop_1_or_2 = 1)

    #Running the generations of wave 1's journy to the new deme#
    gen1 = 0
    if sd != None:
        wave_sizes = []
        while len(wave_sizes) < t1:
            size = round(np.random.normal(wave_1_size,wave_1_size*sd))
            if size >= min_wave_size:
                wave_sizes += [size]
            else:
                wave_sizes += [min_wave_size] #See if we want that later
        print(wave_sizes)
        while gen1 < t1:
            wave1 = next_gen(wave1, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity = wave_sizes[gen1])
            gen1 += 1
        print("Wave 1 arrived at the deme!")
        print("Average fitness =", np.mean(wave1.fitness))
    else:
        while gen1 < t1:
            wave1 = next_gen(wave1, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity = wave_1_size)
            gen1 += 1
        print("Wave 1 arrived at the deme!")
        print("Average fitness =", np.mean(wave1.fitness))
    
    #Running the generations spent alone in the deme#
    gen2 = 0
    deme_pop = wave1
    while gen2 < t2:
        deme_pop = next_gen(deme_pop, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity)
        gen2 += 1
    print("\n", t2, "generatuins spent alone at the deme!")
    print("Average fitness =", np.mean(deme_pop.fitness))

    #Migration waves of population 2 start arriving
    gen3 = 0
    gen_since_last_wave = waves_interval
    freq_matrix = []

    while gen3 < t3:
        if gen_since_last_wave == waves_interval:
            deme_pop.add(new_pop(pop_size = waves_2_size ,n_chunks = n_chunks ,pop_1_or_2 = 2))
            gen_since_last_wave = 0
        deme_pop = next_gen(deme_pop, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity)
        freq_matrix += [allele_freq(deme_pop)]
        gen3 += 1
        gen_since_last_wave += 1
        
    print("\n", t3, "generatuins spent drizzling!")
    print("Average fitness =", np.mean(deme_pop.fitness))

    #Calculating allele frequency#
    allele_frequency = allele_freq(deme_pop)

    #finishing sound#
    winsound.Beep(400,200)
    winsound.Beep(400,200)
    winsound.Beep(600,400)

    print ("\nDONE!")
    print("\nAllele frequencies are:\n", allele_freq(deme_pop)[0], "from the first population\n", allele_freq(deme_pop)[1], "from the second population")
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Finish Time =", current_time)
    
    return freq_matrix

def panel(t1, t2, t3, n_chunks, wave_1_size, capacity, homo_harm, waves_2_size,
         waves_interval, del_mut_lambda, del_mut_harm, recomb_rate):
    #Simulates the entire process#

    now = datetime.now()

    current_time = now.strftime("%H:%M:%S")
    print("Start Time =", current_time)
    
    #Creating the waves of migration#
    wave1 = new_pop(pop_size = wave_1_size ,n_chunks = n_chunks ,pop_1_or_2 = 1)

    #Running the generations of wave 1's journy to the new deme#
    gen1 = 0
    while gen1 < t1:
        wave1 = next_gen(wave1, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity = wave_1_size)
        gen1 += 1
    print("Wave 1 arrived at the deme!")
    print("Average fitness =", np.mean(wave1.fitness))
    
    #Running the generations spent alone in the deme#
    gen2 = 0
    deme_pop = wave1
    while gen2 < t2:
        deme_pop = next_gen(deme_pop, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity)
        gen2 += 1
    print("\n", t2, "generatuins spent alone at the deme!")
    print("Average fitness =", np.mean(deme_pop.fitness))

    #Migration waves of population 2 start arriving
    gen3 = 0
    gen_since_last_wave = waves_interval
    freq_matrix = []

    while gen3 < t3:
        if gen_since_last_wave == waves_interval:
            deme_pop.add(new_pop(pop_size = wave_1_size ,n_chunks = n_chunks ,pop_1_or_2 = 2))
            gen_since_last_wave = 0
        deme_pop = next_gen(deme_pop, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity)
        freq_matrix += [allele_freq(deme_pop)]
        gen3 += 1
        gen_since_last_wave += 1
        
    print("\n", t3, "generatuins spent drizzling!")
    print("Average fitness =", np.mean(deme_pop.fitness))

    #Calculating allele frequency#
    allele_frequency = allele_freq(deme_pop)

    #finishing sound#
    winsound.Beep(400,200)
    winsound.Beep(400,200)
    winsound.Beep(600,400)

    print ("\nDONE!")
    print("\nAllele frequencies are:\n", allele_freq(deme_pop)[0], "from the first population\n", allele_freq(deme_pop)[1], "from the second population")
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Finish Time =", current_time)
    
    return graph(freq_matrix)

def exp_homogenous(t1 = 60, t2 = 40, t3 = 100, n_chunks = 1, wave_1_size = 25, capacity = 2500, homo_harm = 0.9, waves_2_size = 20,
         waves_interval = None, recomb_rate = 0.5):

    now = datetime.now()
    del_mut_lambda = 0
    del_mut_harm = 1

    current_time = now.strftime("%H:%M:%S")
    print("Start Time =", current_time)
    
    #Creating the waves of migration#
    deme_pop = new_homogenous_pop(pop_size = capacity ,n_chunks = n_chunks,  homo_harm =  homo_harm ,pop_1_or_2 = 1)

    print("Average fitness =", np.mean(deme_pop.fitness))

    #Migration waves of population 2 start arriving
    gen3 = 0
    gen_since_last_wave = waves_interval
    freq_matrix = []

    while gen3 < t3:
        if gen_since_last_wave == waves_interval:
            deme_pop.add(new_homogenous_pop(pop_size = wave_1_size ,n_chunks = n_chunks, homo_harm = homo_harm ,pop_1_or_2 = 2))
            gen_since_last_wave = 0
        deme_pop = next_gen(deme_pop, homo_harm, del_mut_lambda, del_mut_harm, recomb_rate, capacity)
        freq_matrix += [allele_freq(deme_pop)]
        gen3 += 1
        gen_since_last_wave += 1
        
    print("\n", t3, "generatuins spent drizzling!")
    print("Average fitness =", np.mean(deme_pop.fitness))

    #Calculating allele frequency#
    allele_frequency = allele_freq(deme_pop)

    #finishing sound#
    winsound.Beep(400,200)
    winsound.Beep(400,200)
    winsound.Beep(600,400)

    print ("\nDONE!")
    print("\nAllele frequencies are:\n", allele_freq(deme_pop)[0], "from the first population\n", allele_freq(deme_pop)[1], "from the second population")
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Finish Time =", current_time)
    
    return (freq_matrix)

def plot_single(t1 = 60, t2 = 40, t3 = 100, n_chunks = 100, wave_1_size = 25, capacity = 25000, homo_harm = 0.975, waves_2_size = 25,
         waves_interval = None, del_mut_lambda = 2.2, del_mut_harm = 1, recomb_rate = 0.01, sd = None, min_wave_size = 2):
    freq_matrix = exp(t1, t2, t3, n_chunks, wave_1_size, capacity, homo_harm, waves_2_size, waves_interval, del_mut_lambda, del_mut_harm, recomb_rate, sd, min_wave_size)
    matrix = np.array(freq_matrix)
    from_1 = matrix[:,0]
    from_2 = matrix[:,1]
    x=range(len(freq_matrix))
    plt.stackplot(x,from_1, from_2, labels=['Austronesian','Papuan'])
    plt.legend(loc='upper right')
    
    plt.show()
    return None

    

def figure_1(t1 = 60, t2 = 40, t3 = 100, n_chunks = 100, wave_1_size = 25, capacity = 25000, homo_harms = [1,0.99,0.95], waves_2_size = 25,
         waves_interval = None, del_mut_lambdas = [0,1,2.2], del_mut_harm = 0.99, recomb_rate = 0.1):
    fig, axs = plt.subplots(len(homo_harms), len(del_mut_lambdas), sharex = True)
    fig.supylabel("Inbreeding Depression")
    fig.suptitle("Deleterious Mutations λ")
    axs[0,0].set_title("0", size = 10)
    axs[0,1].set_title("1", size = 10)
    axs[0,2].set_title("2.2", size = 10)
    axs[0,0].set_ylabel("0\n\nAncestry composition")
    axs[1,0].set_ylabel("1%\n\n")
    axs[2,0].set_ylabel("5%\n\n")
    axs[2,0].set_xlabel("Generation")
    for i in range(len(homo_harms)):
        for j in range(len(del_mut_lambdas)):
            print(i,j)
            
            homo_harm = homo_harms[i]
            del_mut_lambda = del_mut_lambdas[j]

            freq_matrix = exp(t1, t2, t3, n_chunks, wave_1_size, capacity, homo_harm, waves_2_size, waves_interval, del_mut_lambda, del_mut_harm, recomb_rate)
            matrix = np.array(freq_matrix)
            from_1 = matrix[:,0]
            from_2 = matrix[:,1]
            x=range(len(freq_matrix))
            axs[i,j].stackplot(x,from_1, from_2, labels=['Austronesian','Papuan'])
            if i == 0 and j == len(del_mut_lambdas)-1:
                axs[i,j].legend(loc='upper right')
            if j > 0:
                axs[i,j].set_yticks([])
             
    plt.show()
    return None

def figure_2(t1 = 60, t2 = 40, t3 = 100, n_chunks = 100, wave_1_sizes = [25,100,250,1000], sds = [0.25,0.5,0.75], capacity = 25000, homo_harm = 0.95, waves_2_size = 25,
         waves_interval = None, del_mut_lambda = 2.2, del_mut_harm = 0.99, recomb_rate = 0.1, min_wave_size = 2):
    fig, axs = plt.subplots(len(sds), len(wave_1_sizes), sharex = True)
    fig.supylabel("Standard deviation")
    fig.suptitle("Mean wave size")
    axs[0,0].set_title(str(wave_1_sizes[0]), size = 10)
    axs[0,1].set_title(str(wave_1_sizes[1]), size = 10)
    axs[0,2].set_title(str(wave_1_sizes[2]), size = 10)
    axs[0,3].set_title(str(wave_1_sizes[3]), size = 10)
    axs[0,0].set_ylabel(str(sds[0])+"\n\nAncestry composition")
    axs[1,0].set_ylabel(str(sds[1])+"\n\n")
    axs[2,0].set_ylabel(str(sds[2])+"\n\n")
    axs[2,0].set_xlabel("Generation")
    for i in range(len(sds)):
        for j in range(len(wave_1_sizes)):
            print(i,j)
            
            sd = sds[i]
            wave_1_size = wave_1_sizes[j]

            freq_matrix = exp(t1, t2, t3, n_chunks, wave_1_size, capacity, homo_harm, waves_2_size, waves_interval, del_mut_lambda, del_mut_harm, recomb_rate, sd, min_wave_size)
            matrix = np.array(freq_matrix)
            from_1 = matrix[:,0]
            from_2 = matrix[:,1]
            x=range(len(freq_matrix))
            axs[i,j].stackplot(x,from_1, from_2, labels=['Austronesian','Papuan'])
            if i == 0 and j == len(wave_1_sizes)-1:
                axs[i,j].legend(loc='upper right')
            if j > 0:
                axs[i,j].set_yticks([])
             
    plt.show()
    return None
