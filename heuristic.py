
from xml.dom import minidom
import time
from collections import defaultdict
import random

class Graph:

    def __init__(self, vertices):
        self.V = vertices
        self.graph = defaultdict(list)

    def read_answer(self, k_list):
        sequence = ''
        for i in k_list:
            sequence += i[0]
        sequence += i[1:]
        return sequence


    def print_all_paths_util(self, u, visited, path, all):
        global N, S0
        
        visited[u] -= 1
        path.append(u)

        # wybierz wierzchołek
        to_choose_from = []
        for j in self.graph[u]:
            if visited[j] > 0:
                to_choose_from.append(j)
        
        if len(to_choose_from) > 0: 
            # wybierz losowy wierzchołek
            chosen_j = random.choice(to_choose_from)
            self.print_all_paths_util(chosen_j, visited, path, all)
        else :
            # jeżeli nie ma następnika, to dodaj ścieżkę
            answer = self.read_answer(path)
            if answer not in all:
                all.append(answer)
        
        visited[u] += 1

    
    def print_all_paths(self, path, all):
        global S0, visited_with_counter
        self.print_all_paths_util(S0, visited_with_counter, path, all)
        
    def create_graph(self, olis):
        global K
        for oli in olis:
            self.graph[oli] = []
        
        # dodaj każdy oli gdzie wszedzie jest następnikiem
        for oli in olis:
            for key in self.graph.keys():
                if oli[0:(K-1)] == key[1:K] and oli != key:
                    # dodaj krawędź
                    self.graph[key].append(oli)

    def print_graph(self):
        for key, value in self.graph.items():
            print("-->", key, ":  ", value, end="\n")
        

def greedy_algorithm():
    global N, olis
    path = []
    all = []
    g = Graph(N)
    g.create_graph(olis)
    g.print_all_paths(path, all) # create answer
    return all[0]

# --------------- HEURISTIC -----------------

# increase lmer count (after inserting a solution or deleting it)
def count_lmers(olis_set):
    global lmers
    for oli in olis_set:
        if oli in lmers:
            lmers[oli] += 1
        else:
            lmers[oli] = 1


def extend_move_insert(mer):
    global best_solution, tabu, K
    
    # check if can be improved
    if best_solution[(len(best_solution)-K):] != mer[:-1]:
        # return false if no improvement
        return False
    
    # insert it in solution
    best_solution += mer[-1]
    tabu.append(mer)
    return True
    

def restart(reference_set, greedy_solution):
    global lmers

    # swap solutions
    min_length_solution = float('inf')
    to_swap = -1
    for i, solution in enumerate(reference_set):
        if len(solution) < min_length_solution:
            to_swap = i
            min_length_solution = len(solution)

    # delete olis from lmers from to swap
    for i in range(len(reference_set[to_swap])-1):
        to_delete = reference_set[to_swap][i:(i+K)]
        if to_delete in lmers:
            lmers[to_delete] -= 1

    # add olis to lmers from greedy_solution
    for i in range(len(greedy_solution)-1):
        to_add = greedy_solution[i:(i+K)]
        if to_add in lmers:
            lmers[to_add] += 1

    reference_set[to_swap] = greedy_solution
    return reference_set


def main_tabu():
    global best_solution, olis, K
    reference_set = [best_solution, '']

    # 3) while not all restarts done
    i = 0

    greedy_solution = greedy_algorithm()
    if len(reference_set) > 1:
            reference_set = restart(reference_set, greedy_solution)
    while i < 5:
        # replace the worst solution from reference_set by greedy solution
        #greedy_solution = greedy_algorithm()
        #if len(reference_set) > 1:
           # reference_set = restart(reference_set, greedy_solution)
        
        reference_set = sorted(reference_set, key=lambda x: len(x), reverse=True)
        best_solution = reference_set[0]
        print("reference_set after sorting ", reference_set)

        # 4)
        # 12, 13)
        is_inserted = False
        olis = sorted(lmers, key=lambda lmer: lmers[lmer])
        print("\nbest_solution before extending/deleting: ", best_solution)
        for oli in olis:
            is_inserted = extend_move_insert(oli)
            print("I checked addition for ", oli, " and it is ", is_inserted)
            if is_inserted:
                lmers[oli] += 1
                print("counter + ", oli)
                print("best_solution after adding: ", best_solution)
                break

        if not is_inserted:
            fragment = best_solution[len(best_solution)-K:]
            
            if fragment in lmers:
                if lmers[fragment] > 0:
                    lmers[fragment] -= 1
                    print("counter - ", fragment)
                    best_solution = best_solution[:-1]
                    print("best_solution after delete: ", best_solution)

        print("LMERS: ", lmers)
        reference_set[0] = best_solution
        print("reference_set[0] ", reference_set[0])
        olis = sorted(lmers, key=lambda lmer: lmers[lmer])
        i += 1
       
# --------------------------- UTIL ----------------------------------
        
def read_instance():
    global file, dna, N, S0, probe, K, olis, visited_with_counter, lmers
    # parsuj xml
    file = minidom.parse(input())
    dna = file.firstChild
    N = int(dna.getAttribute('length')) #długość badanej sekwencji
    S0 = dna.getAttribute('start') # początkowy fragment długości k
    probe = dna.getElementsByTagName('probe')[0]
    K = len(probe.getAttribute('pattern')) # długość sond oligonukleotydowych
   
    for oli in probe.getElementsByTagName('cell'):
        olis.append(oli.firstChild.nodeValue)

        # add to lmers set
        lmers[oli.firstChild.nodeValue] = 0

        for _ in range(int(oli.getAttribute('intensity'))):
            if oli.firstChild.nodeValue not in visited_with_counter:
                visited_with_counter[oli.firstChild.nodeValue] = int(oli.getAttribute('intensity'))
        

# ------------------------ MAIN --------------------------------

if __name__ == '__main__':
    
    file = -1
    dna = -1
    N = -1
    S0 = ''
    probe = ''
    K = -1
    olis = []
    visited_with_counter = {}
    lmers = {}

    read_instance()
    best_solution = greedy_algorithm()
    #print("greedy: ", best_solution)

    # HEURISTIC
    tabu = []
    main_tabu()
    # moves - add / delete / shift of nucleotide
    # tabu tab - list of inserted / deleted / shifted nucleotides remebered for certain amount of iterations
    
    # global evaluation -> length of sequence (for making the sequence longer)

    # (when algorithm is stuck) 
    # count number of times oligonucleotide is included in solutions, 
    # add l-mer with lowest value to solution 
    # else if no feasable insertion then deletion is applied
    # e.g. frequency value of l-mer = number of iterations
