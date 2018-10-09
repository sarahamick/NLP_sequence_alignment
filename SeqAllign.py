import re
    #PARAMETERS:
        #VOCAB_MAP: maps each letter in the protein vocabulary to its column/row
        #index in the cost_matrix f.x. if the protein at position i is A
        #in one string, and R in the other, you can find the cost of that
        #difference at position cost_matrix[As index in vocab_map][Rs index in vocab_map]
        #
        #COST_MATRIX: a 2D array that maps every letter in the vocabularly with every
        #other letter in the vocabulary. fx. in this cost matrix:   A  B
        #the value of a mismatch between A and B is -1.           A 2 -1
        #                                                         B -1 2
        #
        #ANIMAL_PROTEIN_MAP: a map of an animal name to its protein string


def allignment_loop(vocab_map, cost_matrix, sequence_1, sequence_2):

    sequence_1_aligned = ""
    sequence_2_aligned = ""

    A = [[0 for x in range(len(sequence_2)+1)] for y in range(len(sequence_1)+1)]

    delt = cost_matrix[vocab_map[sequence_1[1]]][vocab_map["*"]]
    for i in range(len(sequence_1)+1):
        A[i][0] = i * delt
    for j in range(len(sequence_2)+1):
        A[0][j] = j * delt

    for i in range(len(sequence_1)):
        for j in range(len(sequence_2)):

            alpha = cost_matrix[vocab_map[sequence_1[i]]][vocab_map[sequence_2[j]]] #cost of difference between letters
            delta1= cost_matrix[vocab_map["*"]][vocab_map[sequence_2[j]]] #cost of a gap in sequence1
            delta2 = cost_matrix[vocab_map[sequence_1[i]]][vocab_map["*"]] #cost of a gap in sequence2

            cost_diff = alpha + A[i][j]
            cost_gap_seq1 = delta2 + A[i][j+1]
            cost_gap_seq2 = delta1 + A[i+1][j]

            minimum = min(cost_diff, cost_gap_seq1, cost_gap_seq2)

            A[i+1][j+1] = minimum

    print("A[][]:")
    for line in A:
        print(line)
    print("\nCost:", A[len(sequence_1)][len(sequence_2)])


def allignment(vocab_map, cost_matrix, animal_protein_map):

    animal1_count = 0
    animal2_count = 0
    for animal1 in animal_protein_map:
        for animal2 in animal_protein_map:
            animal2_count = animal2_count + 1
            if animal2_count > animal1_count and not animal1 == animal2:
                print("\nCalculating sequence alignment for:")
                print(animal1, "(", animal_protein_map[animal1], ")")
                print(animal2, "(", animal_protein_map[animal2], ")\n")
                allignment_loop(vocab_map, cost_matrix, animal_protein_map[animal1], animal_protein_map[animal2])
        animal2_count = 0
        animal1_count = animal1_count + 1


def read_in_input():

    animal_to_protein = dict()

    regex = re.compile(r">")
    regex_linebreak = re.compile(r"\n")
    lines = []

    with open("data/Toy_FASTAs-in.txt") as f:
        lines = f.readlines()

    animal = ""
    protein = ""
    first_loop = True
    for line in lines:
        if regex.search(line) is not None:
            if not first_loop:
                animal_to_protein[animal] = protein
                animal = ""
                protein = ""
            animal = animal + line[1:].split(" ")[0]
            first_loop = False
        else:
            if regex_linebreak.search(line) is not None:
                protein = protein + line[:-1]
            else:
                protein = protein + line

    animal_to_protein[animal] = protein

    print("\nANIMAL->PROTEIN MAPPINGS:")
    for line in animal_to_protein:
        print(line, ">>>", animal_to_protein[line])
    return animal_to_protein

def generate_cost_matrix():
    regex = re.compile(r"\n")
    cost_matrix = []

    with open("data/BLOSUM62_switchedsigns.txt") as f:
        skip_first_line = f.readline()
        for line in f:
            cost_matrix.append(line.split(" "))

    for line in cost_matrix:

        for i in range(len(line)):
            if regex.search(line[i]) is not None:
                aux = line[i][:-1]
                line[i] = int(aux)
            else:
                aux = int(line[i])
                line[i] = aux

    print("COST MATRIX: ")
    for line in cost_matrix:
        print(line)
    return cost_matrix

def generate_vocab_mapping():
    regex = re.compile(r"\n")
    vocabulary = dict()

    with open("data/BLOSUM62.txt") as f:
        vocab = f.readline().split(" ")
        for i in range(len(vocab)):
            if regex.search(vocab[i]) is not None:
                vocabulary[vocab[i][:-1]] = i
            else:
                vocabulary[vocab[i]] = i

        print("\nPROTEIN LETTERS->COST MATRIX INDEX MAPPING")
        print(vocabulary, '\n')
    return vocabulary

def run():
    vocab_map = generate_vocab_mapping()
    cost_matrix = generate_cost_matrix()
    animal_protein_map = read_in_input()
    allignment(vocab_map, cost_matrix, animal_protein_map)

run()
