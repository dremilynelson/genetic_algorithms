from enum import Enum
from typing import Dict, List
import random
import math
import copy
import time
import os

GeneInstruction = Enum('Instruction', ['UP', 'DOWN', 'RIGHT', 'LEFT'])
MutationType = Enum('Mutation', ['AddGene', 'DeleteGene', 'MutateGene', 'DuplicateGene', 'Shuffle'])
ScoringTypes = Enum("ScoringType", ["TotalDistance", "DistanceFromTarget", 
                                    "DistanceMinusGenomeLength", "TotalDistanceWithGenomePenalty",
                                    "DistanceToTargetPlusDistanceBonus"])

### CONFIG ####
GRID_DIM = 40
START_CHAR = "O"
TARGET_CHAR = "X"
PASSED_CHAR = "+"
CRAWLER_CHAR = "@"
EMPTY_CHAR = "-"
MAX_DISTANCE = 5
STARTING_GENOME_SIZE = 7
MAX_GENOME_SIZE = 25
GENERATION_SIZE = 25
TOTAL_GENERATIONS = 50
SLEEP_TIME = 1
GENOME_LENGTH_PENALTY = 3
GENOME_LENGTH_BONUS = 3
DISTANCE_TRAVELED_BONUS = 3

class Genome:
    # holds a list of genes
    genes: List 

    def __init__(self, genes = None):
        # if not provided with a gene list, randomizes a genome
        if genes is None:
            self.genes = []
            for i in range(0, STARTING_GENOME_SIZE):
                self.genes.append(Gene())
        # otherwise, uses provided list
        else:
            self.genes = genes

    def __str__(self):
        s = "Genome:\n"
        for i in range(0, len(self.genes)):
            s += self.genes[i].__str__() + "\n"

        return(s)
    
    def mutate(self):
        # handles mutation event

        new_gene_list = copy.deepcopy(self.genes)
        # pick mutation class 
        mutation_type = random.choice(list(MutationType))
        if mutation_type == MutationType.AddGene:
            if len(new_gene_list) >= MAX_GENOME_SIZE:
                pass # nothing happens, genome is max size
            else:
                added_gene = Gene() # random gene
                new_gene_list.append(added_gene)
        elif mutation_type == MutationType.DeleteGene:
            if len(self.genes) > 1:
                gene_to_delete = random.choice(new_gene_list)
                new_gene_list.remove(gene_to_delete)
        elif mutation_type == MutationType.MutateGene:
            random_gene = random.choice(new_gene_list)
            random_gene.mutate()
        elif mutation_type == MutationType.DuplicateGene:
            if len(new_gene_list) >= MAX_GENOME_SIZE:
                pass
            else:
                random_gene = random.choice(new_gene_list)
                new_gene_list.append(Gene(random_gene.instruction, random_gene.distance))
        elif mutation_type == MutationType.Shuffle:
            random.shuffle(new_gene_list)

        return(Genome(genes=new_gene_list))

class Position:
    # class to contain info about a point on a grid
    x: int
    y: int 
    has_been_passed: bool

    def __init__(self, x=0, y=0):
        self.x=x
        self.y=y
        self.has_been_passed = False

    def __str__(self):
        return(f"({self.x}, {self.y})")

    def equals(self, other_position):
        if self.x == other_position.x and self.y == other_position.y:
            return True 
        else: 
            return False

    # given a chance in x and y, calculate coordinates of new position 
    def calc_new_coords(self, x=0, y=0):
        new_x = self.x + x
        new_y = self.y + y
        if new_x >= GRID_DIM or new_x < 0:
            # torus
            new_x = (new_x) % (GRID_DIM)

        if new_y >= GRID_DIM or new_y < 0:
            new_y = (new_y) % (GRID_DIM)
    
        return(new_x, new_y)
    
    # calculate distance between 2 positions
    def distance_from(self, other_position):
        return math.dist(
            (self.x, self.y),
            (other_position.x, other_position.y)
        )
    
class Screen:
    # class to hold Positions in a grid
    # grid is toroidal
    target: Position
    start: Position
    grid: Dict

    def __init__(self, target, start=Position(0, 0)):
        self.target = target
        self.start = start

        self.reset()

    def reset(self):
        grid = {}

        for i in range(0, GRID_DIM):
            grid[i] = {}
            for j in range(0, GRID_DIM):
                pos = Position(i, j)
                if pos.equals(self.target):
                    grid[i][j] = self.target 
                elif pos.equals(self.start):
                    grid[i][j] = self.start
                else:
                    grid[i][j] = pos
        self.grid = grid


    # use configured symbols to draw ASCII picture of path
    def draw_with_crawler_path(self, crawler = None):
        for i in range(0, GRID_DIM):
            row_string = ""
            for j in range(0, GRID_DIM):
                pos = self.grid[i][j]
                if pos.equals(self.target):
                    row_string += TARGET_CHAR + " "
                elif pos.equals(self.start):
                    row_string += START_CHAR + " "
                elif crawler is not None and crawler.position.equals(pos):
                    row_string += CRAWLER_CHAR + " "
                elif pos.has_been_passed:
                    row_string += PASSED_CHAR + " "
                else:
                    row_string += EMPTY_CHAR + " "
            
            print(row_string)

class Crawler:
    # an organism with a genome. moves over a grid.
    genome: Genome
    position: Position
    grid: Screen
    distance_traveled: int
    score: int

    def __init__(self, genome, parent = None):
        self.genome = genome 
        self.grid = Screen(target = TARGET_POSITION, start = START_POSITION)
        self.position = self.grid.start
        self.grid = self.grid
        self.distance_traveled = 0
        self.score = -1
        self.parent = parent

    def __str__(self):
        s = f"Crawler: score={self.score}"
        return(s)

    def act(self):
        # activate each gene in order
        for gene in self.genome.genes:
            new_pos_list = gene.act_on_position(self.position)
            if new_pos_list == []: # did not move
                pos = (self.position.x, self.position.y)
            # fetch new position from grid
            for pos in new_pos_list:
                self.grid.grid[pos[0]][pos[1]].has_been_passed = True
                self.distance_traveled += 1
            updated_pos = self.grid.grid[pos[0]][pos[1]]
            self.position = updated_pos

    def reset(self):
        self.position = START_POSITION
        self.grid = Screen(target = TARGET_POSITION, start = START_POSITION)

    def mutate(self):
        # use genome to generate a child genome with a single, random mutation
        new_genome = self.genome.mutate()

        baby = Crawler(genome = new_genome, parent = self)
        return(baby)
    
    def draw_path(self):
        # call grid to draw path
        self.grid.draw_with_crawler_path(self)
        self.reset()

    def score_self(self, method):
        # score Crawler based on configured scoring type
        if method == ScoringTypes.TotalDistance:
            score = self.distance_traveled
        elif method == ScoringTypes.DistanceFromTarget:
            score = 100 - self.position.distance_from(self.grid.target)
        elif method == ScoringTypes.DistanceMinusGenomeLength:
            score = 100 - self.position.distance_from(self.grid.target) - GENOME_LENGTH_PENALTY*len(self.genome.genes)
        elif method == ScoringTypes.TotalDistanceWithGenomePenalty:
            score = self.distance_traveled - GENOME_LENGTH_PENALTY*len(self.genome.genes)
        elif method == ScoringTypes.DistanceToTargetPlusDistanceBonus:
            score = (100 - self.position.distance_from(self.grid.target)) + DISTANCE_TRAVELED_BONUS * self.distance_traveled


        self.score = score
        return(score)
        
# config for start & target positions
TARGET_POSITION = Position(20, 20)
START_POSITION = Position(2, 2)

class Gene:
    # class to hold a gene. has a direction & a distance
    instruction: GeneInstruction
    distance: int

    def __init__(self, instruction=None, distance=None):
        # if no arguments are provided, randomize a gene
        if instruction is None: # randomize
            instruction = random.choice(list(GeneInstruction))
        if distance is None:
            distance = random.randint(0, MAX_DISTANCE)

        self.instruction = instruction
        self.distance = distance

    def mutate(self):
        # mutate random.
        rand1 = random.random()
        rand2 = random.random()
        if rand1 > 0.5:
            # mutate distance 
            if rand2 > 0.5 or self.distance==0:
                self.distance += 1
            else:
                self.distance -= 1

        else:
            # mutate direction
            self.instruction = random.choice(list(GeneInstruction))

    def act_on_position(self, position: Position):
        # given a position, calculates effect of gene on that position
        position_list = []

        for i in range(0, self.distance):
            if self.instruction == GeneInstruction.UP:
                new_pos = position.calc_new_coords(x=-1)
            elif self.instruction == GeneInstruction.DOWN:
                new_pos = position.calc_new_coords(x=1)
            elif self.instruction == GeneInstruction.LEFT:
                new_pos = position.calc_new_coords(y=-1)
            elif self.instruction == GeneInstruction.RIGHT:
                new_pos = position.calc_new_coords(y=1)

            position = Position(x = new_pos[0], y = new_pos[1])

            position_list.append(new_pos)
        return(position_list)
    
    def __str__(self):
        return(f"Gene: {self.instruction} {self.distance} positions")

def main():

    generation = 0
    generations = {}
    bests = {}
    
    # generate initial generation
    generations[0] = []
    for i in range(0, GENERATION_SIZE):
        random_genome = Genome()
        crawler = Crawler(random_genome)
        generations[0].append(crawler)

    os.system('clear')

    for i in range(0, TOTAL_GENERATIONS):
        print(f"GENERATION {i}")
        print("--------------------")


        # each crawler acts & is scored
        for crawler in generations[i]:
            crawler.act()
            crawler.score_self(ScoringTypes.DistanceToTargetPlusDistanceBonus)
            crawler.reset()


        # select winner
        best = -1
        best_crawler = None
        for crawler in generations[i]:
            if crawler.score > best:
                best = crawler.score
                best_crawler = crawler
        
        # act again so we can draw the path of the winner on the screen
        best_crawler.act()
        best_crawler.draw_path()
        best_crawler.reset()

        # print some info about the winner
        print(f"best score is {best_crawler.score}")
        print(f"with a genome of {best_crawler.genome}")
        print(f"final position: {best_crawler.position}")

        # sleep so screen is human readable
        time.sleep(SLEEP_TIME)

        # clear screen if this is not the last generation
        if not i == TOTAL_GENERATIONS-1:
            os.system('clear')

        # save the best crawler for future reference
        bests[i] = best_crawler

        # mutate children
        # but include this gen's winner in next gene
        # (guess he cloned himself)
        generations[i+1] = [Crawler(genome=best_crawler.genome)]
        for n in range(0, GENERATION_SIZE - 1):
            child = best_crawler.mutate()
            generations[i+1].append(child)


main()



