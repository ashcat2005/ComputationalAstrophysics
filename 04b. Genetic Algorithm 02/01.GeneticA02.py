'''
###########################################################
 _______________________________________________
|        _         _      ____      _           |
|       / \    ___| |__  / ___|__ _| |_   2023  |
|      / _ \  / __| '_ \| |   / _` | __|        |
|     / ___ \_\__ \ | | | |__| (_| | |_         |
|    /_/   \_\____/_| |_|\____\__,_|\__|        |
|_______________________________________________|

###########################################################

Genetic Algorithm 02
by Edward Larrañaga - 2023

###########################################################

'''

import numpy as np
import pygame, sys, random

def generate_points(N_points):
    '''
    Generates a set of N_points with random coordinates
    in the region with -1<x<1 and -1<y<1
    '''
    return np.random.rand(N_points,2)*2 - 1


class GeneticAlgorithm():
    def __init__(self, N, selection_rate, mutation_rate, file, N_points=10):
        '''
        Inizialates the class creating the points and the first 
        generation path
        '''
        self.N = N
        self.selection_rate = selection_rate
        self.mutation_rate = mutation_rate
        if file:
            self.points = np.genfromtxt(file, delimiter=',')
            self.N_points = len(self.points)
        else:
            self.points = generate_points(N_points)
            self.N_points = N_points
        self.s_points = self.scale_points(self.points)
        self.points_list = list(range(1,self.N_points))
        self.best_score = 0
        self.top_score = 0
        self.best_path = None
        self.top_path = None
        self.gen_paths()
        self.sample = self.select_sample()
        
    def scale_points(self, points):
        '''
        Scales positions of points to plot
        '''
        scaled_points = np.copy(points)
        scaled_points[:,0] = 150*scaled_points[:,0] + 250 
        scaled_points[:,1] = 150*scaled_points[:,1] + 350
        return scaled_points 
    
    def gen_paths(self):
        '''
        Generates a random closed path
        '''
        self.paths = []
        for _ in range(self.N):
            random.shuffle(self.points_list)
            self.paths.append([0]+ self.points_list +[0])
    
    def distance(self, p1,p2):
        '''
        Calculates the distance berween two points
        '''
        return np.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)

    def calc_scores(self, path):
        '''
        Calculates the score of a path as the sum of the 
        distances between all points
        '''
        score = 0
        for i in range(self.N_points): 
            score += self.distance(self.points[path[i]], self.points[path[i+1]])
        return score
    
    def select_sample(self):
        '''
        Selects the sample with the best paths in the generation
        '''
        number_selected_samples = int(self.N*self.selection_rate)
        unsorted_selection = {}
        for path in self.paths:
            score = 100/self.calc_scores(path)
            unsorted_selection[score] = path
        
        sorted_selection = dict(sorted(unsorted_selection.items(), reverse=True))
        self.best_score = list(sorted_selection.keys())[0]
        self.best_path = sorted_selection[self.best_score]

        if self.best_score > self.top_score:
            self.top_score = self.best_score
            self.top_path = self.best_path
        
        sorted_scores = list(sorted_selection.keys())[0:number_selected_samples]
        sorted_samples = list(sorted_selection.values())[0:number_selected_samples]
        sum_scores = sum(sorted_scores)
        selection = []
        for i in range(len(sorted_scores)): 
            sorted_scores[i] = int(sorted_scores[i]*100/sum_scores)
            for _ in range(sorted_scores[i]):
                selection.append(sorted_samples[i])
        return selection
    
    def crossover(self, selection):
        '''
        Generates the crossover of paths in the selected sample
        '''
        new_paths = []
        for _ in range(self.N):
            sample01 = random.choice(selection)[1:-1]
            sample02 = random.choice(selection)[1:-1]
            point_1 = random.choice(range(len(sample01)))
            point_2 = random.choice(range(point_1,len(sample01)))
            heritance = sample01[point_1: point_2]
            new_sample = []
            i,j = 0,0
            while i < point_1:
                if sample02[j] not in new_sample:
                    new_sample.append(sample02[j])
                    i += 1
                    j += 1
                else:
                    j += 1
            new_sample += heritance
            i = point_2 
            while i < len(sample01):
                if sample02[j] not in new_sample:
                    new_sample.append(sample02[j])
                    i += 1
                    j += 1
                else:
                    j += 1
            for k in range(len(new_sample)):
                if random.random() < self.mutation_rate:
                    y = random.choice(new_sample)
                    new_sample[k],y = y,new_sample[k]
            new_sample = [0] + new_sample + [0]
            new_paths.append(new_sample)
        self.paths = new_paths


    def train(self):
        '''
        Trains the genetic algorithm
        '''
        self.gen_paths()
        self.sample = self.select_sample()
        self.crossover(self.sample)




class GUI():
    def __init__(self):
        self.reset = False
        big_font =  pygame.font.Font("freesansbold.ttf", 40)
        text_font =  pygame.font.Font("freesansbold.ttf", 20)
        self.title = big_font.render('Genetic Algorithm', False, 'black')
        self.begin = text_font.render('<< Press SPACE to begin >>', False, 'black')
        self.begin_rect = self.begin.get_rect(center=(SCREEN_WIDTH//2, 19*SCREEN_HEIGHT//20))
        self.again = text_font.render('<< Press SPACE to reset >>', False, 'black')
        self.again_rect = self.again.get_rect(center=(SCREEN_WIDTH//2, 19*SCREEN_HEIGHT//20))
        self.n_samples = text_font.render('Number of Samples: '+ str(ga.N), False, 'black')
        self.mutation_rate = text_font.render('Mutation rate: '+ str(ga.mutation_rate), False, 'black')

    def text_update(self):
        '''
        Updates the text labels for the GUI
        '''
        big_font =  pygame.font.Font("freesansbold.ttf", 30)
        text_font =  pygame.font.Font("freesansbold.ttf", 20)
        self.generation = text_font.render('Generation: '+ str(generation), False, 'black')
        self.points = big_font.render('Number of Points: '+ str(ga.N_points), False, 'black')
        self.top_score = big_font.render('Minimum distance: '+ str(round(100/ga.top_score,3)), False, 'black')
    
    def screen_update(self):
        '''
        Updates the complete set of text labels
        '''
        self.text_update()
        screen.blit(self.title, (SCREEN_WIDTH//2,SCREEN_HEIGHT//10))
        if self.reset == False:
            screen.blit(self.begin, self.begin_rect)
        else:
            screen.blit(self.again, self.again_rect)
        screen.blit(self.n_samples, (SCREEN_WIDTH//2, 350))
        screen.blit(self.mutation_rate, (SCREEN_WIDTH//2, 380))
        screen.blit(self.generation, (SCREEN_WIDTH//2, 450))
        screen.blit(self.points, (SCREEN_WIDTH//2, SCREEN_HEIGHT//4))
        screen.blit(self.top_score, (SCREEN_WIDTH//2, SCREEN_HEIGHT//4 + 35))

    def draw_paths(self):
        '''
        Draw the paths to visualize the best result
        '''
        '''for path in ga.sample:
            for i in range(ga.N_points):
                pygame.draw.line(screen, 'gray', 
                                 ga.s_points[path[i]], ga.s_points[path[i+1]])
        '''
        if ga.best_path:
            for i in range(ga.N_points):
                pygame.draw.line(screen, 'cornflowerblue', 
                                 ga.s_points[ga.best_path[i]], ga.s_points[ga.best_path[i+1]])

        if ga.top_path:
            for i in range(ga.N_points):
                pygame.draw.line(screen, 'crimson', 
                                 ga.s_points[ga.top_path[i]], ga.s_points[ga.top_path[i+1]],
                                 width=2)

    def draw_points(self):
        '''
        Draws the points
        '''
        for point in ga.s_points: 
            pygame.draw.circle(screen, 'red', point, 3)
        pygame.draw.circle(screen, 'blue', ga.s_points[0], 3)
    
    def update(self):
        '''
        Updates the GUI
        '''
        screen.blit(background,[0,0])
        self.draw_paths()
        self.draw_points()
        self.screen_update()
        pygame.display.flip()


#-----------------------------------------------------------------------------#
# MAIN

# Inizialization
pygame.init()
clock = pygame.time.Clock()

# Screen settings
SCREEN_WIDTH = 960
SCREEN_HEIGHT = 720
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
pygame.display.set_caption('Genetic Algorithm')

background = pygame.Surface([SCREEN_WIDTH, SCREEN_HEIGHT])
background.fill(pygame.Color('white'))

file = 'points2.csv'
N_points = 15

N = 3000
selection_rate = 0.03
mutation_rate = 0.02


ga =GeneticAlgorithm(N, selection_rate, mutation_rate, file, N_points=N_points)
gui = GUI()

begin = False
generation = 1
max_generations = 1500
#-----------------------------------------------------------------------------#

# Main loop
while True:
    for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    if not gui.reset:
                        begin = True
                    else:
                        ga.__init__(N, selection_rate, mutation_rate, file, N_points=N_points)
                        generation = 1
                        gui.reset = False
                
    if not begin:
        gui.update()
        clock.tick(60)
    else:
        while generation < max_generations and begin==True:
            for event in pygame.event.get():
                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_q:
                        begin = False
                        gui.reset = True
            ga.train()
            gui.update()
            generation += 1
        begin=False
        gui.reset = True
        clock.tick(60)