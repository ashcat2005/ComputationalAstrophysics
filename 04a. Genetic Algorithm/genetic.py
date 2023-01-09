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

Genetic Algorithm 01
by Edward Larrañaga - 2023

###########################################################

'''

import random
import string
import pygame, sys


class genetic_algorithm():
    '''
    Genetic Algorithm
    Uses genetic evolution to obtain a target string
    '''
    def __init__(self, target, N, mutation_rate):
        '''
        Generation of the first generation population 
        '''
        self.target = target
        self.n_g = len(self.target)
        self.N = N
        self.mutation_rate = mutation_rate
        self.pop = []
        self.found = False
        self.begin = False

        for _ in range(self.N):
            self.pop.append(self.get_random_string())
        self.selection = self.select()
    
    def get_random_string(self):
        rand_str = ''.join(random.choice(letters) for i in range(self.n_g))
        return rand_str

    def select(self):
        self.mean_score = 0.
        self.max_score = 0.
        self.best_sample = ''
        selection = []
        for sample in self.pop:
            score = 0
            for i in range(self.n_g):
                if sample[i] == self.target[i]:
                    score += 1
            score = 100*score/self.n_g
            self.mean_score += score/self.N
            if score > self.max_score:
                self.max_score = score
                self.best_sample = sample
            for _ in range(round(score)):
                selection.append(sample)  
        if self.max_score == 100:
            self.found = True
            self.begin = False
        return selection

    def generate(self, selection):
        new_pop = []
        for i  in range(self.N):
            sample01 = random.choice(selection)
            sample02 = random.choice(selection)
            new_sample = ''
            for i in range(self.n_g):
                if random.random() < self.mutation_rate:
                    c = random.choice(letters)
                else: 
                    c = random.choice([sample01,sample02])[i]          
                new_sample += c
            new_pop.append(new_sample)
        return new_pop
    
    def update(self):
        self.pop = self.generate(self.selection)
        self.selection = self.select()
        gui.screen_update()
        pygame.display.flip()
        screen.blit(background,[0,0])


class GUI():
    def __init__(self):
        big_font =  pygame.font.Font("freesansbold.ttf", 40)
        text_font =  pygame.font.Font("freesansbold.ttf", 20)
        self.title = big_font.render('Genetic Algorithm', False, 'black')
        self.begin = text_font.render('<< Press SPACE to begin >>', False, 'black')
        self.begin_rect = self.begin.get_rect(center=(SCREEN_WIDTH//2, 19*SCREEN_HEIGHT//20))
        self.best = text_font.render('Best result:', False, 'black')
        self.n_samples = text_font.render('Number of Samples: '+str(ga.N), False, 'black')
        self.mutation_rate = text_font.render('Mutation rate: '+str(ga.mutation_rate), False, 'black')
        
    
    def text_update(self):
        big_font =  pygame.font.Font("freesansbold.ttf", 35)
        text_font =  pygame.font.Font("freesansbold.ttf", 20)
        samples_font =  pygame.font.Font("freesansbold.ttf", 15)
        self.best_sample = big_font.render(ga.best_sample, False, 'black')
        self.generation = text_font.render('Generation: '+ str(generation), False, 'black')
        self.mean_score = text_font.render('Mean Score: '+ str(round(ga.mean_score))+"%", False, 'black')
        self.samples = []
        for i in range(21):
            self.samples.append(samples_font.render(ga.pop[i], False, 'black'))
                   
    def screen_update(self):
        self.text_update()
        screen.blit(self.title, (50,SCREEN_HEIGHT//10))
        screen.blit(self.begin, self.begin_rect)
        screen.blit(self.best, (50, 150))
        screen.blit(self.best_sample, (50, 200))

        screen.blit(self.n_samples, (50, 350))
        screen.blit(self.mutation_rate, (50, 380))

        screen.blit(self.generation, (50, 450))
        screen.blit(self.mean_score, (50, 480))

        for i in range(21):
            screen.blit(self.samples[i], (600, 100+25*i))






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

letters = string.ascii_lowercase + ' ' + string.ascii_uppercase


N = 1000 
target = 'Vainilla y Zeus son novios'
mutation_rate = 0.01

ga = genetic_algorithm(target, N, mutation_rate)
gui = GUI()

#-----------------------------------------------------------------------------#

# Main loop
generation = 1
while generation<10000:
    for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    if ga.found==True:
                        generation = 1
                        ga.__init__(target, N, mutation_rate)
                    ga.begin = True
    
    if ga.begin == True:
        ga.update()
        generation += 1
    else:
        gui.screen_update()
        pygame.display.flip()
        screen.blit(background,[0,0])
    clock.tick(10)
        
