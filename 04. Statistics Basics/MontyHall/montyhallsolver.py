'''
###############################################################################
 ______________________________________
|   ___  ___  ___  ___  ___  _____     |
|  | . || __|| | ||  _|| . ||_   _|    |
|  |   ||__ ||   || |_ |   |  | |      |
|  |_|_||___||_|_||___||_|_|  |_| 2022 |
|______________________________________|

Monty Hall Game Solver
by Edward Larrañaga - 2022

###############################################################################
'''

import random
import matplotlib.pyplot as plt

class Game():
    def __init__(self):
        self.wins = 0
        self.loses = 0
    
    def new_game(self):
        '''
        Defines the initial aprameters for the game
        '''
        self.true_door = random.choice([0,1,2])
        if self.true_door == 0:
            self.empty_A = 1
            self.empty_B = 2
        elif self.true_door == 1:
            self.empty_A = 0
            self.empty_B = 2
        else:
            self.empty_A = 0
            self.empty_B = 1
    
    def choose_door(self):
        self.selected_door = random.choice([0,1,2])
    
    def show_empty_door(self):
        if self.selected_door == self.true_door:
            self.empty_door = self.empty_A
            self.change_door = self.empty_B
        elif self.selected_door == self.empty_A:
            self.empty_door = self.empty_B
            self.change_door = self.true_door
        else:
            self.empty_door = self.empty_A
            self.change_door = self.true_door
    
    def swap_doors(self):
        self.selected_door, self.change_door = self.change_door, self.selected_door
    
    def evaluate_result(self):
        if self.selected_door == self.true_door:
            self.wins += 1
        else:
            self.loses +=1




if __name__ == '__main__':
    N = 10000
    # Swapping Doors
    game_ws = Game()
    for i in range(N):
        game_ws.new_game()
        game_ws.choose_door()
        game_ws.show_empty_door()
        game_ws.swap_doors()
        game_ws.evaluate_result()
    
    # No-Swapping Doors
    game_ns = Game()
    for i in range(N):
        game_ns.new_game()
        game_ns.choose_door()
        game_ns.show_empty_door()
        game_ns.evaluate_result()
    
    fig, ax = plt.subplots(1,2, figsize=(10,6))
    ax[0].bar(['Win','Lose'],[game_ws.wins,game_ws.loses], color=['cornflowerblue','crimson'])
    ax[0].set_ylabel('Attempts')
    ax[0].set_ylim(0,N)
    ax[0].grid()
    ax[0].set_title('Swapping Doors')
    ax[1].bar(['Win','Lose'],[game_ns.wins,game_ns.loses], color=['cornflowerblue','crimson'])
    ax[1].set_ylabel('Attempts')
    ax[1].set_ylim(0,N)
    ax[1].grid()
    ax[1].set_title('No Swapping Doors')
    plt.show()


