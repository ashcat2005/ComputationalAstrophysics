'''
###############################################################################
 ______________________________________
|   ___  ___  ___  ___  ___  _____     |
|  | . || __|| | ||  _|| . ||_   _|    |
|  |   ||__ ||   || |_ |   |  | |      |
|  |_|_||___||_|_||___||_|_|  |_| 2022 |
|______________________________________|

Monty Hall Game!
by Edward Larrañaga - 2022

###############################################################################

Sounds acknowledgement:

Open Treasure Chest 8 Bit.wav
Mrthenoronha. 
https://freesound.org/people/Mrthenoronha/sounds/519630/

Game Over
wolderado. 
https://freesound.org/people/wolderado/sounds/415096/

###############################################################################
'''


import pygame, sys, random

#-----------------------------------------------------------------------------#
# SCENES

class IntroScene():
    def __init__(self, screen):
        self.screen = screen
        self.background = pygame.Surface([SCREEN_WIDTH, SCREEN_HEIGHT])
        self.background.fill(pygame.Color('white'))

    def update(self):
        # Handling events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
        # Update
        pygame.display.flip()
        self.screen.blit(self.background,[0,0])


class Game():
    def __init__(self, screen):
        self.newgame()
        self.screen = screen
        self.background = pygame.Surface([SCREEN_WIDTH, SCREEN_HEIGHT])
        self.background.fill(pygame.Color('white'))

    def texts(self, tx01, tx02):
        text_font =  pygame.font.Font("freesansbold.ttf", 20)
        win_font =  pygame.font.Font("freesansbold.ttf", 40)
        light_grey = (100,100,100)
        texts = { 'welcome'   : 'Welcome to the Monty Hall Game!',
                  'choose'    : '<< Choose one of the doors (1, 2 or 3) >>',
                  'selection' : 'You have selected door No. '+str(self.selection+1),
                  'continue'  : '<< Press SPACE to continue >>',
                  'empty'      : 'Behind door No. ' + str(self.empty_door+1) + ' there is ... nothing!',
                  'yn'        : '<< Do you want to trade the door? (y/n) >>',
                  'selection2': 'Then, you have selected door No. '+str(self.selection+1),
                  'continue2' : '<< Press SPACE to see the result >>',
                  'prize'     : 'The prize was behind door No. '+str(self.true_door+1),
                  'none'     : ''
                 }
        self.tx01 = text_font.render(texts[tx01], False, light_grey)
        self.tx01_rect = self.tx01.get_rect(center=(SCREEN_WIDTH//2, SCREEN_HEIGHT//4))
        self.tx02 = text_font.render(texts[tx02], False, light_grey)
        self.tx02_rect = self.tx02.get_rect(center=(SCREEN_WIDTH//2, 3*SCREEN_HEIGHT//4))
        if self.selection == self.true_door:
            self.textresult = win_font.render('YOU WIN!', False, light_grey)
        else:
            self.textresult = win_font.render('YOU LOSE!', False, light_grey)
        self.textresult_rect = self.textresult.get_rect(center=(SCREEN_WIDTH//2, 3*SCREEN_HEIGHT//4))

    def newgame(self):
        self.true_door = random.choice([0,1,2])
        print(self.true_door+1)
        # Sprites
        self.doors_group = pygame.sprite.Group()
        door01 = Door(SCREEN_WIDTH//4,SCREEN_HEIGHT//2, '1','e')
        door02 = Door(SCREEN_WIDTH//2,SCREEN_HEIGHT//2, '2','e')
        door03 = Door(3*SCREEN_WIDTH//4,SCREEN_HEIGHT//2, '3','e')
        if self.true_door == 0:
            door01 = Door(SCREEN_WIDTH//4,SCREEN_HEIGHT//2, '1','f')
            self.empty_doorA = 1
            self.empty_doorB = 2
        elif self.true_door == 1:
            door02 = Door(SCREEN_WIDTH//2,SCREEN_HEIGHT//2, '2','f')
            self.empty_doorA = 0
            self.empty_doorB = 2
        else:
            door03 = Door(3*SCREEN_WIDTH//4,SCREEN_HEIGHT//2, '3','f')
            self.empty_doorA = 0
            self.empty_doorB = 1
        self.doors_group.add(door01, door02, door03)
        self.selection = 0
        self.empty_door = 0
        self.texts('welcome','choose')
        self.phase = 1
           
    def update(self):
        pick_phase = {1 : self.phase1,
                      2 : self.phase2,
                      3 : self.phase3,
                      4 : self.phase4,
                      5 : self.phase5
                     }
        pick_phase[self.phase]()
        # Update
        pygame.display.flip()
        self.screen.blit(self.background,[0,0])
        self.doors_group.draw(self.screen)
        self.screen.blit(self.tx01, self.tx01_rect)
        self.screen.blit(self.tx02, self.tx02_rect)

    def phase1(self):
        # Handling events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_1:
                    self.door_selection(0)
                if event.key == pygame.K_2:
                    self.door_selection(1)
                if event.key == pygame.K_3:
                    self.door_selection(2)
    
    def phase2(self):
        # Handling events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    self.empty_selection()

    def phase3(self):
        # Handling events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_y:
                    self.change_selection()
                if event.key == pygame.K_n:
                    self.texts('selection2','continue2')
                    self.phase = 4

    def phase4(self):
        # Handling events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    self.show_results()
    
    def phase5(self):
        # Handling events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    self.newgame()
        self.screen.blit(self.tx01, self.tx01_rect)
        self.screen.blit(self.textresult, self.textresult_rect)


    def door_selection(self, selection):
        self.selection = selection
        self.doors_group.sprites()[self.selection].select_door()
        self.texts('selection','continue')
        self.phase = 2

    def empty_selection(self):
        if self.selection == self.true_door:
            self.empty_door = self.empty_doorA
            self.change_door = self.empty_doorB
        elif self.selection == self.empty_doorA:
            self.empty_door = self.empty_doorB
            self.change_door = self.true_door
        else:
            self.empty_door == self.empty_doorA
            self.change_door = self.true_door
        self.doors_group.sprites()[self.empty_door].current_sprite = 0
        for i in range(50):
            self.doors_group.sprites()[self.empty_door].open_door()
            self.update()
        self.texts('empty', 'yn')
        self.phase = 3
        
    def change_selection(self):
        self.doors_group.sprites()[self.selection].unselect_door()
        self.selection, self.change_door = self.change_door, self.selection
        self.doors_group.sprites()[self.selection].select_door()
        self.texts('selection2','continue2')
        self.phase = 4
    
    def show_results(self):
        if self.selection == self.true_door:
            pygame.mixer.Sound.play(win_sound)
        else:
            pygame.mixer.Sound.play(fail_sound)
        self.doors_group.sprites()[self.selection].current_sprite = 0
        self.doors_group.sprites()[self.change_door].current_sprite = 0
        for i in range(50):
            self.doors_group.sprites()[self.selection].open_door()
            self.doors_group.sprites()[self.change_door].open_door()
            self.update()
        self.texts('prize','none')
        self.phase = 5



#-----------------------------------------------------------------------------#
# SPRITES

class Door(pygame.sprite.Sprite):
    '''
    General sprite class used for the doors in the game
    '''
    def __init__(self, pos_x, pos_y, door_num, state):
        super().__init__()
        self.sprites =[]
        self.sprites.append(pygame.image.load('sprites/'+door_num+'door01.png'))
        self.sprites.append(pygame.image.load('sprites/'+door_num+state+'door02.png'))
        self.sprites.append(pygame.image.load('sprites/'+door_num+state+'door03.png'))
        self.sprites.append(pygame.image.load('sprites/'+door_num+state+'door04.png'))
        self.sprites.append(pygame.image.load('sprites/'+state+'door05.png'))
        self.sprites.append(pygame.image.load('sprites/'+state+'opendoor.png'))
        self.sprites.append(pygame.image.load('sprites/'+door_num+'doorsel.png'))
        self.current_sprite = 0
        self.image = self.sprites[self.current_sprite]
        self.rect = self.image.get_rect(center = (pos_x, pos_y))
    
    def select_door(self):
        self.image = self.sprites[6]

    def unselect_door(self):
        self.image = self.sprites[0]
    
    def open_door(self):
        self.current_sprite += 0.11
        if self.current_sprite > 5:
            self.current_sprite = 5
        self.image = self.sprites[int(self.current_sprite)]
        

#-----------------------------------------------------------------------------#
# MAIN

# Inizialization
pygame.init()
clock = pygame.time.Clock()

# Screen settings
SCREEN_WIDTH = 960
SCREEN_HEIGHT = 720
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
pygame.display.set_caption('Monty Hall Game!')

# Game scenes
intro = IntroScene(screen)
game = Game(screen)

# Sounds
win_sound = pygame.mixer.Sound('sounds/win.wav')
fail_sound = pygame.mixer.Sound('sounds/fail.wav')


#-----------------------------------------------------------------------------#

# Main loop
while True:
    game.update()
    clock.tick(60)