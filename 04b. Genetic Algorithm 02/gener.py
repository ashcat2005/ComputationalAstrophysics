import numpy as np
import random

denom = 8   
for i in range(2*denom):
    angle = i*np.pi/denom
    print(np.cos(angle),",", np.sin(angle))