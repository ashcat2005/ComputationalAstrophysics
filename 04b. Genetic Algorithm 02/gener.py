import numpy as np
import random

X = np.array([10,20,3,4])


for i in range(20):
    j = random.choice(range(len(X)))
    print((j+1)//len(X))