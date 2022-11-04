'''
Linear Regressor showing the parameter
optimization in real time
'''



import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


class LinearRegression():
    '''
    Linear regression class
    '''
    def __init__(self):
        # Initial random parameters
        np.random.seed(413)
        self.a1 = np.random.rand()
        self.a2 = np.random.rand() 
    
    def pri(self):
        print(self.a1)

    def Y(self, x):
        '''
        Function to fit
        '''
        return self.a1 + self.a2*x
    
    def xi2(self,x,y):
        '''
        Cost function
        '''
        delta = np.zeros_like(y)
        delta[:] = self.Y(x[:]) - y[:]
        return np.sum(delta**2)
    
    def fit(self, x, y):
        '''
        Optimization function
        '''
        alpha= 0.0001 # Learning rate
        tol = 1e-13 # Tolerance
        d1 = np.zeros_like(y)
        d2 = np.zeros_like(y)
        
        print('Función de costo inicial = ', self.xi2(x, y))
        xrange= np.linspace(10,20,40)
        plt.figure()
        self.plot(x, y, xrange)
        
        n= 0 # Epochs
        epsilon = 1
        while epsilon>tol and n<200000:
            # Gradients
            d1[:] = self.Y(x[:]) - y[:]
            grad_Xi_a1 = 2*np.sum(d1)
            d2[:] = d1[:]*x[:]
            grad_Xi_a2 = 2*np.sum(d2)

            Xi2_before = self.xi2(x, y)
            self.a1 = self.a1 - alpha*grad_Xi_a1
            self.a2 = self.a2 - alpha*grad_Xi_a2
            Xi2_after = self.xi2(x, y)
            epsilon = abs(Xi2_before - Xi2_after)
            n +=1
        
        print('Función de costo final = ', self.xi2(x, y))
        print('Number of epochs = ',n)
        self.plot(x, y, xrange)
        plt.show()
        return self.a1, self.a2

    def plot(self, x, y, xrange):
    	plt.cla()
    	plt.scatter(x, y,label='observational data')
    	plt.plot(xrange, self.Y(xrange), '--k', label='linear fit')
    	plt.xlabel(r'apparent magnitude')
    	plt.ylabel(r'logarithm of velocity')
    	plt.legend()
    	plt.pause(0.001)
    
def main():
	# Reading and preparing the data
	df = pd.read_csv("data/hubble.csv")
	df['log10_velocity'] = np.log10(df['velocity'])
	X_train = np.asarray(df['mean_m'])
	y_train = np.asarray(df['log10_velocity'])

	lr = LinearRegression()
	a1, a2 = lr.fit(X_train, y_train)
	print('\nThe optimized parameters are')
	print('a1 = ',a1)
	print('a2 = ',a2)



if __name__ == '__main__':
	main()