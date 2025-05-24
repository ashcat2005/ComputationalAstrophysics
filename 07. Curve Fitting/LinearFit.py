'''
Eduard Larrañaga
Computational Astrophysics 
2020

Linear Fit
'''
import numpy as np

def linear_fit(x_data, y_data, sigma_y, sigma_x = 0., dydx = 0.):
     '''
     ------------------------------------------------------------
     linear_fit(x_data, y_data, sigma_y, sigma_x = 0., dydx = 0.)
     Returns the linear fit coeffcients for a given set of data
     ------------------------------------------------------------
     Arguments:
     x_data: numpy array with x-data
     y_data: numpy array with y-data
     sigma_y: numpy array with error data in the y component 
     sigma_x: numpy array with error data in the x component 
              (optional) 
     dydx: slope or derivative information. This is needed when 
           sigma_x is present (optional)
     -----------------------------------------------------------
     '''
     # Total error
     sigma2 = sigma_y**2 + (dydx*sigma_x)**2
     Si = 1./sigma2
     S = Si.sum()
     xi = x_data * sigma2
     yi = y_data * sigma2
     xi2 = x_data *x_data * sigma2
     xiyi =x_data * y_data * sigma2
     Sx = xi.sum()
     Sy = yi.sum()
     Sx2 = xi2.sum()
     Sxy = xiyi.sum()
     
     # Parameters in the linear fit 
     a_1 = (Sy*Sx2 - Sx*Sxy)/(S*Sx2 - Sx*Sx)
     a_2 = (S*Sxy - Sx*Sy)/(S*Sx2 - Sx*Sx)
     
     # Associated error to the parameters
     sigma_a1 = np.sqrt(Sx2/(S*Sx2 - Sx*Sx))
     sigma_a2 = np.sqrt(S/(S*Sx2 - Sx*Sx))
     
     # Chi_squared parameter
     deltai = sigma2 * (a_1 + a_2*x_data - y_data)**2
     chi2 = deltai.sum()
    
     # R^2 Score
     R2 = 1 - chi2/(y_data.mean())
     return a_1, a_2, sigma_a1, sigma_a2, chi2, R2