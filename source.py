# Liner-Quadratic Regulator (LQR)
# Coded by Takuro Tokunaga
# Last updated: March 21st 2026

# Reference:
# https://qiita.com/taka_horibe/items/5d053cdaaf4c886d27b5
# https://qiita.com/trgkpc/items/8210927d5b035912a153

import numpy as np
import time
import scipy.linalg as linalg
from scipy.integrate import quad

np.set_printoptions(precision=2)

class Comments:
    def __init__(self, comment1, comment2):
        self.com1 = comment1
        self.com2 = comment2

    def begin(self):
        print ("### " + str(self.com1)  + "     ###")

    def end(self):
        print ("### " + str(self.com2)  + "   ###")

    pass

class TimeCounter:
    def __init__(self, time_counter):
        self.time_counter = time_counter        
        return None

    def time_return(self):
        return self.time_counter

class LQR:
    # parameters 
    SIZE = 4
    CRITERIA = np.power(10., -7)
    DT = np.power(10., -3)

    def __init__(self):        

        #  matrix and vectors
        self.A = np.zeros((self.SIZE, self.SIZE), dtype='float64')
        self.B = np.zeros((self.SIZE, 1), dtype='float64')
        self.Q = np.zeros((self.SIZE, self.SIZE), dtype='float64')
        self.R = np.array([[1.0]]) # R is scalar but defined as 1 x 1 matrix

        # To be solved:
        self.P = np.zeros((self.SIZE, self.SIZE), dtype='float64')

        # Initialization
        self.A[0][1] = 1.0
        self.A[1][1] = -15.0
        self.A[1][2] = 10.0
        self.A[2][3] = 1.0
        self.A[3][3] = -15.0

        self.B[1] = 10.0
        self.B[3] = 1.0

        self.Q[0][0] = 1.0        
        self.Q[2][2] = 1.0        

        #
        print("### A matrix: ")
        print(self.A)

        print("### B vector: ")
        print(self.B)

        #
        print("### Q matrix: ")
        print(self.Q)

        self.P = self.Q
        self.P_past = self.P

        self.FirstIteration = True        
    
    def Iteration(self):        

        while True:            
            temp1 = np.dot(self.P, self.A)
            temp2 = np.dot(self.A.T, self.P)
            temp3 = -(self.P @ self.B @ linalg.inv(self.R) @ self.B.T @ self.P) # R is skipped because it is 1
            temp4 = self.Q
            self.P = self.P_past + self.DT*(temp1 + temp2 + temp3 + temp4)
            
            # Skip 1st iteration
            if self.FirstIteration:
                self.FirstIteration = False                
                continue

            difference = abs(np.max(self.P) - np.max(self.P_past))
            #print("Difference: " + str(difference))

            if difference < self.CRITERIA:
                break
            
            self.P_past = self.P.copy()

        return self.P
    
    def ArimotoPotter(self):

        print("### Implementing: ")
        
        return 0

def main():    
    ## start
    start = TimeCounter(time.time())
    msg = Comments('Calculation started.', 'Calculation completed.')        
    msg.begin()    
    ##

    lqr = LQR()

    ans = lqr.Iteration()
    print("### P matrix: ")
    print(ans)

    ans = lqr.ArimotoPotter()
    print("### P matrix: ")
    print(ans)

    ## end
    msg.end()                      # message method
    end = TimeCounter(time.time()) # time display
    print(f"### elapsed_time: {end.time_return() - start.time_return():.2f} [sec] ###")
    ##

# --- main routine ---
if __name__ == '__main__':
    main()
