# https://qiita.com/trgkpc/items/8210927d5b035912a153
# Last updated:
# March 29th, 2026

import numpy as np
import scipy.linalg as linalg

class Solver:
    def __init__(self, A, B, Q, R):        
        self.A = A
        self.B = B
        self.Q = Q
        self.R = R
        
    def MatrixOperation(self):
        # 1. Hamilton matrix
        H = np.block([[self.A.T, -self.B @ linalg.inv(self.R) @ self.B.T], [-self.Q , -self.A]])

        # 2. Eigen values
        eigenvalue, w = linalg.eig(H)

        # 3. Sub-matrix
        Y_, Z_ = [], []
        n = len(w[0])//2

        # Sorting of eigen values
        index_array = sorted([i for i in range(2*n)], key = lambda x:eigenvalue[x].real)
    
        # Choose n from the smallest real part
        for i in index_array[:n]:
            Y_.append(w.T[i][:n])
            Z_.append(w.T[i][n:])
        Y = np.array(Y_).T
        Z = np.array(Z_).T

        # 4. You get P matrix
        if linalg.det(Y) != 0:
            return Z @ linalg.inv(Y)
        else:
            print("Warning: Y is not regular matrix. Result may be wrong!")
            return Z @ linalg.pinv(Y)