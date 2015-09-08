""" 
Copyright (C) 2015  Nicola Dileo
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>


Module: carl_algebra.py
--------

This module contain the code of two iterative methods for solving
sparse linear systems: method of Jacobi and method of Gauss-Seidel.

A detailed discussion of the iterative methods can be found in the 
following link:
http://www-users.cs.umn.edu/~saad/books.html
https://en.wikipedia.org/wiki/Iterative_method
"""

import numpy as np
import scipy 
import scipy.linalg as algebra
from scipy import sparse
from scipy.sparse import linalg
from time import clock

np.set_printoptions(precision = 3)


def compute(A, method, MAXITER, TOLL):
	results = {}
	b = np.ones(A.get_shape()[0])
	if method == 'JACOBI':
		tic = clock()
		solution,iterations = Jacobi(A,b,MAXITER,TOLL)
		toc = clock()	
		results['solution'] = solution
		results['iterations'] = iterations
		results['elapsed'] = toc - tic
	else:
		tic = clock()
		solution,iterations = GaussSeidel(A,b,MAXITER,TOLL)
		toc = clock() 
		results['solution'] = solution
		results['iterations'] = iterations
		results['elapsed'] = toc - tic

	return results


def Jacobi(A, b, MAXITER, TOLL):    
    n = len(b)
    xk = np.ones(shape = n,dtype = float)
      
    D = sparse.diags(A.diagonal(), 0, format = 'csc',)
    L = sparse.tril(A, format = 'csc')
    U = sparse.triu(A, format = 'csc')     
     
    T = -(linalg.inv(D))*(L+U)
    c = (linalg.inv(D))*b
     
    i = 0
    err = TOLL + 1
    while i < MAXITER and err > TOLL:
        x = T*xk + c
        err = np.linalg.norm(x-xk, 1)/np.linalg.norm(x,1)
        xk = x
        i += 1
       
    return xk, i



def GaussSeidel(A, b, MAXITER, TOLL):
    n = len(b)
    xk = np.ones(shape = n,dtype = float)
     
    D = sparse.diags(A.diagonal(), 0, format = 'csc',)
    L = sparse.tril(A, format = 'csc')
    U = sparse.triu(A, format = 'csc')
    
    T = -(linalg.inv(D+L))* U
    c = (linalg.inv(D+L))* b
    
    i = 0
    err = TOLL + 1
    while i < MAXITER and err > TOLL:
        x = T*xk + c
        err = np.linalg.norm(x-xk, 1)/np.linalg.norm(x,1)
        xk = x
        i += 1

    return xk, i
    


