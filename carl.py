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


Module: carl.py
--------

Carl, in honor of Carl Jacobi and Carl Friedrich Gauss, is a program 
written in python that implements two of the best known iterative methods 
for solving sparse linear systems: Jacobi method and Gauss-Seidel method.
According to algebraic theory, solving a linear system can be taked back to 
the problem of finding the vector x that satisfies the relation 
			   Ax = b
where A, b are parametrs given in input, called respectively
- A coefficient matrix: describe a particular problem. The files containing matrices
		        can be download from http://math.nist.gov/MatrixMarket, a popular
		        repositories that contain a wide collection of matrices.
			Once download a .mtx file, put it in dataset folder
- b vector of known terms: for our purpose b is initialized with a vector filled with 1

In order to run the program you need to use the syntax
python carl.py <mtx file> [gauss|jacobi] <MAXITER> <TOLL>
For default MAXITER is 100 and TOLL is 1e-4
Examples:

python carl.py pde225.mtx gauss
python carl.py pde225.mtx jacobi 2000 
python carl.py pde225.mtx jacobi 3000 1e-8
"""

import carl_algebra
import carl_IO
import sys


if __name__ == '__main__':
	args = sys.argv[1:]
	if len(args) < 2:
		print('[!] Error: No such arguments')
		exit(-1)
	else:
		filename = args[0]
		method = args[1].upper()
		max_iter = 100
		tollerance = 1e-4		

		if len(args) == 3:
			max_iter = int(args[2])
		if len(args) == 4:
			max_iter = int(args[2])
			tollerance = float(args[3])
			
		matrix = None		
		try:
			matrix = carl_IO.read_mtxfile(filename)
			print('[+] Matrix in dataset/%s loaded correctly'%(filename))
		except:
			print('[!] Error: Unable to read matrix in dataset/%s'%(filename))
			exit(-1)
		
		if method != 'JACOBI' and method != 'GAUSS':
			print('[!] Error: Unkown method %s'%(method))
			exit(-1)

		print('\n[+] Start computing %s method with %d total iterations and tollerance = %.10f'%(method,max_iter,tollerance))
		
		results = carl_algebra.compute(matrix,method,max_iter,tollerance)
		print('[+] Solution: %s'%(str(results['solution'])))
		print('[+] Computed %d iterations on maximum of %d'%(results['iterations'],max_iter))
		print('[+] Elapsed time: %.4f seconds'%(results['elapsed']))
			
