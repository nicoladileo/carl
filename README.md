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

A detailed discussion of the iterative methods can be found in the 
following link:

http://www-users.cs.umn.edu/~saad/books.html

https://en.wikipedia.org/wiki/Iterative_method

In order to run the program you need to use the syntax

python carl.py <mtx file> [gauss|jacobi] <MAXITER> <TOLL>

For default MAXITER is 100 and TOLL is 1e-4
Examples:

python carl.py pde225.mtx gauss

python carl.py pde225.mtx jacobi 2000 

python carl.py pde225.mtx jacobi 3000 1e-8
