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


Module: carl_IO.py
--------

This module contain the method read_mtxfile that read the content of
the file labeled filename and stored it in a sparse matrix using 
the compress sparse column format (csc format), a popular format for
storing sparse matrices
"""

from scipy import io


def read_mtxfile(filename):
	mtx = io.mmread('dataset/' + filename)
	return mtx.asformat('csc')
