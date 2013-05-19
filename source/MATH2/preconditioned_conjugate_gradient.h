
/* MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file preconditioned_conjugate_gradient.h

Copyright, 2009 Stefano Endrizzi, Emanuele Cordano, Matteo Dall'Amico and Riccardo Rigon

This file is part of MATH2.
 MATH2 is free software: you can redistribute it and/or modify
    it under the terms of the GNU  General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MATH2 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU  General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 *
 * \file preconditioning.h
 *
 * \author Emanuele cordano
 */

typedef int (*t_Matrix) (DOUBLEVECTOR *x,DOUBLEVECTOR *y);
typedef double (*t_Matrix_element)(long i,DOUBLEVECTOR *eta);


int get_diagonal(DOUBLEVECTOR *diagonal, t_Matrix_element Matrix);

long jacobi_preconditioned_conjugate_gradient_search(long kkm, double epsilon,  DOUBLEVECTOR *x, DOUBLEVECTOR *b, t_Matrix_element funz);
