
/*! MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file preconditioned_conjugate_gradient.c

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
 * \file pre_conditioning.c
 *
 * \author Emanuele Cordano
 *
 *
 */
#include "turtle.h"
#include "t_utilities.h"
#include "doublevector_utilities.h"
#include "preconditioned_conjugate_gradient.h"

#define MAX_VALUE_DIAG 1e-10
#define MAX_REITERTION 1 //50
// CG_CONVERGENCE 1e-10

t_Matrix Functional_M;

int get_diagonal(DOUBLEVECTOR *diagonal, t_Matrix_element Matrix) {
	/*
	 *
	 * \author Emanuele Cordano
	 * \date May 2008
	 *
	 *\param diagonal (DOUBLEVECTOR *) - the diagonal doublevector of the matrix
	 *\param Matrix (t_Matrix) - a matrix from which the diagonal is extracted
	 *
	 *\brief it saved the square root of diagonal of Matrix in a doublevector
	 *
	 *\return 0 in case of success , -1 otherwise
	 *
	 */

	long i;
	DOUBLEVECTOR *x_v;


	x_v=new_doublevector(diagonal->nh);
//	y_v=new_doublevector(diagonal->nh);
	for (i=x_v->nl;i<=x_v->nh;i++) {
		x_v->element[i]=0.0;

	}
	for (i=x_v->nl;i<=x_v->nh;i++) {
			x_v->element[i]=1.0;
			diagonal->element[i]=(*Matrix)(i,x_v);

			x_v->element[i]=0.0;
	}

	free_doublevector(x_v);
	//free_doublevector(y_v);


	return 0;
}

long jacobi_preconditioned_conjugate_gradient_search(long icnt, double epsilon,  DOUBLEVECTOR *x, DOUBLEVECTOR *b, t_Matrix_element funz){

	/*!
	 *\param icnt  - (long) number of reiterations
	 *\param epsilon - (double) required tollerance (2-order norm of the residuals)
	 *\param x     - (DOUBLEVECTOR *) vector of the unknowns x in Ax=b
	 *\param b     - (DOUBLEVECTOR *) vector of b in Ax=b
	 *\param funz  - (t_Matrix_element) - (int) pointer to the application A (x and y doublevector y=A(param)x ) it return 0 in case of success, -1 otherwise.
	 *
	 *
	 *\brief algorithm proposed by Jonathan Richard Shewckuck in http://www.cs.cmu.edu/~jrs/jrspapers.html#cg and http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
	 *
	 * \author Emanuele Cordano
	 * \date June 2009
	 *
	 *\return the number of reitarations
	 */


	double delta,delta_new,alpha,beta,delta0;
	DOUBLEVECTOR *r, *d,*q,*y,*sr,*diag;

	int sl;
	long icnt_max;
	long j;
	double p;

	r=new_doublevector(x->nh);
	d=new_doublevector(x->nh);
	q=new_doublevector(x->nh);
	y=new_doublevector(x->nh);
	sr=new_doublevector(x->nh);
	diag=new_doublevector(x->nh);

	icnt=0;
	icnt_max=x->nh;


	for (j=x->nl;j<=x->nh;j++){
		y->element[j]=(*funz)(j,x);

	}

    get_diagonal(diag,funz);
//    print_doublevector_elements(diag,PRINT);
//    stop_execution();

    delta_new=0.0;

    for (j=y->nl;j<=y->nh;j++) {

    	r->element[j]=b->element[j]-y->element[j];
    	if (diag->element[j]<0.0) {
    		diag->element[j]=1.0;
    		printf("\n Error in jacobi_preconditioned_conjugate_gradient_search function: diagonal of the matrix (%lf) is negative at %ld \n",diag->element[j],j);
    		stop_execution();
    	}
  //  	diag->element[j]=fmax(diag->element[j],MAX_VALUE_DIAG*fabs(r->element[j]));

    	//d->element[j]=r->element[j]/diag->element[j];
    	if (diag->element[j]==0.0) {  //ec 20100315
    		d->element[j]=0.0;
    	} else {
    		d->element[j]=r->element[j]/(diag->element[j]);
    	}

    	delta_new+=r->element[j]*d->element[j];
    }


  //  printf("delta0 =%le", delta0);

  //  double epsilon0=epsilon;
//    double pe=5.0;


	while ((icnt<=icnt_max) && (max_doublevector(r)>epsilon)) {

		delta=delta_new;
	//	s=(* funz)(q,d);
		p=0.0;

		for(j=q->nl;j<=q->nh;j++) {
			q->element[j]=(*funz)(j,d);
			p+=q->element[j]*d->element[j];

		}
		alpha=delta_new/p;
		for(j=x->nl;j<=x->nh;j++) {
			x->element[j]=x->element[j]+alpha*d->element[j];
		}


	    delta_new=0.0;
	    sl=0;
	    for (j=y->nl;j<=y->nh;j++) {
	    	if (icnt%MAX_REITERTION==0) {
					y->element[j]=(*funz)(j,x);
					r->element[j]=b->element[j]-y->element[j];
	    	} else {
					r->element[j]=r->element[j]-alpha*q->element[j];
	    	}
	    	if (diag->element[j]==0.0) { // ec_20100315
	    		sr->element[j]=0.0;
	    		d->element[j]=0.0;
	    	} else {
	    		sr->element[j]=r->element[j]/diag->element[j];
	    	}

	    	delta_new+=sr->element[j]*r->element[j];
/*    	if (((j==y->nl) && sl==0 ) || (sl==1)) {
	    		printf("delta_new =%le (j=%ld)  ",delta_new,j);
	//    		if (delta_new==0.0) sl=1;
	    	}*/
	    }
	    beta=delta_new/delta;
	   // double aa=1.0e-21;
	// ec  Initial residual:   printf("delta_new =%le p=%le alpha=%le beta=%le delta_max=%le\n",delta_new,p,alpha,beta,max_doublevector(r));

		for (j=d->nl;j<=d->nh;j++) {
			 d->element[j]=sr->element[j]+beta*d->element[j];
		}

		icnt++;


	}
	free_doublevector(diag);
	free_doublevector(sr);
	free_doublevector(r);
	free_doublevector(d);
	free_doublevector(q);
	free_doublevector(y);



	return icnt;

}




