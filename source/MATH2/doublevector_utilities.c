
/*! MATH2 CONTAINS ALGEBRAIC ROUTINES FOR GEOtop AND OTHER MODELS
MATH2 Version 0.9375 KMackenzie

file doublevector_utilities.c

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
 * \file doublevector_utilities.c
 *
 * \author Emanuele cordano
 */



#include "turtle.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "linear_span.h"

#include "doublevector_utilities.h"

#define MAX_REITERATON 1

#define DELTA_MIN epsilon

int linear_comb_doublevector(DOUBLEVECTOR *result,DOUBLEVECTOR *a, DOUBLEVECTOR *b, double ca, double cb) {
	/*!
	 * \param result - (DOUBLEVECTOR *) result vector (r=a*ca+b*cb)
	 * \param a      - (DOUBLECTOR *) vector
	 * \param b      - (DOUBLECTOR *) vector
	 * \param ca      - coefficient for "a" vector
	 * \param cb      - coefficient for "b" vector
	 *
	 * \author Emanuele Cordano
	 * \date November 2008
	 *
	 *\return if r is solved correctly or -1 otherwise
	 *
	 */
	long i;

	if((result->nh!=a->nh) || (result->nh!=b->nh)  || (b->nh!=a->nh) ) {
		printf("Error: in linear_comb vectors result and a and b do not have the same number of elements!!  /n");
		return -1;
	}

	for (i=a->nl;i<=a->nh;i++){
		result->element[i]=ca*a->element[i]+cb*b->element[i];

	}

	return 0;

}




double max_doublevector(DOUBLEVECTOR *v) {
	/*!
	 * \date March 2009
	 * \author Emanuele Cordano
	 *
	 * \brief maximum value in a doublevector
	 *
	 */
	 long j;
	 double MK=fabs(v->element[v->nl]);

	 for (j=v->nl+1;j<=v->nh;j++){
		 MK=fmax(MK,fabs(v->element[j]));
	 }

	 return MK;
}

double min_doublevector(DOUBLEVECTOR *v) {
	/*!
	 * \date February 2008
	 * \author Emanuele Cordano
	 *
	 * \brief minimum value in a doublevector
	 *
	 */
	 long j;
	 char *function_name="min_doublevector";
	 double MK=fabs(v->element[v->nl]);
	 for (j=v->nl;j<=v->nh;j++) {
		 if (v->element[j]<0.0) {
			 printf("Warning in function %s, the at %ld position is %le (negative) \n",function_name,j,v->element[j]);
		 }
	 }

	 for (j=v->nl+1;j<=v->nh;j++){
		 MK=fmin(MK,fabs(v->element[j]));
	 }

	 return MK;
}




