/*
 * boussinesq_matrix.c
 *
 *  Created on: Jan 18, 2010
 *      Author: ecor
 */



#include "turtle.h"
#include "t_utilities.h"
#include "doublevector_utilities.h"
#include "geometry.h"
#include "geometry_attribute.h"
#include "preconditioned_conjugate_gradient.h"
#include "boussinesq_matrix.h"
#define NULL_VALUE_bm -9999.0

DOUBLEBIN *get_t_Matrix_elements(t_Matrix_element operator, polygon_connection_attribute_array *pca,long boundary) {
	/*!
	 *
	 *\date January 2010
	 *\author Emanuele Cordano
	 *
	 *\param operator (t_Matrix_element) - The operator representing the sparse matrix
	 *\param pca ( *polygon_connection_attribute_array) - the polygon connetion ttributes pca
	 *
	 *\return the sparse matrix as a DOUBLEBIN
	 *
	 *
	 */

	LONGVECTOR *index;
	long i,r,j,kp;
	DOUBLEBIN *dbin;
	DOUBLEVECTOR *x_test;

	index=new_longvector(pca->nh);


	for (i=index->nl;i<=index->nh;i++) {
		index->element[i]=1+pca->element[i]->connections->nh;
	}

	x_test=new_doublevector(pca->nh);

	for (i=x_test->nl;i<=x_test->nh;i++) {
		x_test->element[i]=0.0;
	}

	dbin=new_doublebin(index);

	for (r=dbin->index->nl;r<=dbin->index->nh;r++) {
		x_test->element[r]=1.0;
		dbin->element[r][1]=(*operator)(r,x_test);
		x_test->element[r]=0.0;
		for (j=2;j<=dbin->index->element[r];j++) {
			if (pca->element[r]->connections->element[j-1]==boundary) {
				dbin->element[r][j]=NULL_VALUE_bm;
			} else {
				kp=pca->element[r]->connections->element[j-1];
				x_test->element[kp]=1.0;
				dbin->element[r][j]=(*operator)(r,x_test);
				x_test->element[kp]=0.0;
			}
		}

	}

	free_doublevector(x_test);
	free_longvector(index);

	return dbin;
}
