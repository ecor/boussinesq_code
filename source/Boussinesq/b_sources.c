/*!
 * \file b_sources.c
 *
 * \author Emanuele Cordano
 *
 * \attention
 *
 This file is part of Boussinesq.


    Boussinesq is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Boussinesq is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Boussinesq.  If not, see <http://www.gnu.org/licenses/>.

 */


#include "turtle.h"

#include "t_alloc.h"
#include "t_io.h"
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "rw_maps.h"
#include "geometry.h"
#include "g_raster2plvector.h"
#include "bigcells2.h"
#include "geometry2.h"
#include "b_utilities.h"
#include "b_solver.h"
#include "b_volumes.h"
#include "doublevector_utilities.h"
#include "preconditioned_conjugate_gradient.h"
#include "keywords_file_b.h"
#include "b_sources.h"

extern STRINGBIN *filenames;
extern DOUBLESQUARE_GRID *dsq;
extern DOUBLE_GRID *dgrid; // ec 20100413
extern DOUBLERASTER_MAP *draster;


extern PARAM *param;
extern FLAG *flag;
extern DOUBLEVECTOR *elevation_bottom_fine; /*! map of the bottom elevation (fine grid) */
extern DOUBLEVECTOR *elevation_bottom_coarse; /*! map of the bottom elevation (coarse grid) */
extern DOUBLEVECTOR *elevation_bottom_flines; /*! map of the bottom elevation defined in the lines of the fine grid*/
extern DOUBLEVECTOR *elevation_bottom_bottom; /*! map of the elvetion of the bottom of he coarse cells */
extern DOUBLEVECTOR *porosity_fine;   /*! map of porosity defined on the pixels of the fine grid */
extern DOUBLEVECTOR *water_surface_elevation; /*!   map of instantaneous water surface on the pixels of the coarse grid  */
extern DOUBLEVECTOR *water_mass_error; /*! map of instantaneous  water mass error   */

extern S_TIMES *s_times, *dirichlet_times;
//fra una funzione di lettura tempi suffissi;



int get_sources(double t,DOUBLEVECTOR *sources){
	/*!
	 *
	 *
	 * \author Emanuele Cordano
	 * \date May 2009
	 *
	 *\param t (double) time;
	 *\param sources - (DOUBLEVECTOR *) sources
	 *
	 *
	 */

	static double t0;
	static double t1;
	static long k;
	static DOUBLEVECTOR *s0;
	static DOUBLEVECTOR *s1;
	char *map0=NULL;
	char *map1=NULL; /* modified by Emanuele Cordano on October 2009 */
	long j;



	if ((t<=s_times->times->element[s_times->times->nl]) || (s1==NULL) || (s0==NULL)) {

		k=s_times->times->nl;
		t0=s_times->times->element[k];
		t1=s_times->times->element[k+1];
		map0=join_strings(filenames->element[I_SOURCEMAPSERIES_COARSE]+1,s_times->s_suffixes->element[k]+1);
		if (s0==NULL) s0=read_doublevector_from_raster(A_FLAG,map0,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
		map1=join_strings(filenames->element[I_SOURCEMAPSERIES_COARSE]+1,s_times->s_suffixes->element[k+1]+1);
		if (s1==NULL) s1=read_doublevector_from_raster(A_FLAG,map1,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);

	} else while ((t>t1) && (k<s_times->times->nh)) {


		copy_doublevector(s1,s0);
		t0=t1;
		free_doublevector(s1);
        k++;
        free(map1);
        map1=join_strings(filenames->element[I_SOURCEMAPSERIES_COARSE]+1,s_times->s_suffixes->element[k+1]+1);
        s1=read_doublevector_from_raster(A_FLAG,map1,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
        t1=s_times->times->element[k+1];
	}

	if ((s1->nh!=s0->nh) || (sources->nh!=s0->nh)) printf("Error in get_sources: vector s0 [%ld]and s1 [%ld] and s[%ld] have different size!! \n",s1->nh,s0->nh,sources->nh);
	for(j=sources->nl;j<=sources->nh;j++) {
		sources->element[j]=(t-t0)/(t1-t0)*s1->element[j]+(t1-t)/(t1-t0)*s0->element[j];

	}

	if (map0!=NULL) free(map0);
	if (map1!=NULL) free(map1);

	return 0;
}


S_TIMES  *get_s_times(char *filename,short print) {
	/*!
	 * \author Emanuele Cordano
	 * \date 29 May 2009
	 *
	 *\param filename (char *)
	 *\param print (short)
	 *
	 */
	S_TIMES *s_t;
	FILE *f;
// leggere s_times
	s_t=(S_TIMES *)malloc(sizeof(S_TIMES));
	if (!s_t) t_error("In get_s_times: s_t was not allocated! ");

	f=t_fopen(filename,"r");
	int ifile=read_index(f,print);
	s_t->times=read_doublearray(f,print);
	s_t->s_suffixes=read_stringarray(f,print);

	t_fclose(f);

	/* verifies*/

	if (s_t->times->nh!=s_t->s_suffixes->index->nh) {
		printf("Error in get_sorces: times and s_suffixes has different size %ld ad %ld respectively!! \n",s_t->times->nh,s_t->s_suffixes->index->nh);
	}
	if (s_t->times->nh<2) {
		printf("Error in get_sorces: times  %ld are less than two values!! \n",s_t->times->nh);
	}

	long k;
	for (k=s_t->times->nl;k<s_t->times->nh;k++) {
		if (s_t->times->element[k+1]<s_t->times->element[k]) printf("Error in get_sources: times[%ld]=%le is lower than times[%ld]=%le !! \n",k+1,s_t->times->element[k+1],k,s_t->times->element[k]);
	}
	return s_t;

}

void free_s_times(S_TIMES* s_t) {
	/*!
	 * \author Emanuele Cordano
	 * \date 29 May 2008
	 *
	 */

	if (s_t->times!=NULL) free_doublevector(s_t->times);
	if (s_t->s_suffixes!=NULL) free_stringbin(s_t->s_suffixes);

	if (s_t!=NULL) free(s_t);

}



int get_dirichletsnode(double t,DOUBLEVECTOR *dirichlet){
	/*!
	 *
	 *
	 * \author Emanuele Cordano
	 * \date May 2009
	 *
	 *\param t (double) time;
	 *\param dirichlet_node map  - (DOUBLEVECTOR *) dirichlet
	 *
	 *
	 */

	static double t0,tp0;
	static double t1,tp1;
	static long k;
	static DOUBLEVECTOR *d0;
	static DOUBLEVECTOR *d1;
	char *map0,*map1;
	double tp;
	long j;

	if ((t<=dirichlet_times->times->element[dirichlet_times->times->nl]) || (d1==NULL) || (d0==NULL)) {
//		printf("h %ld t=%lf",dirichlet_times->times->nl,t);

		k=dirichlet_times->times->nl;
		t0=dirichlet_times->times->element[k];
		t1=dirichlet_times->times->element[k+1];
		map0=join_strings(filenames->element[I_DIRICHLETMAPSERIES_COARSE]+1,dirichlet_times->s_suffixes->element[k]+1);
		if (d0==NULL) d0=read_doublevector_from_raster(A_FLAG,map0,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
		map1=join_strings(filenames->element[I_DIRICHLETMAPSERIES_COARSE]+1,dirichlet_times->s_suffixes->element[k+1]+1);
		if (d1==NULL) d1=read_doublevector_from_raster(A_FLAG,map1,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
		tp0=pow(t0,param->exp_dirichlet);
		tp1=pow(t1,param->exp_dirichlet);
	} else while ((t>t1) && (k<dirichlet_times->times->nh)) {


		copy_doublevector(d1,d0);
		t0=t1;
		free_doublevector(d1);
        k++;
        map1=join_strings(filenames->element[I_DIRICHLETMAPSERIES_COARSE]+1,dirichlet_times->s_suffixes->element[k+1]+1);
        d1=read_doublevector_from_raster(A_FLAG,map1,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
        t1=dirichlet_times->times->element[k+1];
        tp0=pow(t0,param->exp_dirichlet);
        tp1=pow(t1,param->exp_dirichlet);
    	//tp0=t0; //ec_20100511
    	//tp1=t1; //ec_20100511

	}

	if ((d1->nh!=d0->nh) || (dirichlet->nh!=d0->nh)) printf("Error in get_sources: vector d0 [%ld]and d1 [%ld] and dirichlet[%ld] have different size!! \n",d1->nh,d0->nh,dirichlet->nh);
	for(j=dirichlet->nl;j<=dirichlet->nh;j++) {
		if ((d0->element[j]<=param->null_dirichlet) || (d1->element[j]<=param->null_dirichlet)) {

			dirichlet->element[j]=param->null_dirichlet;
		} else {

			tp=pow(t,param->exp_dirichlet);
		//	tp=t; //ec_20100511
			dirichlet->element[j]=(tp-tp0)/(tp1-tp0)*d1->element[j]+(tp1-tp)/(tp1-tp0)*d0->element[j];
	//		printf("tp0=%le tp1=%le t1=%le tp=%le val0=%le val1=%le val=%le \n",tp0,tp1,t1,tp,d0->element[j],d1->element[j],dirichlet->element[j]);
	//		stop_execution();
		}

//		printf("s[%ld]=%lf (system %lf) \n",j,sources->element[j],t);
	}

	return 0;
}


