/*!
 * \file b_solver.c
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
#include "t_datamanipulation.h"
#include "t_utilities.h"
#include "rw_maps.h"

#include "doublevector_utilities.h"
#include "preconditioned_conjugate_gradient.h"
#include "keywords_file_b.h"

#include "geometry.h"
#include "g_raster2plvector.h"
#include "bigcells2.h"
#include "geometry2.h"
#include "b_utilities.h"

#include "b_volumes.h"
#include "b_sources.h"
#include "b_v_advection.h"
#include "b_solver.h"

#define MIN_T_DIAG_NEGATIVE_INDEX -20

#define MIN_T_DIAG 1.0e-10
extern STRINGBIN *filenames;
extern DOUBLESQUARE_GRID *dsq;
extern DOUBLE_GRID *dgrid; // ec 20100413
extern DOUBLERASTER_MAP *draster;
extern char *wpath;

extern PARAM *param;
extern FLAG *flag;
extern DOUBLEVECTOR *elevation_bottom_fine; /*! map of the bottom elevation (fine grid) */
extern DOUBLEVECTOR *elevation_bottom_coarse; /*! map of the bottom elevation (coarse grid) */
extern DOUBLEVECTOR *elevation_bottom_flines; /*! map of the bottom elevation defined in the lines of the fine grid*/
extern DOUBLEVECTOR *elevation_bottom_bottom; /*! map of the elvetion of the bottom of he coarse cells */
extern DOUBLEVECTOR *porosity_fine;   /*! map of porosity defined on the pixels of the fine grid */
extern DOUBLEVECTOR *outlet_coefficient; /*!<map of the coefficient of the rating curve in the outlet q_discharge=C*h_sup^m */
extern DOUBLEVECTOR *water_surface_elevation; /*!   map of instantaneous water surface on the pixels of the coarse grid  */
extern DOUBLEVECTOR *water_mass_error; /*! map of instantaneous  water mass error   */
extern DOUBLEVECTOR *water_depth;

extern DOUBLEVECTOR *surface_water_velocity; /*! map of the surface velocity defined in the lines (velocity is positive when going kp1 to kp2 polygons with kp1>kp2 )  */
extern DOUBLEVECTOR *F1_wet_vert_area; /*! coefficient F1 for surface flow as reported in Casulli,2008 */
DOUBLEVECTOR *eta_v, *sources, *dirichlet, *t_diagonal,*t_diagonal_no_dirichlet; /*! eta_v and sources and dirichlet  are  defined as a global value, t_diagonal is the diagonal of t_st_element tensor */
S_TIMES *s_times, *dirichlet_times;
//double Ks=10;

double t_st_operator_element(long i,DOUBLEVECTOR *eta){
	/*!
	 *
	 *
	 * \author Emanuele Cordano
	 * \date October 2009
	 */
	double cond_dirichlet=MAX_ELEVATION_VALUE;

	return (t_st_operator_element_subs(i,eta,cond_dirichlet)+t_st_advection_operator_element(i,eta,cond_dirichlet)); // ec 20100517 add double cond_dirichlet // corrected by ec 20100321 +1.0e-10*eta->element[i]);// (thirdtrial_2) corrected by ec 20100319
}

double t_st_operator_element_no_dirichlet(long i,DOUBLEVECTOR *eta){
	/*!
	 *
	 *
	 * \author Emanuele Cordano
	 * \date May 2010
	 */
	double cond_dirichlet=param->null_dirichlet;

	return (t_st_operator_element_subs(i,eta,cond_dirichlet)+t_st_advection_operator_element(i,eta,cond_dirichlet)); // ec 20100517 add double cond_dirichlet // corrected by ec 20100321 +1.0e-10*eta->element[i]);// (thirdtrial_2) corrected by ec 20100319
}

double t_st_operator_element_subs(long i,DOUBLEVECTOR *eta,double cond_dirichlet) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date June 2009
	 *
	 * \param i - index of the polygon at which t_st_operator_element is calculated
	 * \param eta - (DOUBLEVECTOR *) vector of the unknowns
	 * \param cond_dirichlet - (double) (=MAX_ELEVATION_VALUE) in case dirichlet value are neglected or (=param->null_dirichlet) othrwise (added by ec on 20100517)
	 *
	 *\details
	 * This function calculates the i-th element of the following vector
	 * where eta is the vector of the independent variables
	 * \f[ \left[ \mathbf{T} \cdot \eta^{n+1} \right]_{i}= -\Delta t^{n} \, \sum_{j \in
S_i} \, K_S  \, A_j(\eta^n_{i,m(i,j)})\frac{\eta_{m(i,j)}^{n+1}-\eta_{i}^{n+1}}{\delta_{i,m(i,j)}} \f]
	 * Please see the pdf documentation with the algebraic passages for major details.
	 *
	 */
	double sum=0.0;
	double eta_previous;
	long j,je,line_index,polygon2_index;
	double eta1,eta2,t_coeff,coeff,dist;


	for (j=dgrid->coarse->polygons->element[i]->edge_indices->nl;j<=dgrid->coarse->polygons->element[i]->edge_indices->nh;j++) {

			je=j;
			if (dgrid->coarse->links->element[i]->connections->element[je]!=dgrid->coarse->boundary_indicator) {

				line_index=dgrid->coarse->polygons->element[i]->edge_indices->element[je];
				polygon2_index=dgrid->coarse->links->element[i]->connections->element[je];
				if (dirichlet->element[polygon2_index]<=cond_dirichlet) { // ec 20100517
					dist=dgrid->coarse->links->element[i]->d_connections->element[je];
					eta1=water_surface_elevation->element[i];
					eta2=water_surface_elevation->element[polygon2_index];
					eta_previous=water_surface_elevation_mean(eta1,eta2);
					t_coeff=param->Ks*vertical_area_subs(eta_previous,line_index); // MODIFY HERE  param->p_outlet
			//		t_coeff=param->Ks*fmax(vertical_area_subs(eta1,line_index),vertical_area_subs(eta2,line_index));
		//			t_coeff=param->Ks*(vertical_area(eta1,line_index)+vertical_area(eta2,line_index))/2.0; /* modified on 28 July 2008 */
					coeff=-param->dt*t_coeff*(eta->element[polygon2_index]-eta->element[i])/dist;
					sum=sum+coeff;
				}
			}
		}

	return sum;
}
/*
int T_st_operator(DOUBLEVECTOR *y, DOUBLEVECTOR *eta) {
	*!
	 * \author Emanuele Cordano
	 * \date 19 April 2009 modified on 29 june 2009
	 *
	 * \param y - (DOUBLEVECTOR *) output
	 * \param eta - (DOUBLEVECTOR *) eta
	 *
	 *
	 *

	long i;
//	long je;  modified by Emanuele Cordano on 24 Jun 2009
//	double sum,coeff;
//	long i_cell=(eta->nh+1)/2;
//	double t_coeff;
//	double dist,eta1,eta2;
//	long line_index,polygon2_index;

	if ((y->nh!=eta->nh) ||(eta->nh!=dgrid->coare->big->polygons->nh)) {
		printf("Error in T_st_operator function y and etahas different size [%ld,%ld expected: %ld!\n]",y->nh,eta->nh,dsq->big->grid->polygons->nh);
		return -1;
	}

	for (i=eta->nl;i<=eta->nh;i++) {
		y->element[i]=y->element[i]+t_st_operator_element(i,eta);
	}

	return 0;

}
*/

/*
int wet_area_operator(DOUBLEVECTOR *y, DOUBLEVECTOR *x) {
	*!
	 *\author Emanuele Cordano
	 *\date 19 April 2009
	 *
	 *	 * \param y - (DOUBLEVECTOR *) output
	 * \param x - (DOUBLEVECTOR *) input
	 *
	 *
	long i;
	double area;


	if ((y->nh!=x->nh) ||(x->nh!=dsq->big->grid->polygons->nh)) {
		printf("Error in wet_area_operator function y and etahas different size [%ld,%ld expected: %ld!\n]",y->nh,x->nh,dsq->big->grid->polygons->nh);
		return -1;
	}

	for (i=x->nl;i<=x->nh;i++)  {
//		vol=volume(eta_v->element[i],i);
		area=wet_area(eta_v->element[i],i,param->deta);
		y->element[i]=y->element[i]+area*x->element[i];
	}

	return 0;
}
*/
/*int volume_operator(DOUBLEVECTOR *y, DOUBLEVECTOR *eta) {
	*!
	 *\author Emanuele Cordano
	 *\date 19 April 2009
	 *
	 *	 * \param y - (DOUBLEVECTOR *) output
	 * \param eta - (DOUBLEVECTOR *) eta
	 *
	 *
	long i;
	double vol;


	if ((y->nh!=eta->nh) ||(eta->nh!=dsq->big->grid->polygons->nh)) {
		printf("Error in volume_operator function y and eta have different size [%ld,%ld expected: %ld!\n]",y->nh,eta->nh,dsq->big->grid->polygons->nh);
		return -1;
	}

	for (i=eta->nl;i<=eta->nh;i++)  {
		vol=0.0;

		vol=volume(eta->element[i],i);
//		printf("\nDebug y[%ld]=%le vol=%le area=%lf \n",i,y->element[i],vol,area);
		y->element[i]=y->element[i]+vol;
//		printf("\nDebug 2 y[%ld]=%le vol=%le \n",i,y->element[i],vol);

	}

	return 0;
}
*/
/*
int volume_operator_minus(DOUBLEVECTOR *y, DOUBLEVECTOR *eta) {
	*!
	 *\author Emanuele Cordano
	 *\date 19 April 2009 modified on 1 June 2009
	 *
	 * \param y - (DOUBLEVECTOR *) output
	 * \param eta - (DOUBLEVECTOR *) eta
	 *
	 *
	 *
	long i;
	double vol,area;


	if ((y->nh!=eta->nh) ||(eta->nh!=dsq->big->grid->polygons->nh)) {
		printf("Error in volume_operator_minus function y and eta have different size [%ld,%ld expected: %ld!\n]",y->nh,eta->nh,dsq->big->grid->polygons->nh);
		return -1;
	}
	double t_inv=param->t-param->dt/2.0;
	//printf("t_inv= %lf \n",t_inv);
	get_sources(t_inv,sources);
	for (i=eta->nl;i<=eta->nh;i++)  {
		vol=volume(eta->element[i],i)+sources->element[i]*dsq->big->grid->polygons->element[i]->area2D*param->dt;
//		vol=0.0;
//		printf("\nDebug y[%ld]=%le vol=%le \n",i,y->element[i],vol);
		y->element[i]=y->element[i]-vol;
//		printf("\nDebug  2 y[%ld]=%le vol=%le \n",i,y->element[i],vol);

	}

	return 0;
}
*/

double b_smatrix_element (long i,DOUBLEVECTOR *x){

	/*!
	 * \author Emanuele Cordano
	 * \date 29 June 2009
	 *
	 *
	 * \param i - index of the polygon at which t_st_operator_element is calculated
	 * \param eta - (DOUBLEVECTOR *) vector of the unknowns
	 *
	 *\details Function inverted and solved with a preconditioned conjugate gradient method obtained as follows:
	 * \f[ \left[ \mathbf{T} \cdot \eta^{n+1} \right]_{i}= -\Delta t^{n} \, \sum_{j \in
S_i} \, K_S  \, A_j(\eta^n_{i,m(i,j)})\frac{x_{m(i,j)}^{n+1}-x_{i}^{n+1}}{\delta_{i,m(i,j)}}+ Ws(^m\eta^{n+1})x_i \f]
	 * which is  the  t_st_operator_element() function plus x_i multiplied the wet area of the i-th pixel. This function has the same signature of t_st_operator_element()
	 * Please see the documentation in pdf with the algebraic passages for further details.
	 *
	 */
	if ((t_diagonal->element[i]!=MIN_T_DIAG_NEGATIVE_INDEX) && (dirichlet->element[i]<=param->null_dirichlet)) { // ec 10100323
		return (t_st_operator_element_no_dirichlet(i,x)+wet_area(eta_v->element[i],i,param->deta)*x->element[i]-t_diagonal_no_dirichlet->element[i]*x->element[i]+t_diagonal->element[i]*x->element[i]); //corrected by ec 20100323 //corrected by ec 20100321 for thridtrial_3
	} else if (dirichlet->element[i]>param->null_dirichlet) {

		return (wet_area(eta_v->element[i],i,param->deta)*x->element[i]);
	} else {
		return (x->element[i]);
	}

}

/*

int b_smatrix(DOUBLEVECTOR *y, DOUBLEVECTOR *x) {
	*!
	 * \author Emanuele Cordano
	 * \date  29 june 2009
	 *
	 * \param y - (DOUBLEVECTOR *) output
	 * \param x - (DOUBLEVECTOR *) input
	 *
	 * \details It fills the vector y using the function b_smatrix_element() for each element of y.
	 *

	long i;

	if ((y->nh!=x->nh) ||(x->nh!=dsq->big->grid->polygons->nh)) {
		printf("Error in b_smatrix function y and x has different size [%ld,%ld expected: %ld!\n]",y->nh,x->nh,dsq->big->grid->polygons->nh);
		return -1;
	}

	for (i=x->nl;i<=x->nh;i++) {
		y->element[i]=b_smatrix_element(i,x);
	}

	return 0;

}
*/
//int b_knownterm_0(DOUBLEVECTOR *be0)
/*int b_knownterm(DOUBLEVECTOR *be) {
	*!
	 *
	 *\author Emanuele Cordano
	 *\date April 2009
	 *
	 *\param be : (DOUBLEVECTOR *)
	 *
	 * /DA METTERE A POSTO!!!
	 *


	int l,s,m,p;

	l=zeros(be);
//	print_doublevector_elements(eta_v,100);
//	print_doublevector_elements(water_surface_elevation,100);

    m=volume_operator(be,eta_v);

    s=volume_operator_minus(be,water_surface_elevation);

	p=T_st_operator(be,eta_v);


	if ( (p==0) && (m==0) && (s==0)) {
	//	printf("b_knownterm exit 0 \n");
		return 0;
	} else {
		printf("b_knownterm exit -1 \n");
	}

	return -1;

}
*/
long Newton_convergence(DOUBLEVECTOR *x_temp,DOUBLEVECTOR *be, DOUBLEVECTOR *be0) { //doublee ?????
	/*!
	 *\author Emanuele Cordano
	 *\date 20 April 2009
	 *
	 *\param x_term (DOUBLEVECTOR *) x_term - doublevector of residuals (unkowns)
	 *\param be (DOUBLEVECTOR *) be - known term of the equation in Newton Itaration Method (balance equation);
	 *\param be0 (DOUBLEVECTOR *) be - known term of the equation in .... (balance equation);
	 *
	 * \details - This functions contains the iterative the Newton-like iterative according to Casulli, 2008 :
\f[
^{m+1}\eta^{n+1}= ^{m}\eta^{n+1}-\left[\mathbf{Ws}(
^{m}\eta^{n+1})+\mathbf{T}\right]^{-1}\left[\mathbf{V}(
^{m}\eta^{n+1})+\mathbf{T}\cdot^{m}\eta^{n+1} -\mathbf{b3}\right]
 \f]
 where $m$ is the reiteration level. Please see the pdf documentation of the analytical passages for the documentation.


	 *
	 */
/*
 * The introduced vector is $\mathbf{Ws}( ^{m}\eta^{n+1})$ is the vector of the wet area for each cell,  in fact from (\ref{def_Vw}), it is:
\begin{equation} \label{def_Vw}
 Ws_{i}(\psi_i) =\frac {\mathrm{d} V_{i}(\eta_i)}{\mathrm{d}  \eta_i}
\end{equation}
The  solution  of (\ref{def_solver}) is unique and can be implemented   \citep{Casulli2008} .
 */
	char *function_name="Newton_convergence";
	double a=0;
	double max_b,max_x;
	long i;
	double x_temp_max,x_temp_adm=param->x_temp_adm;
	long kkm,kkmm=0;
	long time_0=0;
	long time_1=0;
	double elapsed;
	double volume_term=0.0,be_term=0.0;
#ifdef WRITE_ITERATION_NUMBER
	char *LOG_FILE_NAME=join_strings(wpath,"boussinesq_iterations.log");
	long kcnt=0;
	FILE *FILE_LOG;

	FILE_LOG=fopen(LOG_FILE_NAME,"a");
#endif
	double initial_residual;
	long newton_cnt=0;
	/* initialization of known term */

	for(i=eta_v->nl;i<=eta_v->nh;i++) {
		if (dirichlet->element[i]<=param->null_dirichlet) {
			volume_term=volume(eta_v->element[i],i);
			be_term=be0->element[i];
			be->element[i]=volume_term+t_st_operator_element(i,eta_v)-be_term;
			if ((t_diagonal->element[i]<=MIN_T_DIAG) && (volume_term<=MIN_T_DIAG) && (be_term<=MIN_T_DIAG)) { //ec 20100323
				t_diagonal->element[i]=MIN_T_DIAG_NEGATIVE_INDEX;
				be->element[i]=eta_v->element[i]-water_surface_elevation->element[i];
			} //end ec 20100323 PUT AN ELSE CONDITION HERE!!!
		} else {
			be->element[i]=volume(eta_v->element[i],i)-volume(dirichlet->element[i],i); //eta_v->element[i]-dirichlet->element[i];
			//		printf ("be[%ld]=%lf\n",i,be->element[i]);
				//		stop_execution();Initial residual
		}

	}
	initial_residual=max_doublevector(be);
	//printf ("Initial residual: %le",initial_residual);
	//stop_execution();
	/* END initialization of known term */
	double newton_toll=param->max_errors;

	do {
		for (i=eta_v->nl;i<=eta_v->nh;i++) {
		//	if (dirichlet->element[i]>param->null_dirichlet) printf("eta[%ld]=%lf \n",i,eta_v->element[i]);
				//eta_v->element[i]=dirichlet->element[i];
		} // added_by_ec_20100517
	//	zeros(x_temp);//fvector
		for(i=eta_v->nl;i<=eta_v->nh;i++) {
			if (dirichlet->element[i]<=param->null_dirichlet) {
				be->element[i]=volume(eta_v->element[i],i)+t_st_operator_element(i,eta_v)-be0->element[i];//here_is_the errorr ec 20100518
			} else {
				be->element[i]=volume(eta_v->element[i],i)-volume(dirichlet->element[i],i); //eta_v->element[i]-dirichlet->element[i];
			//	printf ("be[%ld]=%lf,eta_v->element[%ld]=%lf \n",i,be->element[i],i,eta_v->element[i]);
			//	stop_execution();
			}

		}

	//	b_knownterm(be);

		kkm=0;
	//	kkmm=kkmm+1;
		newton_cnt++;
		//if (max_doublevector(be)>param->max_errors) {
		if (max_doublevector(be)>newton_toll) {
			time_0=clock();
		//	kkm=conjugate_gradient_search(kkm,param->max_errors,x_temp,be,b_smatrix); /* please uncomment for the unconditioned CG solver */
			//kkm=jacobi_preconditioned_conjugate_gradient_search(kkm,param->max_errors, x_temp,be,(*b_smatrix_element));
			kkm=jacobi_preconditioned_conjugate_gradient_search(kkm,1e-6, x_temp,be,(*b_smatrix_element));
			time_1=clock();
			elapsed=(double)((time_1-time_0)/CLOCKS_PER_SEC);
//		max_residual	printf ("Number of reiterations %ld (execution time : %lf )\n",kkm,elapsed);
#ifdef	WRITE_ITERATION_NUMBER
			fprintf(FILE_LOG,"%ld,",kkm);
#endif
			if (newton_cnt>1) {
				for (i=x_temp->nl;i<=x_temp->nh;i++) {
					if (x_temp->element[i]<0.0) {
					//	printf("WARNING in Newton_covergence: x_temp[%ld] is negative (%le) (be=%le) (eta0=%lf, eta1=%lf) at %ld iteration !! \n",i,x_temp->element[i],be->element[i],water_surface_elevation->element[i],eta_v->element[i],kkmm);
					//stop_execution();
				// controllare  errato!!!
						x_temp->element[i]=0.0;
					}
				}

			}
			for(i=x_temp->nl;i<=x_temp->nh;i++) {

				eta_v->element[i]=eta_v->element[i]-x_temp->element[i];
				//if (dirichlet->element[i]>param->null_dirichlet)  eta_v->element[i]=dirichlet->element[i];

			}
			//printf ("xtemo %le",max_doublevector(x_temp));
			DOUBLEVECTOR* residual=new_doublevector(eta_v->nh);
//			double sum_tmp=0.0;
			for(i=eta_v->nl;i<=eta_v->nh;i++) {
					//	residual->element[i]=be->element[i]-b_smatrix_element(i,x_temp);
						if (dirichlet->element[i]<=param->null_dirichlet) {
							residual->element[i]=volume(eta_v->element[i],i)-be0->element[i]+t_st_operator_element(i,eta_v);
						//	residual->element[i]=volume(eta_v->element[i],i)-be0->element[i]+t_st_operator_element_no_dirichlet(i,eta_v);
					//		if ((be0->element[i]<=0.0) && (t_diagonal->element[i]<1.0e-10)) {
					//			be->element[i]=0.0;
		//						t_diagonal->element[i]=;
					//		}
						} else {
							residual->element[i]=volume(eta_v->element[i],i)-volume(dirichlet->element[i],i);
							//		x_temp->element[i];
							//double residual2=eta_v->element[i]-dirichlet->element[i];
				//			printf("residual[%ld]=%le residual2=%le \n",i,residual->element[i],residual2);
				//			stop_execution();
						}
				//		 sum_tmp=sum_tmp+residual->element[i]+be0->element[i];
				//		if (sum_tmp<=0.0) {
					//		printf("res[%ld]=%le be[%ld]=%le sum=%le\n",i,residual->element[i],i,be0->element[i],sum_tmp);

					//	stop_execution();
					//	}
			}
	//		if (sum_tmp<=0.0) {
			//	printf("sum=%le \n",sum_tmp);
			//if (sum_tmp<=0.0)	{printf("sum_tmp=%f",sum_tmp);stop_execution();}
	//		}
		//	free_doublevector(residual);
			max_b=max_doublevector(residual);
			free_doublevector(residual);
			// ec 20100324 time 00:48
			max_x=max_doublevector(x_temp);
			if (max_x==0.0) {

				max_b=0.0;
			}
			// end ec 20100324 time 00:48
		} else {
	//		kkm=-1;
	//		printf ("Mass balance verified with error=%le (max %le)\n",max_doublevector(be),param->max_errors);
			max_b=0.0;
		}

	//	x_temp_max=max_doublevector(x_temp);


	//	printf("max_residual=%le adm=%le \n",max_b,newton_toll);

		}while (max_b>newton_toll); //(max_b>1.0e-6);

	//while (max_b>param->max_errors);//((kkm>0) || (x_temp_max>x_temp_adm));
#ifdef WRITE_ITERATION_NUMBER
	fprintf(FILE_LOG,"  %ld   \n",newton_cnt);
	fclose(FILE_LOG);
#endif
	return newton_cnt; //  kkmm;
}

int time_loop(short print,int (*write_output)(void *v1, void *v2)){
	/*!
	 *
	 *
	 * \author Emanuele Cordano
	 * \date 22 April 2009
	 *
	 *
	 *\param (short) - print FLuidTurtle Parametesrs
	 *\param int (*write_output)(void *v1, void *v2) - function which writes the output for every requested time steps
	 *
	 *
	 *\return value 0 in case of success, -1 otherwise
	 *\details In this routin the discretized system equation:
	 * \f[ \mathbf{V}(\eta^{n+1})+\mathbf{T}  \cdot \eta^{n+1} = \mathbf{b3} \f]
	 *  is solved. Please see the meanings of the symbols in the pdf  documentation with algebraic passages.
	 *  The nonlinear equation system is solved with a Newtton-like
	 *
	 */
	double t=param->t;
	double t_start=param->t_start;
	double t_end=param->t_end;
	double dt_print=param->dt_print;
	double dt=param->dt;
	double t_inv;
	OUTPUT_FILENAMES *fn;
	long i;
	int s;
	FILE *fd;

	char SSSS[ ]={"SSSS"}; /* string containing th number of the printed map */
	long kks,lt;
	DOUBLEVECTOR *x,*be,*be0;
	x=new_doublevector(water_surface_elevation->nh);
	for (i=x->nl;i<=x->nh;i++) { // ec 20110419
		x->element[i]=0.0;
	} // end ec 20110419
	be=new_doublevector(water_surface_elevation->nh);
	be0=new_doublevector(water_surface_elevation->nh);
	fn=new_output_filenames(print);

	eta_v=new_doublevector(water_surface_elevation->nh);
	sources=new_doublevector(water_surface_elevation->nh);
	dirichlet=new_doublevector(water_surface_elevation->nh);
	t_diagonal=new_doublevector(water_surface_elevation->nh); // ec 20100322
	t_diagonal_no_dirichlet=new_doublevector(water_surface_elevation->nh); // ec 20100518
	for (i=t_diagonal->nl;i<=t_diagonal->nh;i++) {// ec 20100322
		t_diagonal->element[i]=0.0;
		t_diagonal_no_dirichlet->element[i]=0.0;
	}// end ec 20100322

	fd=fopen(filenames->element[I_TIMES_FT_FILE]+1,"r");
	if (fd!=NULL) {
		fclose(fd);
		s_times=get_s_times(filenames->element[I_TIMES_FT_FILE]+1,print);
		t=param->t_start;
		get_sources(t,sources);
	} else {
		for (i=sources->nl;i<=sources->nh;i++) {
			sources->element[i]=0.0;
		}
	}
	fd=fopen(filenames->element[I_DIRICHLETTIMES_FT_FILE]+1,"r");

	if (fd!=NULL) {
			fclose(fd);
			printf("cxxx\n");
			dirichlet_times=get_s_times(filenames->element[I_DIRICHLETTIMES_FT_FILE]+1,print);
			t=param->t_start;
			get_dirichletsnode(t,dirichlet);
	} else {
		printf("%s\n",filenames->element[I_DIRICHLETTIMES_FT_FILE]+1);
		for (i=dirichlet->nl;i<=dirichlet->nh;i++) {
			dirichlet->element[i]=param->null_dirichlet;
		}
	}


	printf("water_table_surface_elevation cells:%ld\n",water_surface_elevation->nh);


	while ((t>=t_start) && (t<=t_end)) {
		t=t+dt;
		param->t=t;
		t_inv=t-dt/2.0;

		if (s_times!=NULL) get_sources(t_inv,sources);
		if (dirichlet_times!=NULL) get_dirichletsnode(t,dirichlet);

		update_F1_wet_vert_area();

		for (i=water_surface_elevation->nl;i<=water_surface_elevation->nh;i++){
			eta_v->element[i]=water_surface_elevation->element[i]; /* modified by Emanuele Cordano on may 20, 2009 */
	//		printf(" %le",sources->element[i]);
	//		stop_execution();
			if (sources->element[i]>0.0) eta_v->element[i]=elevation_bottom_bottom->element[i]+sources->element[i]*dt;  /* modified by Emanuele Cordano on may 11, 2010 */

			if (dirichlet->element[i]>param->null_dirichlet) eta_v->element[i]=dirichlet->element[i];


			if (water_depth==NULL) {
				be0->element[i]=volume(water_surface_elevation->element[i],i)+sources->element[i]*dgrid->coarse->polygons->element[i]->area2D*dt-q_discharge_from_outlet_cell(water_surface_elevation->element[i],i)*dt;

			}else {
				be0->element[i]=(water_depth->element[i]+sources->element[i]*dt)*dgrid->coarse->polygons->element[i]->area2D-q_discharge_from_outlet_cell(water_surface_elevation->element[i],i)*dt;

			}
			/* known term for advection in surface flow */

			be0->element[i]=be0->element[i]+b_advection(i); // modified by Emanuele Cordano on 29 Oct 2009

			/* added by Emanuele Cordano   */
		}
		//printf ("stop here11");
		//stop_execution();
		if (water_depth!=NULL) {
			printf("Warning: volumes at the nodes were calculated according to water depth value. Water depth values are now deleting \n");
			free_doublevector(water_depth);
			water_depth=NULL;

		}
		get_diagonal(t_diagonal,t_st_operator_element); //get digonal of t_st_element ec-20100323
		get_diagonal(t_diagonal_no_dirichlet,t_st_operator_element_no_dirichlet); //get digonal of t_st_element_no_dirichlet ec-20100518
		kks=Newton_convergence(x,be,be0);

		/* update of wet_area_surface_velocity */

		update_velocity(eta_v); // added by Emanuele Cordano on 29 Oct 2009
		for (i=water_surface_elevation->nl;i<=water_surface_elevation->nh;i++){

			water_surface_elevation->element[i]=eta_v->element[i];
			water_mass_error->element[i]=be->element[i];
		}
		if (print==1) printf("Time_loop function: time %le [%le,%le] executed with %ld iterations \n",t,t_start,t_end,kks);
		if ((long)t%(long)dt_print==0) {
			//if (print==1) printf("Time_loop function: time %le [%le,%le] executed with %ld iterations \n",t,t_start,t_end,kks);

			lt=(long)(t/dt_print);
			//		SSSS="0000";

	//		printf(" t= %ld %le \n ",lt,t);
//			stop_execution();
			write_suffix( SSSS, (long)t/(long)dt_print,0);
			s=(*write_output)((void *)fn,(void *)SSSS);
			if (s!=0) printf("Error in time_loop: results at %s was not correctly written (exit %d)",SSSS,s);
		}


	}

	if (dirichlet_times!=NULL) free_s_times(dirichlet_times);
	if (s_times!=NULL) free_s_times(s_times);
	free_doublevector(eta_v);
	free_doublevector(be);
	free_doublevector(be0);
	free_doublevector(x);
	free_doublevector(dirichlet);
	free_doublevector(sources);
	free_doublevector(t_diagonal);
	free_OUTPUT_FILENAMES(fn);
	if (print==1) printf ("Function time_loop was successfully executed!! ");
	//free(SSSS);
	return 0;
}


double water_surface_elevation_mean(double eta1,double eta2) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date June 2010
	 *
	 */

	double val=eta1;

	if (flag->arithmetic_mean0==1) {
		val=(eta1+eta2)/2.0;
	} else {
		val=fmax(eta1,eta2);
	}

	return val;
}


//int write_map_results(void *file_results, void *file_error,void *argument3)
int write_map_results(void *output, void *SSSS) {
	/*
	 *\author Emanuele Cordano
	 *\date 21 April 2009
	 *
	 *\param (void *) file_results - name of files where results are plotted
	 *\param (void *) file_error - name of files where errors are plotted
	 *\param (void *) argument3 - pointer which points to information about the time step
	 *
	 *
	 *
	 */
	//char *filename_results=(char*)file_results;
	//char *filename_error=(char*)file_error;
    //char *SSSS;


    OUTPUT_FILENAMES *fu=(OUTPUT_FILENAMES *)output;
    int l,s,g;
    char *results;

    char *error;
    char *velocity_name;

    results=join_strings(fu->file_result,(char *)SSSS);
    error=join_strings(fu->file_error,(char *)SSSS);
    velocity_name=join_strings(results,"_velocity.ft");

    g=write_doublevector_in_a_ascii_ft_format_single_file(velocity_name,surface_water_velocity);

    l=write_raster_from_doublevector(results,water_surface_elevation,draster->coarse->UV,dsq->big->indices_pixel,draster->coarse->layer[DTM_MASK]);
    s=write_raster_from_doublevector(error,water_mass_error,draster->coarse->UV,dsq->big->indices_pixel,draster->coarse->layer[DTM_MASK]);
    if ((l==0) && (s==0)) {
    	free(results);
    	free(error);
    	return 0;
    } else {
    	printf ("Error in write_map_results l=%d s=%d \n ",l,s);
    }
    free(results);
    free(error);
    free(velocity_name);
   // free(SSSS);
    return -1;
}
/* output are seved and then written in ascii format files through the following functions. */

OUTPUT_FILENAMES *new_output_filenames(short print){
/*
 * \author Emanuele Cordano
 * \data April 2006
 *
 */
	OUTPUT_FILENAMES *fn;

	fn=(OUTPUT_FILENAMES *)malloc((size_t)(sizeof(OUTPUT_FILENAMES)));
	if (!fn) t_error("Fn (OUTPUT_FILENAMES) was not allocated");

	fn->file_result=copy_string(filenames->element[O_COARSEMAP_WATERSURFACE_ELEVATION]+1);
	fn->file_error=copy_string(filenames->element[O_COARSEMAP_WATERMASS_ERROR]+1);

	if (print==1) printf("OUTPUT_FILENAMES was successfully initialized file_results: %s file_error: %s !! \n",fn->file_result,fn->file_error);

	return fn;
}

void free_OUTPUT_FILENAMES(OUTPUT_FILENAMES *fn){
/*
 *
 * \author Emanuele Cordano
 * \date April 2006
 *
 */

//if (fn->SSSS!=NULL) free(fn->SSSS);
if (fn->file_error!=NULL) free(fn->file_error);
if (fn->file_result!=NULL) free(fn->file_result);

free(fn);
}


