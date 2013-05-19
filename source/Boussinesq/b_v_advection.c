/*!
 * \file b_v_advection.c
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
#include "t_utilities.h"
#include "rw_maps.h"
#include "geometry.h"
#include "g_raster2plvector.h"
#include "bigcells2.h"
#include "geometry2.h"
#include "b_utilities.h"
#include "b_solver.h"
#include "keywords_file_b.h"
#include "b_volumes.h"
#include "b_v_advection.h"

extern STRINGBIN *filenames;


extern DOUBLERASTER_MAP *draster;
//	extern DOUBLESQUARE_GRID *dsq;
extern DOUBLE_GRID *dgrid; // ec 20100413
extern PARAM *param;
extern FLAG *flag;

extern DOUBLEVECTOR *elevation_bottom_fine; /*! map of the bottom elevation (fine grid) */
extern DOUBLEVECTOR *elevation_bottom_coarse; /*! map of the bottom elevation (coarse grid) */
extern DOUBLEVECTOR *elevation_bottom_flines; /*! map of the bottom elevation defined in the lines of the fine grid*/
extern DOUBLEVECTOR *elevation_bottom_bottom; /*! map of the elevation of the bottom of he coarse cells */
extern DOUBLEVECTOR *elevation_terrain_fine; /*! map of the surface terrain elevation (fine grid) */
extern DOUBLEVECTOR *elevation_terrain_flines; /*! map of the surface terrain elevation (fine grid lines) */
extern DOUBLEVECTOR *outlet_coefficient; /*!<map of the coefficient of the rating curve in the outlet q_discharge=C*h_sup^m */
extern DOUBLEVECTOR *porosity_fine;   /*! map of porosity defined on the pixels of the fine grid */
extern DOUBLEVECTOR *runoff_coefficient; /*!<map of the coefficient of the dissipation runoff coefficient  C*(velocity)^p ,*/

///extern DOUBLEVECTOR *runoff_coefficient; /*!<map of the coefficient of the rating curve in the outlet q_discharge=C*h_sup^m */
extern DOUBLEVECTOR *water_surface_elevation; /*!   map of instantaneous water surface on the pixels of the coarse grid  */
extern DOUBLEVECTOR *water_mass_error; /*! map of instantaneous  water mass error   */


extern DOUBLEVECTOR *surface_water_velocity; /*! map of the surface velocity defined in the lines (velocity is positive when going kp1 to kp2 polygons with kp1>kp2 )  */


extern DOUBLEVECTOR *F1_wet_vert_area; /*! coefficient F1 for surface flow as reported in Casulli,2008 */
extern DOUBLEVECTOR  *dirichlet;

//#define grav param->gravity  // gravity acceleration

#define VELOCITY_INITIALIZATION_VALUE -9999

double F1_coefficient(long j, double velocity, double wet_vert_area) {
	/*!
	 *
	 *\param j - (long) index of the line
	 *\param velocity - (double) eulerian velocity thourg the line
	 *\param wat_vert_area - (double) - water vertical area over the j-th line
	 * \author Emanuele Cordano
	 * \date October 2009
	 *
	 */

	double coef1=wet_vert_area;
	double coef2=coef1+runoff_coefficient->element[j]*pow(fabs(velocity),param->p_runoff-1.0)*param->dt*dgrid->coarse->lines->element[j]->length2d;
//	double coef2=coef1+runoff_coefficient->element[j]*param->dt*dsq->big->grid->lines->element[j]->length2d; // veriation from giuseppe formeta version
	if (coef2==0.0) {
		if (coef1==0.0)  {
			return 1.0;
		} else {
			printf("Error in F1_coefficient: divide by 0  function return 0! \n");
			return 0.0;
		}
	}
//	printf("result is %le",coef1/coef2);
//	stop_execution();

	return coef1/coef2;
}

double b_advection(long i) {
/*!
 *
 *\param j - (long) index of the polygon
 *\param (global)velocity - (DOUBLEVECTOR *) eulerian velocity thourg the line
 *\param (global)wat_vert_area - (DOUBLEVECTOR *) - water vertical area over the j-th line multiplied by F1
 *
 * \author Emanuele Cordano
 * \date October 2009
 */
	double bterm=0.0;
	double dt=param->dt;
	double eta_previous;
	double area_vert;
	long j,kl,kp;
	char *function_name="b_advection()";

	long boundary=dgrid->coarse->boundary_indicator;

	if ((i<dgrid->coarse->polygons->nl) || (i>dgrid->coarse->polygons->nh)) printf("Error in %s : polygon %ld ( %ld %ld) does not exist \n",function_name,i,dgrid->coarse->polygons->nl,dgrid->coarse->polygons->nh);
	if ((surface_water_velocity->nh!=F1_wet_vert_area->nh) || (surface_water_velocity->nh!=dgrid->coarse->lines->nh)) printf("Error in %s: velocity elements %ld F1_wet_area_elements %ld  lines %ld \n",function_name,surface_water_velocity->nh,F1_wet_vert_area->nh,dgrid->coarse->lines->nh);
	for (j=dgrid->coarse->polygons->element[i]->edge_indices->nl;j<=dgrid->coarse->polygons->element[i]->edge_indices->nh;j++) {
		// for each edge of i-th polygon
		kl=dgrid->coarse->polygons->element[i]->edge_indices->element[j];
		kp=dgrid->coarse->links->element[i]->connections->element[j];
		if (kp!=boundary) {
			if (kp==i) {
				printf("Error in %s: kp is equal to i (%ld) \n",function_name,i);
			} else if ((kp<dgrid->coarse->polygons->nl) || (kp>dgrid->coarse->polygons->nh)) {
				printf("Error in %s : polygon %ld (near %ld) ( %ld %ld) does not exist \n",function_name,kp,i,dgrid->coarse->polygons->nl,dgrid->coarse->polygons->nh);
			} else if ((kl<dgrid->coarse->lines->nl) || (kl>dgrid->coarse->lines->nh)) {
				printf("Error in %s : line %ld (of polygon %ld) ( %ld %ld) does not exist \n",function_name,kl,i,dgrid->coarse->lines->nl,dgrid->coarse->lines->nh);
			} else {
				// WARNING: if kp>i u is oriented from i to kp, otherwise from kp to i
				eta_previous=water_surface_elevation_mean(water_surface_elevation->element[i],water_surface_elevation->element[kp]);
				area_vert=vertical_area_surf(eta_previous,kl);
			//	if (i>100) { printf("area_vert: %le",area_vert);
					//	stop_execution();
				//		}

				//mettere aree
				//area_vert=vert
				if (kp>i) {
					bterm=bterm-dt*area_vert*asymmetric_surface_velocity(kl,eta_previous);
				} else {
					bterm=bterm+dt*area_vert*asymmetric_surface_velocity(kl,eta_previous);
				}
				//funzione
			}

		}

	}
	//long

	return bterm;
}

double t_st_advection_operator_element(long i, DOUBLEVECTOR *eta,double cond_dirichlet) {
	/*
	 *
	 *\author Emanuele Cordano, Giuseppe Formetta
	 *\date October 2009
	 *
	 *\param i - (long) index of the polygon \
	 *\param eta - (DOUBLEVECTOR *) water surface elevation
	 *\param cond_dirichlet - (double) (=MAX_ELEVATION_VALUE) in case dirichlet value are neglected or (=param->null_dirichlet) othrwise (added by ec on 20100517)
	 *
	 *\return T tensor (for advection law)

	 */
	double dt=param->dt;
	long j,kp,kl;


	double forcing,val=0;
	long boundary=dgrid->coarse->boundary_indicator;
	char *function_name="t_st_advection_operator_element";
	double eta_previous,area_vert;
	if ((i<dgrid->coarse->polygons->nl) || (i>dgrid->coarse->polygons->nh)) printf("Error in %s : polygon %ld ( %ld %ld) does not exist \n",function_name,i,dgrid->coarse->polygons->nl,dgrid->coarse->polygons->nh);
	if (F1_wet_vert_area->nh!=dgrid->coarse->lines->nh) printf("Error in %s:  F1_wet_area_elements %ld  lines %ld \n",function_name,F1_wet_vert_area->nh,dgrid->coarse->lines->nh);

	for (j=dgrid->coarse->polygons->element[i]->edge_indices->nl;j<=dgrid->coarse->polygons->element[i]->edge_indices->nh;j++) {
			// for each edge of i-th polygon
			kl=dgrid->coarse->polygons->element[i]->edge_indices->element[j];
			kp=dgrid->coarse->links->element[i]->connections->element[j];
			if (kp!=boundary) {
				if (kp==i) {
					printf("Error in %s: kp is equal to i (%ld) \n",function_name,i);
				} else if ((kp<dgrid->coarse->polygons->nl) || (kp>dgrid->coarse->polygons->nh)) {
					printf("Error in %s : polygon %ld (near %ld) ( %ld %ld) does not exist \n",function_name,kp,i,dgrid->coarse->polygons->nl,dgrid->coarse->polygons->nh);
				} else if ((kl<dgrid->coarse->lines->nl) || (kl>dgrid->coarse->lines->nh)) {
					printf("Error in %s : line %ld (of polygon %ld) ( %ld %ld) does not exist \n",function_name,kl,i,dgrid->coarse->lines->nl,dgrid->coarse->lines->nh);
				} else if (dirichlet->element[kp]<=cond_dirichlet){ // added by ec 20100517
					//mettere descrizione del tensore
					forcing=(eta->element[kp]-eta->element[i])/dgrid->coarse->links->element[i]->d_connections->element[j];
					eta_previous=water_surface_elevation_mean(water_surface_elevation->element[i],water_surface_elevation->element[kp]);
					area_vert=vertical_area_surf(eta_previous,kl);
		//			if (i>100) { printf("area_vert: %le",area_vert);
		//			stop_execution();
		//			}
					val=val+dt*area_vert*symmetric_surface_velocity(kl,forcing);

				}

			}

		}

	return val;

}

double symmetric_surface_velocity(long j, double forcing) {
/*!
 *
 * \author Emanuele Cordano
 * \date October 2009
 *
 *\param j - (long) index of the line \
 *\param eta - (DOUBLEVECTOR *) water surface elevation
 *
 */
//	double vel=
	double vel=-param->gravity*param->dt*F1_wet_vert_area->element[j]*forcing; // ci va moltplicat per l'area!!!
//	printf ("%le %le %ld",F1_wet_vert_area->element[j],vel,j);
//	stop_execution();
	return vel;

}

double asymmetric_surface_velocity(long j,double eta_previous) {
	/*!
	 *
	 *\author Emanuele Cordano
	 *\date October 2009
	 *
	 *\param j (long) line index
	 */

//	double area_vert=vertical_area_surf(eta_previous,j); //correggere qui !!!!
	double vel=0.0;

//	if (area_vert>0.0)
	vel=F1_wet_vert_area->element[j]*surface_water_velocity->element[j]; ///area_vert;


	return vel;

}

int update_velocity(DOUBLEVECTOR *eta) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date October 2009
	 *
	 *NOO \param (DOUBLEVECTOR *) - eta water surface
	 *
	 */

	long i,j,kp,kl;
	double forcing,eta_previous;
	double simm_vel,asimm_vel;

	long boundary=dgrid->coarse->boundary_indicator;
	char *function_name="update_velocity";
	DOUBLEVECTOR *velocity_temp;

	velocity_temp=new_doublevector(surface_water_velocity->nh);

	if (eta->nh!=dgrid->coarse->polygons->nh) printf("Error in %s: inconsistent number of eta elements %ld (%ld cells) \n",function_name,eta->nh,dgrid->coarse->polygons->nh);
	for (j=surface_water_velocity->nl;j<=surface_water_velocity->nh;j++) {
		velocity_temp->element[j]=VELOCITY_INITIALIZATION_VALUE;

	}


	for (i=eta->nl;i<=eta->nh;i++) {
		for (j=1;j<=dgrid->coarse->links->element[i]->connections->nh;j++) {
			kp=dgrid->coarse->links->element[i]->connections->element[j];
			kl=dgrid->coarse->polygons->element[i]->edge_indices->element[j];

			if (kp!=boundary) {
				if (kp==i) {
					printf("Error in %s: kp is equal to i (%ld) \n",function_name,i);
					return -1;
				} else if ((kp<dgrid->coarse->polygons->nl) || (kp>dgrid->coarse->polygons->nh)) {
					printf("Error in %s : polygon %ld (near %ld) ( %ld %ld) does not exist \n",function_name,kp,i,dgrid->coarse->polygons->nl,dgrid->coarse->polygons->nh);
					return -1;
				} else if ((kl<dgrid->coarse->lines->nl) || (kl>dgrid->coarse->lines->nh)) {
					printf("Error in %s : line %ld (of polygon %ld) ( %ld %ld) does not exist \n",function_name,kl,i,dgrid->coarse->lines->nl,dgrid->coarse->lines->nh);
					return -1;
				} else if (kp>i){ // fa anche le linnee di frontiera!!!
					// WARNING: if kp>i u is oriented from i to kp, otherwise from kp to i
					forcing=(eta->element[kp]-eta->element[i])/dgrid->coarse->links->element[i]->d_connections->element[j];
					simm_vel=symmetric_surface_velocity(kl,forcing);

					eta_previous=water_surface_elevation_mean(water_surface_elevation->element[i],water_surface_elevation->element[kp]);
					asimm_vel=asymmetric_surface_velocity(kl,eta_previous);
					velocity_temp->element[kl]=simm_vel+asimm_vel;
				}
			} else {
				asimm_vel=asymmetric_surface_velocity(kl,water_surface_elevation->element[i]);
				velocity_temp->element[kl]=asimm_vel;



			}
		}
	}

	for (j=surface_water_velocity->nl;j<=surface_water_velocity->nh;j++) {
		surface_water_velocity->element[j]=velocity_temp->element[j];
		if (surface_water_velocity->element[j]==VELOCITY_INITIALIZATION_VALUE) {
			printf("Error in %s:  surface velocity (integrated on the vertical area) not calculated at line %ld (%ld %ld)function return -1 \n",function_name,j,surface_water_velocity->nl,surface_water_velocity->nh);
			return -1;
		}
	}
	free_doublevector(velocity_temp);
//	for (j=surface_water_velocity->nl;j<=surface_water_velocity->nh;j++) {
//		surface_water_velocity->element[j]=asymmetric_surface_velocity(j)+simm_vel->element[j];
//	}
//	free_doublevector(simm_vel);

	// aggiugere v asimmetrica

	return 0;
}

int update_F1_wet_vert_area() {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \data October 2009
	 *
	 *\param eta - (DOUBLEVECTOR *) water srface pressure
	 *
	 */
	char *function_name="update_F1_wet_vert_area";
	long j,i,kp,kl;
	//long novalue=dsq->big->novalue;
	long boundary=dgrid->coarse->boundary_indicator;
	double eta_previous,area_vert,velocity;

//	printf("stop1a");
//	stop_execution();
	if (water_surface_elevation->nh!=dgrid->coarse->polygons->nh) printf("Error in %s: inconsistent number of eta elements %ld (%ld cells) \n",function_name,water_surface_elevation->nh,dgrid->coarse->polygons->nh);
	for (j=F1_wet_vert_area->nl;j<=F1_wet_vert_area->nh;j++) {
		F1_wet_vert_area->element[j]=VELOCITY_INITIALIZATION_VALUE;
	}


	for (i=water_surface_elevation->nl;i<=water_surface_elevation->nh;i++) {
		for (j=1;j<=dgrid->coarse->links->element[i]->connections->nh;j++) {
			kp=dgrid->coarse->links->element[i]->connections->element[j];
			kl=dgrid->coarse->polygons->element[i]->edge_indices->element[j];
			// mettere qui!!!!
			if (kp!=boundary) {
				if (kp==i) {
					printf("Error in %s: kp is equal to i (%ld) \n",function_name,i);
					return -1;
				} else if ((kp<dgrid->coarse->polygons->nl) || (kp>dgrid->coarse->polygons->nh)) {
					printf("Error in %s : polygon %ld (near %ld) ( %ld %ld) does not exist \n",function_name,kp,i,dgrid->coarse->polygons->nl,dgrid->coarse->polygons->nh);
					return -1;
				} else if ((kl<dgrid->coarse->lines->nl) || (kl>dgrid->coarse->lines->nh)) {
					printf("Error in %s : line %ld (of polygon %ld) ( %ld %ld) does not exist \n",function_name,kl,i,dgrid->coarse->lines->nl,dgrid->coarse->lines->nh);
					return -1;
				} else if (kp>i){
					// WARNING: if kp>i u is oriented from i to kp, otherwise from kp to i
					//forcing=(eta->element[kp]-eta->element[i])/dsq->big->grid->links->element[i]->d_connections->element[j];
					//simm_vel=symmetric_surface_water_velocity(kl,forcing);
// fare anche la parte asimmetrica della valocita' ///
					eta_previous=water_surface_elevation_mean(water_surface_elevation->element[i],water_surface_elevation->element[kp]);
					//asimm_vel=asymmetric_surface_velocity(j,eta_previous);
					//surface_water_velocity->element[kl]=simm_vel+asimm_vel;
					area_vert=vertical_area_surf(eta_previous,kl);
				//	velocity=0.0;
				//	if (area_vert>0.0) velocity=surface_water_velocity->element[kl]/area_vert;
					velocity=surface_water_velocity->element[kl];
					F1_wet_vert_area->element[kl]=F1_coefficient(kl,velocity,area_vert);
				//	double sval=F1_coefficient(kl,velocity,area_vert);
				//	printf(" %le %le %ld\n",F1_wet_vert_area->element[kl],sval,kl);
				//	stop_execution();
				}
			} else {
				area_vert=vertical_area_surf(water_surface_elevation->element[i],kl);
				velocity=0.0;
			//	if (area_vert>0.0)
				velocity=surface_water_velocity->element[kl];
				F1_wet_vert_area->element[kl]=F1_coefficient(kl,velocity,area_vert);

			}
		}
	}



	for (j=F1_wet_vert_area->nl;j<=F1_wet_vert_area->nh;j++) {
		if (F1_wet_vert_area->element[j]==VELOCITY_INITIALIZATION_VALUE) {
			printf("Error in %s:  F1 coefficient (integrated on the vertical area) not calculated at line %ld (%ld %ld)function return -1 \n",function_name,j,F1_wet_vert_area->nl,F1_wet_vert_area->nh);
			return -1;
		}
	}

	return 0;

}
