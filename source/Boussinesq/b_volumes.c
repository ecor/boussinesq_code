/*!
 * \file b_volumes.c
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

extern STRINGBIN *filenames;


extern DOUBLERASTER_MAP *draster;
// extern DOUBLESQUARE_GRID *dsq;
extern DOUBLE_GRID *dgrid; // ec 20100413
extern PARAM *param;
extern FLAG *flag;

extern DOUBLEVECTOR *elevation_bottom_fine; /*! map of the bottom elevation (fine grid) */
extern DOUBLEVECTOR *elevation_bottom_coarse; /*! map of the bottom elevation (coarse grid) */
extern DOUBLEVECTOR *elevation_bottom_flines; /*! map of the bottom elevation defined in the lines of the fine grid*/
extern DOUBLEVECTOR *elevation_bottom_bottom; /*! map of the elvetion of the bottom of he coarse cells */
extern DOUBLEVECTOR *elevation_terrain_fine; /*! map of the surface terrain elevation (fine grid) */
extern DOUBLEVECTOR *elevation_terrain_flines; /*! map of the surface terrain elevation (fine grid lines) */
extern DOUBLEVECTOR *outlet_coefficient; /*!<map of the coefficient of the rating curve in the outlet q_discharge=C*h_sup^m */
extern DOUBLEVECTOR *outlet_coefficient_surf; /*!<map of the coefficient of the rating curve in the outlet (surface flow)  q_discharge=C*h_sup^m */
extern DOUBLEVECTOR *porosity_fine;   /*! map of porosity defined on the pixels of the fine grid */

extern DOUBLEVECTOR *water_surface_elevation; /*!   map of instantaneous water surface on the pixels of the coarse grid  */
extern DOUBLEVECTOR *water_mass_error; /*! map of instantaneous  water mass error   */

//double Ks=10;
double volume (double eta, long i) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date October 2009
	 *
	 * \return the volume of the i-th cell of the i-th cells
	 *
	 */

	double vol=volume_subs(eta,i)+volume_surf(eta,i);

	return vol;

}


double volume_surf(double eta, long i) {
	/*!
	 *\author Emanuele Cordano (modified October 2009)
	 *\date 22 October 2009
	 *
	 *\param eta - (double) value of water-table level
	 *\param i   - (long) index of the (coarse) cell
	 *
	 *\details It calculates the (surface water) volume of the i-th cell in function of a value of a free surface elevation and fine-grid topography,
	 *  according to what is specified in the documentation with detailed algebraic passages.
	 *
	 *
	 *
	 *\return (surface water) water volume contained in a cell
	 */

	double vol=0.0;
	long r,c,kp;
//	long novalue=-99.0; //dgrid->fine->novalue;
	double hval,aval;
	char *function_name="volume_surf";

	if ((i<dgrid->small_polygon_content->index->nl) || (i>dgrid->small_polygon_content->index->nh)) printf("Error in volume function: volume %ld does not exist [%ld,%ld]. \n ",i,dgrid->small_polygon_content->index->nl,dgrid->small_polygon_content->index->nh);

	for(r=NL;r<=dgrid->small_polygon_content->index->element[i];r++) {
	// for(r=dsq->small_content_polygon->element[i]->nrl;r<=dsq->small_content_polygon->element[i]->nrh;r++){
	//	for(c=dsq->small_content_polygon->element[i]->ncl;c<=dsq->small_content_polygon->element[i]->nch;c++){
		//	kp=dsq->small_content_polygon->element[i]->element[r][c];
			kp=dgrid->small_polygon_content->element[i][r];
				//if (kp!=novalue){

					if ((kp<elevation_bottom_fine->nl) || (kp>elevation_bottom_fine->nh)) {
						aval=0.0;
						hval=0.0;
	//					printf ("Error in volume 'small' polygon %ld (of %ld) cannot belong to the big polygon %ld (r=%ld,c=%ld)!! novalue=%ld \n",kp,elevation_bottom_fine->nh,i,r,c,novalue);
						printf ("Error in function %s volume 'small' polygon %ld (of %ld) cannot belong to the big polygon %ld (r=%ld,c=%ld)!! \n",function_name,kp,elevation_bottom_fine->nh,i,r,c);
						stop_execution();
					}else {
					//	aval=dsq->fine->grid->polygons->element[kp]->area2D;
						aval=dgrid->fine->polygons->element[kp]->area2D;
						hval=(eta-elevation_terrain_fine->element[kp])*aval;
					}

					if (hval>0.0){
						vol+=hval;

					}
			//	}
	//	}
	}



	return vol;
}

double volume_subs(double eta, long i) {
	/*!
	 *\author Emanuele Cordano (modified October 2009)
	 *\date 17 April 2009
	 *
	 *\param eta - (double) value of water-table level
	 *\param i   - (long) index of the (coarse) cell
	 *
	 *\details It calculates the (subsurface water) volume of the i-th cell in function of a value of a free surface elevation and fine-grid topography,
	 *  according to what is specified in the documentation with detailed algebraic passages.
	 *
	 *
	 *
	 *\return water volume contained in a cell
	 */

	double vol=0.0;
	long r,c,kp;
//	long novalue=dsq->fine->novalue;
	double hval,aval;
	char *function_name="volume_subs";
	if ((i<dgrid->small_polygon_content->index->nl) || (i>dgrid->small_polygon_content->index->nh)) printf("Error in volume function: volume %ld does not exist [%ld,%ld]. \n ",i,dgrid->small_polygon_content->index->nl,dgrid->small_polygon_content->index->nh);

	for(r=NL;r<=dgrid->small_polygon_content->index->element[i];r++){
//	for(r=dsq->small_content_polygon->element[i]->nrl;r<=dsq->small_content_polygon->element[i]->nrh;r++){
	//	for(c=dsq->small_content_polygon->element[i]->ncl;c<=dsq->small_content_polygon->element[i]->nch;c++){
		//	kp=dsq->small_content_polygon->element[i]->element[r][c];
			kp=dgrid->small_polygon_content->element[i][r];
			//	if (kp!=novalue){

					if ((kp<elevation_bottom_fine->nl) || (kp>elevation_bottom_fine->nh)) {
						aval=0.0;
						hval=0.0;
	//					printf ("Error in volume 'small' polygon %ld (of %ld) cannot belong to the big polygon %ld (r=%ld,c=%ld)!! novalue=%ld \n",kp,elevation_bottom_fine->nh,i,r,c,novalue);
						printf ("Error in %s volume 'small' polygon %ld (of %ld) cannot belong to the big polygon %ld (r=%ld,c=%ld)!!  \n",function_name,kp,elevation_bottom_fine->nh,i,r,c);
						stop_execution();
					}else {
				//		aval=dsq->fine->grid->polygons->element[kp]->area2D*porosity_fine->element[kp];
						aval=dgrid->fine->polygons->element[kp]->area2D*porosity_fine->element[kp];
						hval=(fmin(eta,elevation_terrain_fine->element[kp])-elevation_bottom_fine->element[kp])*aval;
					}

					if (hval>0.0){
						vol+=hval;

					}
		//		}
		//}
	}



	return vol;
}

//double vertical_area()

double vertical_area(double eta, long j) {
	/*!
	 *
	 *\author Emanuele Cordano
	 *\author October 2009
	 *
	 *
	 *\return total vertical area over the j-th line in function of water table eta.
	 *
	 */

	double area_vert=vertical_area_subs(eta,j);

	return area_vert;

}


double vertical_area_subs(double eta, long j) {
/*!
 * \param Emanuele Cordano
 * \date  19 April 2009 (modified October 2009)
 *
 *\param eta - (double) water surface elevation;
 *\param j  - (long) line indicator;
 *
 *\detail It calculates the (subsurface wet) vertical area of the j-th line in function of a value of a free surface elevation and fine-grid topography,
 * according to what is specified in the documentation with detailed algebraic passages.

 *\return vertical wet area over a line between two polygons
 */
	long c,kp;
	double area_vert=0.0;
	long novalue=dgrid->novalue;
	long nl_l=dgrid->fine->lines->nl; // dsq->fine->grid->lines->nl;
	long nh_l=dgrid->fine->lines->nh; //dsq->fine->grid->lines->nh;

	if ((j>dgrid->small_line_content->index->nh) || (j<dgrid->small_line_content->index->nl)) printf("Error in vertical_area_subs function: Line %ld does not exist [%ld,%ld] . \n",j,dgrid->small_line_content->index->nl,dgrid->small_line_content->index->nh);

	for(c=NL;c<=dgrid->small_line_content->index->element[j];c++) {

		kp=dgrid->small_line_content->element[j][c];

		if (((kp<nl_l) || (kp>nh_l)) && (kp!=novalue) ) {  // && (kp!=novalue)
			printf ("Error in vertical area lines %ld (of %ld)  cannot belong to the big line %ld (c=%ld)!!  \n",kp,nh_l,j,c);
			stop_execution();
		} else if (kp!=novalue) {
//		} else  {
		//	area_vert+=fmax(fmin(eta,elevation_terrain_flines->element[kp])-elevation_bottom_flines->element[kp],0.0)*dsq->fine->grid->lines->element[kp]->length2d;
			area_vert+=fmax(fmin(eta,elevation_terrain_flines->element[kp])-elevation_bottom_flines->element[kp],0.0)*dgrid->fine->lines->element[kp]->length2d; //dsq->fine->grid->lines->element[kp]->length2d;

		}
	}

	return area_vert;
}


double vertical_area_surf(double eta, long j) {
/*!
 * \param Emanuele Cordano
 * \date  22 October 2009
 *
 *\param eta - (double) water surface elevation;
 *\param j  - (long) line indicator;
 *
 *\detail It calculates the (surface wet) vertical area of the j-th line in function of a value of a free surface elevation and fine-grid topography,
 * according to what is specified in the documentation with detailed algebraic passages.

 *\return (surface) vertical wet area over a line between two polygons
 */
	long c,kp;
	double area_vert=0.0;
	long novalue=dgrid->novalue;
	long nl_l=dgrid->fine->lines->nl;
	long nh_l=dgrid->fine->lines->nh;

	// if ((j>dgrid->coarse->lines->nh) || (j<dgrid->coarse->lines->nl))
	 if ((j>dgrid->small_line_content->index->nh) || (j<dgrid->small_line_content->index->nl)) printf("Error in vertical_area_surf function: Line %ld does not exist [%ld,%ld] . \n",j,dgrid->fine->lines->nl,dgrid->fine->lines->nh);

	for(c=NL;c<=dgrid->small_line_content->index->element[j];c++) {
//	for(c=1;c<=dsq->small_content_line->index->element[j];c++) {
//		kp=dsq->small_content_line->element[j][c];
		kp=dgrid->small_line_content->element[j][c];
		if (((kp<nl_l) || (kp>nh_l)) && (kp!=novalue) ) { //
			printf ("Error in vertical area lines %ld (of %ld)  cannot belong to the big line %ld (c=%ld)!! ",kp,nh_l,j,c);
			stop_execution();
		} else if (kp!=novalue) {
	//	} else  {
			area_vert+=fmax(eta-elevation_terrain_flines->element[kp],0.0)*dgrid->fine->lines->element[kp]->length2d; //dsq->fine->grid->lines->element[kp]->length2d;

		}
	}

	return area_vert;
}





double wet_area(double eta,long i, double deta) {
/*!
 * \author Emanuele Cordano
 * \date 8 May 2009
 *
 * \param eta - (double) value of water-table level
 * \param i   - (long) index of the (coarse) cell
 * \param deta - (double) water surface delta utilized for derivate (area_horiz)
 *
 *
 * \details It calculates the volume of the i-th cell in function of a value of a free surface elevation and fine-grid topography,
	 according to what is specified in the documentation with detailed algebraic passages.
	 The wet area is the derivative respect to the variable eta of the stored volume calculated by function volume()

	\return wet area of the i-th cell
 */

//	double vol1,vol2;

//	vol1=volume(eta-deta/2.0,i);

//	vol2=volume(eta+deta/2.0,i);

//	return (vol2-vol1)/deta;
/*! modified by Emanuele Cordano on 2*/
	long r,c,kp;
	double aval=0.0;
	double area=0.0;
	long novalue=dgrid->novalue; // dsq->fine->novalue;

	if ((i<dgrid->small_polygon_content->index->nl) || (i>dgrid->small_polygon_content->index->nh)) printf("Error in wet_area function: cell %ld does not exist [%ld,%ld]. \n ",i,dgrid->small_polygon_content->index->nl,dgrid->small_polygon_content->index->nh);

	for(r=NL;r<=dgrid->small_polygon_content->index->element[i];r++){
//	for(r=dsq->small_content_polygon->element[i]->nrl;r<=dsq->small_content_polygon->element[i]->nrh;r++){
	//	for(c=dsq->small_content_polygon->element[i]->ncl;c<=dsq->small_content_polygon->element[i]->nch;c++){
			kp=dgrid->small_polygon_content->element[i][r];  //dsq->small_content_polygon->element[i]->element[r][c];
			if (kp!=novalue){
				if ((kp<elevation_bottom_fine->nl) || (kp>elevation_bottom_fine->nh)) {
					aval=0.0;

					printf ("Error in wet_aera 'small' polygon %ld (of %ld) cannot belong to the big polygon %ld (r=%ld,c=%ld)!! novalue=%ld \n",kp,elevation_bottom_fine->nh,i,r,c,novalue);
					stop_execution();
				}else if (eta>=elevation_terrain_fine->element[kp]){
						aval=dgrid->fine->polygons->element[kp]->area2D;   //dsq->fine->grid->polygons->element[kp]->area2D;
				} else if (eta>=elevation_bottom_fine->element[kp]) {
						aval=dgrid->fine->polygons->element[kp]->area2D*porosity_fine->element[kp];//dsq->fine->grid->polygons->element[kp]->area2D*porosity_fine->element[kp];
				} else {
						aval=0.0;
				}

				if (aval>0.0){
						area+=aval;

				}
			}
	//	}
	}

	return area;
}



double min_elevation(long i) {
/*!
 *
 * \author Emanuele Cordano
 * \date 13 May 2009
 *
 * \param i (long) - index of the coarse cell
 *
 * \return the elevation of the lowest point of the cell
 *
 */

	double min_elevation=1.0e+5; /*! initialization of elevation with a very high double float number */
	long r,c,kp;
	long novalue=dgrid->novalue; // dsq->fine->novalue;


	if ((i<dgrid->small_polygon_content->index->nl) || (i>dgrid->small_polygon_content->index->nh)) printf("Error in min_evation: cell %ld does not exist [%ld,%ld]. \n ",i,dgrid->small_polygon_content->index->nl,dgrid->small_polygon_content->index->nh);

	for(r=NL;r<=dgrid->small_polygon_content->index->element[i];r++){
//	for(r=dsq->small_content_polygon->element[i]->nrl;r<=dsq->small_content_polygon->element[i]->nrh;r++){
	//	for(c=dsq->small_content_polygon->element[i]->ncl;c<=dsq->small_content_polygon->element[i]->nch;c++){
	//		kp=dsq->small_content_polygon->element[i]->element[r][c];
			kp=dgrid->small_polygon_content->element[i][r]; //dsq->small_content_polygon->element[i]->element[r][c];
				if (kp!=novalue){

					if ((kp<elevation_bottom_fine->nl) || (kp>elevation_bottom_fine->nh)) {

						printf ("Error in min_elevation 'small' polygon %ld (of %ld) cannot belong to the big polygon %ld (r=%ld,c=%ld)!! novalue=%ld \n",kp,elevation_bottom_fine->nh,i,r,c,novalue);
						stop_execution();
					}else {
						min_elevation=fmin(elevation_bottom_fine->element[kp],min_elevation);
					}


				}
		//}
	}

	return min_elevation;

}

double q_discharge_from_outlet_subs_line(double eta, long j) {
	/*!
	 * \param Emanuele Cordano
	 * \date 14 October 2009
	 *
	 *\param eta - (double) water surface elevation;
	 *\param j  - (long) line indicator;
	 *
	 *\detail It calculates the discharge  of the j-th line in function of a value of a free surface elevation and fine-grid topography,
	 * according to what is specified in the documentation with detailed algebraic passages.

	 *\return outflow discarge through a boundary line
	 */
		long c,kp;
		double q_outlet=0.0,q_outlet_surf=0.0;
		double p=param->p_outlet;
		double p_surf=param->p_outlet_surf;
		long novalue=dgrid->novalue; // dsq->fine->novalue;
		long nl_l=dgrid->fine->lines->nl;
		long nh_l=dgrid->fine->lines->nh;

		if ((j>dgrid->small_line_content->index->nh) || (j<dgrid->small_line_content->index->nl)) printf("Error in q_discharge_from_outlet: Line %ld does not exist [%ld,%ld] . \n",j,dgrid->small_line_content->index->nl,dgrid->small_line_content->index->nh);

		for(c=NL;c<=dgrid->small_line_content->index->element[j];c++) {
	//	for(c=1;c<=dsq->small_content_line->index->element[j];c++) {
		//	kp=dsq->small_content_line->element[j][c];
			kp=dgrid->small_line_content->element[j][c];
			if (((kp<nl_l) || (kp>nh_l)) && (kp!=novalue) ) {
				printf ("Error in q_discharge_from_outlet lines %ld (of %ld)  cannot belong to the big line %ld (c=%ld)!! novalue=%ld \n",kp,nh_l,j,c,novalue);
				stop_execution();
			} else if (kp!=novalue && (outlet_coefficient->element[kp]!=0.0) && (outlet_coefficient!=NULL) && (outlet_coefficient_surf!=NULL)) {
				// calcolo della portata
			//	area_vert+=fmax(eta-elevation_bottom_flines->element[kp],0.0)*dsq->fine->grid->lines->element[kp]->length2d;
				q_outlet+=pow(fmax(fmin(eta,elevation_terrain_flines->element[kp])-elevation_bottom_flines->element[kp],0.0),p)*dgrid->fine->lines->element[kp]->length2d*outlet_coefficient->element[kp]; // dsq->fine->grid->lines->element[kp]->length2d
				q_outlet_surf+=pow(fmax(eta-elevation_terrain_flines->element[kp],0.0),p_surf)*dgrid->fine->lines->element[kp]->length2d*outlet_coefficient_surf->element[kp]; // dsq->fine->grid->lines->element[kp]->length2d
			}
		}

		return q_outlet;
}

double q_discharge_from_outlet_cell(double eta, long i) {
	/*!
	 *
	 * \date 14 October 2009
	 * \author Emanuele Cordano
	 *
	 * \param eta - (double)
	 * \param i - (index of the cell)
	 *
	 * \details   It calculates the discharge through all boundary edges of the i-th cell in function of a value of a free surface elevation and fine-grid topography,
	 * according to what is specified in the documentation with detailed algebraic passages. It calls q_discharge_from_outlet for each boundary edge!

	 *\return  outflow discarge area from the i-th cell
	 */
	double q_out=0;
	long j;
	long nh_edges=dgrid->coarse->polygons->element[i]->edge_indices->nh;   // dsq->big->grid->polygons->element[i]->edge_indices->nh;
	long nl_edges=dgrid->coarse->polygons->element[i]->edge_indices->nl; //dsq->big->grid->polygons->element[i]->edge_indices->nl;
	long kp,kl;
	long kboundary=dgrid->coarse->boundary_indicator;   //dsq->big->grid->boundary_indicator;

	for(j=nl_edges;j<=nh_edges;j++) {
		kp=dgrid->coarse->links->element[i]->connections->element[j]; //dsq->big->grid->links->element[i]->connections->element[j];
		kl=dgrid->coarse->polygons->element[i]->edge_indices->element[j];    //dsq->big->grid->polygons->element[i]->edge_indices->element[j];
		if ((kp==kboundary) && (outlet_coefficient!=NULL)) q_out+=q_discharge_from_outlet_subs_line(eta,kl);
	}

	return q_out;
}


