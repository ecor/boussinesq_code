
/*! \mainpage Boussinesq

The program Boussinesq simulates the dynamics  of water-table surface solved in a hillslope or a small catchment.
The theory is based on the 2D Boussinesq Equation: \f[ s\frac {\partial \eta}  {\partial t} = \nabla \cdot \left[ K_S \, H (\eta,x,y)  \, \nabla \eta \right]+Q \f]
where  \f$ \eta \f$ is the piezometric  elevation (unknown),  \f$ t \f$ is time,  \f$ \nabla \f$ is the space gradient operator,  \f$ H(\eta,x,y) \f$ is the thickness of the aquifer which is a function of \f$ \eta \f$ and space, \f$Q \f$ is a source term which also accounts for boundary conditions, \f$K_S \f$ is the saturated hydraulic conductivity and  \f$s \f$ is porosity .
\warning
The Boussinesq Equation is solved with finite volume numerical methods according to Casulli, 2008 ( <http://www3.interscience.wiley.com/journal/121377724/abstract?CRETRY=1&SRETRY=0> and http://onlinelibrary.wiley.com/doi/10.1002/wrcr.20072/references)
The Maps of distributed quantiaties are distributed as vectors of double float numbers (DOBLEVECTOR data struct type ,in this case)

\author Emanuele Cordano (emanuele.cordano@gmail.com), Riccardo Rigon, Vincenzo Casulli, Stefano Endrizzi, Matteo Dall'Amico

\version 3.2.1

\date 2008-2009 (2013)

\attention
Boussinesq is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
/*!
 * \file main.c
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

#include <sys/stat.h>

#include "turtle.h"
#include "license.h"
#include "get_filenames.h"
#include "rw_maps.h"
#include "linear_span.h"
#include "geometry.h"

#include "geometry_utilities.h"
#include "read_command_line.h"
#include "additional_read_functions.h"
//#include "geometry_io.h"cccc
//#include "geometry_attribute.h"
//#include "geometry_freememory.h"


#include "g_raster2plvector.h"

#include "bigcells2.h"
#include "geometry2.h"
#include "b_solver.h"
#include "b_readgrid.h"
#include "b_volumes.h"
#include "b_utilities.h"
#include "keywords_file_b.h"

#define PRINT print

#define VERBOSE_M "-verbose"
#define CREATE_GRID "-creates_the_grid"
#define WRITE_GRID "-writes_the_grid"

#define MESH_M "-only_mesh"
#define WPATH "-wpath"
#define NOWPATH ""


#define FLOAT_TYPE 0
#define MAP_FORMAT 3

#define BOUNDARY -10
#define EWRES 1
#define NSRES 2

#define INTEGER_NULL -99
#define D_INIT -999.0

#define NO_OPTION "missing_option"
#define ARITHMETIC_MEAN0 "--arithmetic-mean0"
#define MISSING_FILE join_strings(wpath,NO_OPTION)
#define PROGRAM "boussinesq"

#define FILE_LINES "-lines-resume"
#define FILE_POLYGONS "-polygons-resume"
#define FILE_CONNECTIONS "-connections-resume"
#define RESULTS_FILES "-results-rasters"

#define FILE_DATA "data="

#define RFACTOR param->max_error
/* START HEADER FILE *

/* END HEADER FILE */
//char buffer[FILENAME_MAX];


STRINGBIN *filenames;
DOUBLESQUARE_GRID *dsq;
DOUBLE_GRID *dgrid;

DOUBLERASTER_MAP *draster;

 /* all variable which are gridded are declared below as global (dynamic) doublevectors */
DOUBLEVECTOR *elevation_bottom_fine; /*! map of the bottom elevation (fine grid) */
DOUBLEVECTOR *elevation_bottom_coarse; /*! map of the bottom elevation (coarse grid) */
DOUBLEVECTOR *elevation_bottom_flines; /*! map of the bottom elevation defined in the lines of the fine grid*/
DOUBLEVECTOR *elevation_bottom_bottom; /*! map of the elvetion of the bottom of he coarse cells */
DOUBLEVECTOR *elevation_terrain_fine; /*! map of the surface terrain elevation (fine grid) */
DOUBLEVECTOR *elevation_terrain_flines; /*! map of the surface terrain elevation (fine grid lines) */
DOUBLEVECTOR *porosity_fine;   /*! map of porosity defined on the pixels of the fine grid */
DOUBLEVECTOR *outlet_coefficient,*outlet_coefficient_v_fine; /*!<map of the coefficient of the rating curve in the outlet  (subsurface/groundwater flow) q_discharge=C*h_sup^m */
DOUBLEVECTOR *outlet_coefficient_surf,*outlet_coefficient_surf_v_fine;  /*!<map of the coefficient of the rating curve in the outlet (surface flow)  q_discharge=C*h_sup^m */
DOUBLEVECTOR *runoff_coefficient,*runoff_coefficient_v_coarse; /*!<map of the coefficient of the dissipation runoff coefficient  C*(velocity)^p ,*/

DOUBLEVECTOR *water_surface_elevation; /*!   map of instantaneous water surface on the pixels of the coarse grid  */
DOUBLEVECTOR *water_mass_error; /*! map of instantaneous  water mass error   */
DOUBLEVECTOR *water_depth; /* map of water depth */

DOUBLEVECTOR *surface_water_velocity; /*! map of the (wet area) surface velocity defined in the lines (velocity is positive when going kp1 to kp2 polygons with kp1>kp2 )  */
DOUBLEVECTOR *F1_wet_vert_area; /*! coefficient F1 for surface flow as reported in Casulli,2008 */
PARAM *param;
FLAG *flag;
int read_b_parameters_ft(char *filename_data, short print);
char *wpath;

int main(int argc,char *argv[]) {
	/*!
	 *
	 *\author Emanuele Cordano
	 *\date April 2009
	 *
	 */

	short print=read_flag(argc,argv,VERBOSE_M,0);

	flag=(FLAG *)malloc(sizeof(FLAG));
	if (!flag) t_error("flag was not allocated!!");

	flag->arithmetic_mean0=read_flag(argc,argv,ARITHMETIC_MEAN0,print);
	short cgrid=read_flag(argc,argv,CREATE_GRID,print);
	short wgrid=read_flag(argc,argv,WRITE_GRID,print);
	double time_init,time_end;
	long n_cm_layers,n_fm_layers,i,j;

	time_init=clock();

	/* get filanames */
	wpath=read_option_string(argc,argv,WPATH,NOWPATH,print);

	filenames=get_filenames_from_keys(wpath,PROGRAM,print);

		n_cm_layers=N_MAPS; //2;//CM_NLAYERS;
		n_fm_layers=N_MAPS;




	draster=get_doubleraster_map(n_cm_layers,n_fm_layers,filenames->element[I_MORPHO_ELEVATION_COARSE]+1,filenames->element[I_MORPHO_ELEVATION_FINE]+1,print);
//	if (cgrid==1) { //lines commented by Emanuele Cordano on 27-10-2009
	long dc=draster->coarse->layer[DTM_MASK]->nch*2.0;
	long df=draster->fine->layer[DTM_MASK]->nch*2.0; /* temporary modification */
	dsq=get_doublesquare_grid(draster,filenames->element[O_RESUME_FILES]+1,index_pixel_from_a_bin,index_pixel_from_a_bin,dc,df,print);
	if (wgrid==1) write_doublesquare_grid(dsq);
	dgrid=new_double_grid_from_doublesquare_grid(dsq); //ec 20100413
//	} else {
//		dsq=read_doublesquare_grid(draster,filenames->element[O_RESUME_FILES]+1,index_pixel_from_a_bin,index_pixel_from_a_bin,print);
//		if (wgrid==1) write_doublesquare_grid(dsq);
//	}
	/* mask maps for lines */

	draster->coarse->layer[H_MASK]=doublemap_from_longmap(dsq->big->indices_horizontal_lines,dsq->big->novalue,draster->coarse->UV->V->element[2]);
	draster->coarse->layer[V_MASK]=doublemap_from_longmap(dsq->big->indices_vertical_lines,dsq->big->novalue,draster->coarse->UV->V->element[2]);

	draster->fine->layer[H_MASK]=doublemap_from_longmap(dsq->fine->indices_horizontal_lines,dsq->big->novalue,draster->fine->UV->V->element[2]);
	draster->fine->layer[V_MASK]=doublemap_from_longmap(dsq->fine->indices_vertical_lines,dsq->big->novalue,draster->fine->UV->V->element[2]);


	/* abbresses  (indices) */

	 if (wgrid==1) write_cell(filenames->element[O_COARSEMAP_INDEX_CELLS]+1,dsq->big->grid->polygons->nh,draster->coarse->UV,dsq->big->indices_pixel,draster->coarse->layer[DTM_MASK],FLOAT_TYPE,MAP_FORMAT);
	 if (wgrid==1) write_cell(filenames->element[O_FINEMAP_INDEX_CELLS]+1,dsq->fine->grid->polygons->nh,draster->fine->UV,dsq->fine->indices_pixel,draster->fine->layer[DTM_MASK],FLOAT_TYPE,MAP_FORMAT);


	 if (wgrid==1) write_lines(filenames->element[O_COARSEMAP_INDEX_LINES]+1,dsq->big->grid->lines->nh,dsq->big->nhorizontal_lines,draster->coarse->UV,dsq->big->indices_horizontal_lines,dsq->big->indices_vertical_lines,draster->coarse->layer[H_MASK],draster->coarse->layer[V_MASK],FLOAT_TYPE,MAP_FORMAT);
	 if (wgrid==1) write_lines(filenames->element[O_FINEMAP_INDEX_LINES]+1,dsq->fine->grid->lines->nh,dsq->fine->nhorizontal_lines,draster->fine->UV,dsq->fine->indices_horizontal_lines,dsq->fine->indices_vertical_lines,draster->fine->layer[H_MASK],draster->fine->layer[V_MASK],FLOAT_TYPE,MAP_FORMAT);

	/*  allocation and initialization of the gridded variables */



	elevation_bottom_fine=read_doublevector_from_raster(A_FLAG,filenames->element[I_MORPHO_ELEVATION_FINE]+1,draster->fine->layer[DTM_MASK],draster->fine->UV,dsq->fine->indices_pixel);
	if (elevation_bottom_fine==NULL) printf("Error: map corresponding to %s is missing or does not have an acceptable format \n",filenames->element[I_MORPHO_ELEVATION_FINE]+1);

	elevation_bottom_coarse=read_doublevector_from_raster(A_FLAG,filenames->element[I_MORPHO_ELEVATION_COARSE]+1,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
	if (elevation_bottom_coarse==NULL) printf("Error: map corresponding to %s is missing or does not have an acceptable format \n",filenames->element[I_MORPHO_ELEVATION_COARSE]+1);

	elevation_terrain_fine=read_doublevector_from_raster(A_FLAG,filenames->element[I_MORPHO_ELEVATION_TERRAINSURFACE_FINE]+1,draster->fine->layer[DTM_MASK],draster->fine->UV,dsq->fine->indices_pixel);
	if (elevation_terrain_fine==NULL) {
		elevation_terrain_fine=new_doublevector(elevation_bottom_fine->nh);
		for (i=elevation_terrain_fine->nl;i<=elevation_terrain_fine->nh;i++) {
			elevation_terrain_fine->element[i]=elevation_bottom_fine->element[i]+10000.0;
		}
		printf("Error: map corresponding to %s is missing or does not have an acceptable format (tha map is automally calculated!) \n",filenames->element[I_MORPHO_ELEVATION_TERRAINSURFACE_FINE]+1);
	} else {

		for (i=elevation_terrain_fine->nl;i<=elevation_terrain_fine->nh;i++) {
			elevation_terrain_fine->element[i]=fmax(elevation_terrain_fine->element[i],elevation_terrain_fine->element[i]);
		//	elevation_terrain_fine->element[i]=elevation_terrain_fine->element[i]+10000.0;
		}
	}
	/* creating the vector containing the lowest point of each coarse cell */
	elevation_bottom_bottom=new_doublevector(elevation_bottom_coarse->nh);
	for(i=elevation_bottom_bottom->nl;i<=elevation_bottom_bottom->nh;i++) {
		elevation_bottom_bottom->element[i]=min_elevation(i);
	}




	elevation_bottom_flines=interpolete_volume2lines(elevation_bottom_fine,dsq->fine->grid);

	elevation_terrain_flines=interpolete_volume2lines(elevation_terrain_fine,dsq->fine->grid); /* added by EC for surface flow on 22 Oct 2009 */

	/* sorting of internal small cells by Emanueke Cordano on 20100417 */
	/* start TEST ec 20100421
	long nex=2;
	printf(" %ld:",nex);
//	for (i=NL;i<=dgrid->small_polygon_content->index->element[nex];i++) {
//		printf(" %ld, ",dgrid->small_polygon_content->element[nex][i]);
//	}
//	printf(" \n \n elevation:");
	for (i=NL;i<=dgrid->small_polygon_content->index->element[nex];i++) {
	//		printf(" %lf,",elevation_bottom_fine->element[dgrid->small_polygon_content->element[nex][i]]);
			printf(" %ld, %ld ,%lf, \n",i,dgrid->small_polygon_content->element[nex][i],elevation_bottom_fine->element[dgrid->small_polygon_content->element[nex][i]]);
		}

	printf("\n");
//	print_longbin_elements(dgrid->small_polygon_content,100);
	//stop_execution();
	 END start TEST ec 20100421 */
	bubble_sort_elevation_in_longbin(dgrid->small_polygon_content,elevation_bottom_fine);
	/* start TEST ec 20100421
	//nex=100;
	printf(" %ld:",nex);
//	for (i=NL;i<=dgrid->small_polygon_content->index->element[nex];i++) {
	//	printf(" %ld,",dgrid->small_polygon_content->element[nex][i]);
//	}
//	printf(" \n \n elevation:");
	for (i=NL;i<=dgrid->small_polygon_content->index->element[nex];i++) {
	//		printf(" %lf,",elevation_bottom_fine->element[dgrid->small_polygon_content->element[nex][i]]);
		printf(" %ld, %ld ,%lf, \n",i,dgrid->small_polygon_content->element[nex][i],elevation_bottom_fine->element[dgrid->small_polygon_content->element[nex][i]]);
		}
	printf("\n");
//	print_longbin_elements(dgrid->small_polygon_content,100);
	stop_execution();
//	print_longbin_elements(dgrid->small_polygon_content,100);
	stop_execution();
	 END  TEST ec 20100421 */
	 bubble_sort_elevation_in_longbin(dgrid->small_line_content,elevation_bottom_flines);
	 /* end sorting of internal small cells by Emanueke Cordano on 20100417*/
	water_surface_elevation=read_doublevector_from_raster(A_FLAG,filenames->element[I_INITCOND_WATERSURFACE_ELEVATION]+1,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
	if (water_surface_elevation==NULL) printf("Error: map corresponding to %s is missing or does not have an acceptable format \n",filenames->element[I_INITCOND_WATERSURFACE_ELEVATION]+1);

	porosity_fine=read_doublevector_from_raster(A_FLAG,filenames->element[I_POROSITY_FINE]+1,draster->fine->layer[DTM_MASK],draster->fine->UV,dsq->fine->indices_pixel);
	if (water_surface_elevation==NULL) printf("Error: map corresponding to %s is missing or does not have an acceptable format \n",filenames->element[I_POROSITY_FINE]+1);

	water_depth=read_doublevector_from_raster(A_FLAG,filenames->element[I_INITCOND_WATERDEPTH]+1,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
	if ((water_depth==NULL) && (water_surface_elevation==NULL)) printf("Error: map corresponding to %s is missing or does not have an acceptable format \n",filenames->element[I_INITCOND_WATERDEPTH]+1);
	if ((water_depth!=NULL) && (water_surface_elevation==NULL)) { // ec 20100325
		water_surface_elevation=new_doublevector(water_depth->nh);
		printf("Warning: map corresponding to %s is missing and is approximately estimated \n",filenames->element[I_INITCOND_WATERSURFACE_ELEVATION]+1);
		for (i=water_surface_elevation->nl;i<=water_surface_elevation->nh;i++) {
			water_surface_elevation->element[i]=elevation_bottom_bottom->element[i]+water_depth->element[i];  //ec 20100325 approximate relation
		}
	}// end ec 20100325 approximate relation
	/* modified by Emanuele Cordano on 14/10/2009 */
	outlet_coefficient_v_fine=read_doublevector_from_raster(A_FLAG,filenames->element[I_OUTLET_COEFFICIENT_FINE]+1,draster->fine->layer[DTM_MASK],draster->fine->UV,dsq->fine->indices_pixel);
	if (outlet_coefficient_v_fine==NULL) {
		printf("Warning: map corresponding to %s is missing or does not have an acceptable format, there is no outlets \n",filenames->element[I_OUTLET_COEFFICIENT_FINE]+1);
		outlet_coefficient=NULL;
	} else {
		outlet_coefficient=interpolete_volume2lines(outlet_coefficient_v_fine,dgrid->fine);
	}
	/* Emanuele Cordano EC 20100421 end */
	outlet_coefficient_surf_v_fine=read_doublevector_from_raster(A_FLAG,filenames->element[I_OUTLET_COEFFICIENT_SURF_FINE]+1,draster->fine->layer[DTM_MASK],draster->fine->UV,dsq->fine->indices_pixel);
	if (outlet_coefficient_surf_v_fine==NULL) {
		printf("Warning: map corresponding to %s is missing or does not have an acceptable format, there is no outlets \n",filenames->element[I_OUTLET_COEFFICIENT_SURF_FINE]+1);
		outlet_coefficient_surf=NULL;
	} else {
		outlet_coefficient_surf=interpolete_volume2lines(outlet_coefficient_surf_v_fine,dgrid->fine);
	}


	/* modified by Emanuele Cordano on 26/10/2009 */
	runoff_coefficient_v_coarse=read_doublevector_from_raster(A_FLAG,filenames->element[I_RUNOFF_COEFFICIENT_COARSE]+1,draster->coarse->layer[DTM_MASK],draster->coarse->UV,dsq->big->indices_pixel);
	if (runoff_coefficient_v_coarse==NULL) {
			printf("Warning: map corresponding to %s is missing or does not have an acceptable format, coefficients are assumed to be zero \n",filenames->element[I_RUNOFF_COEFFICIENT_COARSE]+1);
		//	outlet_coefficient=NULL;
			runoff_coefficient=new_doublevector(dsq->big->grid->lines->nh);
			for (i=runoff_coefficient->nl;i<=runoff_coefficient->nh;i++) {
				runoff_coefficient->element[i]=0;
			}
	} else {
		runoff_coefficient=interpolete_volume2lines(runoff_coefficient_v_coarse,dsq->big->grid);
	}

	/* end modified by Emanuele Cordano on 14/10/2009 */
	water_mass_error=new_doublevector(water_surface_elevation->nh);
	/* allocation and the initialization of the surface velocity vector */
	surface_water_velocity=read_doublevector_from_a_ft_format_single_file(filenames->element[I_INITIAL_VELOCITY_FT_ARRAY]+1,print);
	/* reading subsurface from ascii file */

	if (surface_water_velocity==NULL) {
		surface_water_velocity=new_doublevector(dsq->big->grid->lines->nh);
		for(j=surface_water_velocity->nl;j<=surface_water_velocity->nh;j++) {
			surface_water_velocity->element[j]=0.0;

		}
	}
	F1_wet_vert_area=new_doublevector(dsq->big->grid->lines->nh); /*! coefficient F1 for surface flow as reported in Casulli,2008 */
	for(j=F1_wet_vert_area->nl;j<=F1_wet_vert_area->nh;j++) {
				F1_wet_vert_area->element[j]=0.0;
	}

	/* end allocation and the infiltration */
	/* read the parameters */
	read_b_parameters_ft(filenames->element[I_PARAM_FT_FILE]+1,print);

	/* time loop */

	int s=time_loop(print,write_map_results);

	free_doublevector(elevation_bottom_fine);
	free_doublevector(elevation_bottom_coarse);
	free_doublevector(elevation_terrain_fine); /* modified by Emanuele on 22/10/2009 */
	free_doublevector(elevation_bottom_flines);
	free_doublevector(elevation_terrain_flines); /* modified by Emanuele on 22/10/2009 */
	free_doublevector(elevation_bottom_bottom);
	free_doublevector(porosity_fine);
	free_doublevector(water_surface_elevation);
	free_doublevector(water_mass_error);

	free_doublevector(surface_water_velocity); /* added by Emauele Cordano on 22/10/2009  */
	free_doublevector(F1_wet_vert_area); /* added by Emauele Cordano on 29/10/2009  */

	/* modified by Emanuele on 14/10/2009 */
	if (outlet_coefficient!=NULL) free_doublevector(outlet_coefficient);
	if (outlet_coefficient_v_fine!=NULL) free_doublevector(outlet_coefficient_v_fine);
	// ec 20100421
	if (outlet_coefficient_surf!=NULL) free_doublevector(outlet_coefficient_surf);
	if (outlet_coefficient_surf_v_fine!=NULL) free_doublevector(outlet_coefficient_surf_v_fine);

	/* modified by Emanuele on 26/10/2009 */
	if (runoff_coefficient!=NULL) free_doublevector(runoff_coefficient);
	if (runoff_coefficient_v_coarse!=NULL) free_doublevector(runoff_coefficient_v_coarse);

	free_stringbin(filenames);
	free_doubleraster_map(draster);
	free_DOUBLE_GRID(dgrid); // ec 20100413
	free_doublesquare_grid(dsq);
	free(param);
	time_end=clock();
	printf("End of simulation, execution time %lf seconds  \n",(time_end-time_init)/CLOCKS_PER_SEC);
	return 0;

}


int read_b_parameters_ft(char *filename_data, short print){
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date 21 April 2009
	 *
	 * \param filename (char*) -
	 * \param print short integer
	 *
	 *\brief read the scalar data from a suitable ascii files writen in a Fluidturtle Formalism
	 */

	FILE *fd;
	DOUBLEVECTOR *scalars;
	/** block 1  - SCALAR DATA*/
	int iks=1;/*#1 ks  hydraulic saturated conductivity */
	int it_start=iks+1; /*#2 t_start - intial time */
	int it_end=it_start+1; /*#3 t_end  - end time */
	int idt=it_end+1; /*#4 dt integration time step */
	int idt_print=idt+1;/*#5 dt results printing time steps time */
	int imaxerror=idt_print+1; /*#6 maximum mass error admitted [m]*/
	int ixtemp_adm=imaxerror+1; /*#7 massimum absolute tollerance on water surface elevation [m] */
	int ideta=ixtemp_adm+1; /* #8 variation of water surce elevation utilzed for the first derivative of cell volume vs eta */
	int i_nondirichlet=ideta+1; /* #9 null value for Dirichlet nodes  */
	int i_exp_dirichlet=i_nondirichlet+1; /* #10 exponent of time in Dirichlet condition according to Lockington, et al., 2000    */
	int i_exp_outlet=i_exp_dirichlet+1; /* #11 exponent  p of the rating curve at the outlet  q_discharge=C*h_sup^p   */
    int i_exp_outlet_surf=i_exp_outlet+1; /* #12 exponent  p of the rating curve at the outlet (surface discharge) q_discharge=C_surf*h_sup^p_surf  */
	int i_exp_runoff=i_exp_outlet_surf+1;  /* #13 exponent for surface flow bottom dissipation (surface runoff) (dimensionless) */
	int i_gravity=i_exp_runoff+1; /* #14 gravity acceleration */
	int iscalars=i_gravity;
	short ifile;


	fd=t_fopen(filename_data,"r");
	ifile=read_index(fd,print);
	scalars=read_doublearray(fd,print);
	if (scalars->nh!=iscalars) {
		printf ("Warning in read_boussinesq_data_ft: Unexpected number of parameters %d instead of %ld \n ",iscalars,scalars->nh);

	}
	param=(PARAM *)malloc(sizeof(PARAM));
	if (!param) t_error("(in read_boussinesq_data_ft) param was not allocated");

//	param->Ks=scalars->element[iks];
	param->Ks=get_value_from_doublevector(iks,scalars);
	param->t_start=get_value_from_doublevector(it_start,scalars);
	//	scalars->element[it_start];
	param->t_end=get_value_from_doublevector(it_end,scalars);
	//	scalars->element[it_end];
	param->dt=get_value_from_doublevector(idt,scalars);
	//	scalars->element[idt_print];
	param->dt_print=get_value_from_doublevector(idt_print,scalars);
	// scalars->element[idt_print];
	param->max_errors=get_value_from_doublevector(imaxerror,scalars);
		// scalars->element[imaxerror];
	param->x_temp_adm=get_value_from_doublevector(ixtemp_adm,scalars);
//		scalars->element[ixtemp_adm];
	param->deta=get_value_from_doublevector(ideta,scalars);
//		scalars->element[ideta];
	param->null_dirichlet=get_value_from_doublevector(i_nondirichlet,scalars);
//		scalars->element[i_nondirichlet];
	param->exp_dirichlet=get_value_from_doublevector(i_exp_dirichlet,scalars);
//		scalars->element[i_exp_dirichlet];
	/* added on 14/10/2009 by Emanuele */
	param->p_outlet=get_value_from_doublevector(i_exp_outlet,scalars);
	// ec 20100421 i_exp_outlet_surf
	param->p_outlet_surf=get_value_from_doublevector(i_exp_outlet_surf,scalars);
	/* added on 22/10/2009 by Emanuele */
	param->p_runoff=get_value_from_doublevector(i_exp_runoff,scalars);
	/* added on 2/12/2009  */
	param->gravity=get_value_from_doublevector(i_gravity,scalars);
	if (param->gravity==-9999.0) param->gravity=9.81;

	t_fclose(fd);

	free_doublevector(scalars);
	/* Initializa time */
	param->t=param->t_start;

	return 0;



}

