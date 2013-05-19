/*!
 * \file b_solver.h
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
typedef struct {
	double dt; /*< time step integration */
	double dt_print;
	double max_errors;
	double t,t_start,t_end;
	double x_temp_adm;
	double Ks; /*!< Conductivity (homogeneous soil) */
	double deta; /*!< eta step for derivative */
	double null_dirichlet; /*< null vale for the dirichlet node */
	double exp_dirichlet;
	double p_outlet;
	double p_outlet_surf;
	double p_runoff;
	double gravity;
} PARAM;

typedef struct {
	short arithmetic_mean0;
} FLAG;

typedef struct {
	char *file_result;
	char *file_error;

//	char *SSSS;
} OUTPUT_FILENAMES;


double t_st_operator_element(long i,DOUBLEVECTOR *eta);
double t_st_operator_element_no_dirichlet(long i,DOUBLEVECTOR *eta);
//double t_st_operator_element(long i,DOUBLEVECTOR *eta,double cond_dirichlet);
double t_st_operator_element_subs(long i,DOUBLEVECTOR *eta,double cond_dirichlet);
double water_surface_elevation_mean(double eta1,double eta2);
// int T_st_operator(DOUBLEVECTOR *y, DOUBLEVECTOR *eta);
//int wet_area_operator(DOUBLEVECTOR *y, DOUBLEVECTOR *eta);

//int volume_operator(DOUBLEVECTOR *y, DOUBLEVECTOR *eta);
//int volume_operator_minus(DOUBLEVECTOR *y, DOUBLEVECTOR *eta);

double b_smatrix_element (long i,DOUBLEVECTOR *x);
int b_smatrix(DOUBLEVECTOR *y,DOUBLEVECTOR *x);

//int b_knownterm(DOUBLEVECTOR *be);

long Newton_convergence(DOUBLEVECTOR *x_temp,DOUBLEVECTOR *be, DOUBLEVECTOR *be0);
int time_loop(short print,int (*write_output)(void *v1, void *v2));
int write_map_results(void *output, void *time_string);


OUTPUT_FILENAMES *new_output_filenames(short print);
void free_OUTPUT_FILENAMES(OUTPUT_FILENAMES *fn);
