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

double F1_coefficient(long j, double velocity, double wet_vert_area);
double b_advection(long i);
double t_st_advection_operator_element(long i, DOUBLEVECTOR *eta,double cond_dirichlet);
double symmetric_surface_velocity(long j, double forcing);
double asymmetric_surface_velocity(long j, double eta_previous);
int update_velocity(DOUBLEVECTOR *eta);
int update_F1_wet_vert_area();
