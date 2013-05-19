/*!
 * \file b_volumes.h
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

double volume (double eta, long i);
double volume_subs(double eta, long i);
double volume_surf(double eta, long i);

double vertical_area(double eta, long j);
double vertical_area_subs(double eta, long j);
double vertical_area_surf(double eta, long j);
double wet_area(double eta,long i, double deta);

double min_elevation(long i);

//added by EC on 14_10_2009
double q_discharge_from_outlet_subs_line(double eta, long j);
double q_discharge_from_outlet_cell(double eta, long i);
