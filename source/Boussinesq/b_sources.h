/*!
 * \file b_sources.h
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
	DOUBLEVECTOR *times;
	STRINGBIN *s_suffixes;
} S_TIMES;


S_TIMES  *get_s_times(char *filename,short print);
void free_s_times(S_TIMES* s_t);

int get_sources(double t,DOUBLEVECTOR *sources);
int get_dirichletsnode(double t,DOUBLEVECTOR *dirichlet);
