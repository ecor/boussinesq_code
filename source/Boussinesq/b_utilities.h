/*!
 * \file b_utilities.h
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



DOUBLEMATRIX *doublemap_from_longmap(LONGMATRIX *lmap, long novalue, double dnovalue);
DOUBLEVECTOR *interpolete_volume2lines(DOUBLEVECTOR *v,GRID *grid);
DOUBLEMATRIX *get_doublemap_spav(LONGMATRIX *lmap, long novalue, double dnovalue,long dr, long dc, DOUBLEMATRIX *Msource);
DOUBLEVECTOR *get_doublevector_for_lines(LONGMATRIX *h_addresses,LONGMATRIX *v_addresses,DOUBLEMATRIX *mh,DOUBLEMATRIX *mv,double novalue);
int write_raster_from_doublevector_v2(char *filename, DOUBLEVECTOR *v, long n0, T_INIT *UVref, LONGMATRIX *indices, DOUBLEMATRIX *Mref, short float_type, short map_format);

int write_line_map(char *filename, DOUBLEVECTOR *v, long n_horizontal, T_INIT *UVref, LONGMATRIX *indices_h, LONGMATRIX *indices_v,DOUBLEMATRIX *Mref_h,DOUBLEMATRIX *Mref_v, short float_type, short map_format);

int write_lines(char *filename,long n_lines, long n_horizontal, T_INIT *UVref, LONGMATRIX *indices_h, LONGMATRIX *indices_v,DOUBLEMATRIX *Mref_h,DOUBLEMATRIX *Mref_v, short float_type, short map_format);

int write_cell(char *filename, long n_cells, T_INIT *UVref, LONGMATRIX *indices, DOUBLEMATRIX *Mref, short float_type, short map_format);
int zeros(DOUBLEVECTOR *v);

int check_matrices(DOUBLEMATRIX *M1,LONGMATRIX *L1, double dnull, long lnull,short print);
double get_value_from_doublevector(long k,DOUBLEVECTOR *v);


int write_doublevector_in_a_ascii_ft_format_single_file(char *fileneme, DOUBLEVECTOR *v);
DOUBLEVECTOR *read_doublevector_from_a_ft_format_single_file (char *filename,short print);
