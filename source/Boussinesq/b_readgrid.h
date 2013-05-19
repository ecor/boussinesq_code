/*!
 * \file b_readgrid.h
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
LINEVECTOR *read_linevector (char *filename, short print);
POLYGON *read_polygon(FILE *fd,short print);
POLYGONVECTOR *read_polygonvector(char *filename,short print);
polygon_connection_attributes *read_connections(FILE *fd,short print);
polygon_connection_attribute_array *read_connection_attribute_array(char *filename,short print);

GRID *read_grid(char *keyname,short print);

SQUARE_GRID *read_square_grid(char *keyname,DOUBLEMATRIX *DTM,long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index),DOUBLEVECTOR *V,int (*check_novalues)(double x, DOUBLEVECTOR *V),short print);
LONGMATRIX_VECTOR *read_fine_indices(char *filename,short print);
LONGBIN *read_line_indices(char *filename,short print);
DOUBLESQUARE_GRID *read_doublesquare_grid (DOUBLERASTER_MAP *draster, char *keyname,long (*index_pixel_from_a_bin_coarse)(long r, long c,LONGVECTOR *s_index),long (*index_pixel_from_a_bin_fine)(long r, long c,LONGVECTOR *s_index),short print);
