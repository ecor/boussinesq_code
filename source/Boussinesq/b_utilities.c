/*!
 * \file b_utilities.c
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
#include "rw_maps.h"
#include "geometry.h"
#include "g_raster2plvector.h"
#include "bigcells2.h"
#include "geometry2.h"
#include "b_utilities.h"

#define D_INIT_VALUE -9999.0
//#include "rw_maps.h"
//#include "gridded.element.input.geotop.h"

//#include "g_raster2plvector.h"
DOUBLEMATRIX *doublemap_from_longmap(LONGMATRIX *lmap, long novalue, double dnovalue){
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 * \param lmap - (LONGMATRIX *) - mask
	 * \param novalue - (long) null value (long)
	 * \param dnovalue - (double) null value (double)
	 *
	 * \return DOUBLEMATRIX
	 */

	long r,c;
	DOUBLEMATRIX *M;
	M=new_doublematrix(lmap->nrh,lmap->nch);

	for (r=lmap->nrl;r<=lmap->nrh;r++){
		for (c=lmap->ncl;c<=lmap->nch;c++)
			if (lmap->element[r][c]==novalue) {
				M->element[r][c]=dnovalue;
			} else {
				M->element[r][c]=(double)lmap->element[r][c];
		}
	}

	return M;

}

DOUBLEVECTOR *interpolete_volume2lines(DOUBLEVECTOR *v,GRID *grid) {
	/*!!
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 * \param v  - (DOUBLEVECTOR *) - vector of the quantiatyies related to the polygons
	 *
	 *
	 * \return a vector of quantity v referred and intarpolated in the lines between two polygons
	 */
	DOUBLEVECTOR *ve;
	long j,i,l,kl,kp,kbond;
	kbond=grid->boundary_indicator;

	if (grid->polygons->nh!=v->nh) printf ("Error in interpolete_volume2lines v (quantity vector) size (%ld) is not equal to the numbers of polygons (%ld) !! \n",v->nh,grid->polygons->nh);

	ve=new_doublevector(grid->lines->nh);

	for (j=ve->nl;j<=ve->nh;j++) {
		ve->element[j]=D_INIT_VALUE;
	}

	for (i=v->nl;i<=v->nh;i++) {
		for(l=grid->polygons->element[i]->edge_indices->nl;l<=grid->polygons->element[i]->edge_indices->nh;l++){
			kl=grid->polygons->element[i]->edge_indices->element[l];
			kp=grid->links->element[i]->connections->element[l];
			if (kp==kbond) {
				if (ve->element[kl]==D_INIT_VALUE) ve->element[kl]=v->element[i];
			} else {
				if (ve->element[kl]==D_INIT_VALUE) ve->element[kl]=(v->element[i]+v->element[kp])/2.0;
			}
		}
	}
	for (j=ve->nl;j<=ve->nh;j++) {
		if (ve->element[j]==D_INIT_VALUE) printf ("Error in interpolete_volume2lines v (quantity vector)  (%ld) at position %ld was not well initialized!! \n",ve->nh,j);
	}

	return ve;
}



DOUBLEMATRIX *get_doublemap_spav(LONGMATRIX *lmap, long novalue, double dnovalue,long dr, long dc, DOUBLEMATRIX *Msource){
	/*!
	 * \author Emanuele Cordano
	 * \date 17 April 2009
	 *
	 * \param lmap - (LONGMATRIX *) - mask
	 * \param novalue - (long) null value (long)
	 * \param dnovalue - (double) null value (double)
	 * \param dr - (long) dr step at row with which function get_doublemap_spav calculates the average (-1 horizontal, 0 vertical)
	 * \param dc - (long) dc step at  column with which function get_doublemap_spav calculates the average (0 horizontal, -1 vertical)
	 * \param Msource - (DOUBLEMATRIX *) Msources
	 *
	 * \return a map whose values are taken from Msource and whose mask is taken from lmap
	 */

	long r,c;
	DOUBLEMATRIX *M;
	M=new_doublematrix(lmap->nrh,lmap->nch);

	for (r=lmap->nrl;r<=lmap->nrh;r++){
			for (c=lmap->ncl;c<=lmap->nch;c++) {
				if ((lmap->element[r][c]==novalue) || (r==lmap->nrl) || (r==lmap->nrh) || (c==lmap->ncl) || (c==lmap->nrh)) {
					M->element[r][c]=dnovalue;
				} else if ((Msource->element[r][c]!=dnovalue) && (Msource->element[r+dr][c+dc])!=dnovalue){
					M->element[r][c]=(Msource->element[r][c]+Msource->element[r+dr][c+dc])/2.0;

				} else if ((Msource->element[r][c]==dnovalue) && (Msource->element[r+dr][c+dc]!=dnovalue)){
					M->element[r][c]=Msource->element[r+dr][c+dc];

				} else  if ((Msource->element[r][c]!=dnovalue) && (Msource->element[r+dr][c+dc]==dnovalue)) {
					M->element[r][c]=Msource->element[r][c];

				} else  {
					printf("Error 1 in get_doublemap_spav no value : address[%ld][%ld]=%ld and ref_map[%ld][%ld]=%lf ref_map[%ld][%ld]=%lf \n",r,c,lmap->element[r][c],r,c,Msource->element[r][c],r+dr,c+dc,Msource->element[r+dr][c+dc]);
					return NULL;
				}
			}
	}

	for (r=M->nrl;r<=M->nrh;r++){
		for (c=M->ncl;c<=M->nch;c++){
			if ((lmap->element[r][c]==novalue) && (M->element[r][c]==dnovalue)) {
			} else if ((lmap->element[r][c]==novalue) || (M->element[r][c]==dnovalue)) {
				printf("Error 2 in get_doublemap_spav no value : (row=%ld col=%ld) address[%ld][%ld]=%ld and ref_map[%ld][%ld]=%lf ref_map[%ld][%ld]=%lf map[%ld][%ld]=%lf\n",M->nrh,M->nch,r,c,lmap->element[r][c],r,c,Msource->element[r][c],r+dr,c+dc,Msource->element[r+dr][c+dc],r,c,M->element[r][c]);
				return NULL;
			}
		}
	}

	return M;

}


DOUBLEVECTOR *get_doublevector_for_lines(LONGMATRIX *h_addresses,LONGMATRIX *v_addresses,DOUBLEMATRIX *mh,DOUBLEMATRIX *mv,double novalue) {
	/*!
	 *\author Emanuele Cordano
	 *\date 17 April 2009
	 *
	 *No\param n_horizontal (long) - number of horizontal lines
	 *No \param n_lins (long) - number of lines
	 *\param h_addresses - (LONGMATRIX *) adresses of horizotal lines;
	 *\param v_addresses - (LONGMATRIX *) adresses of vertical  lines;
	 *\param mh - (DOUBLEMATRIX *) matrix for horizontal lines;
	 *\param mv - (DOUBLEMATRIX *) matrix for vertical lines;
	 *\param novalue - (double) novalue
	 */

	//DOUBLEVECTOR *get_doublevector_from_doublematrix(LONGMATRIX *indices,DOUBLEMATRIX *M, double novalue){
	long r,c,j;
	long n_horizontal;

	DOUBLEVECTOR *v_horizontal,*v_vertical, *v;
	printf ("0) function get_doublevector_for_lines  (novalue=%lf)\n",novalue);
	v_horizontal=get_doublevector_from_doublematrix(h_addresses,mh,novalue);
	n_horizontal=v_horizontal->nh;
	printf ("1) function get_doublevector_for_lines  n_horizontal=%ld \n",n_horizontal);
	for (r=v_addresses->nrl;r<=v_addresses->nrh;r++) {
		for(c=v_addresses->ncl;c<=v_addresses->nch;c++) {
			if (mh->element[r][c]!=novalue)  v_addresses->element[r][c]=v_addresses->element[r][c]-n_horizontal;
		}
	}


	v_vertical=get_doublevector_from_doublematrix(v_addresses,mv,novalue);
	v=new_doublevector(v_horizontal->nh+v_vertical->nh);
	for (j=v->nl;j<=n_horizontal;j++) {
		v->element[j]=v_horizontal->element[j];
	}
	for (j=n_horizontal+1;j<=v->nh;j++) {
		v->element[j]=v_vertical->element[j-n_horizontal];
	}

	for (r=v_addresses->nrl;r<=v_addresses->nrh;r++) {
		for(c=v_addresses->ncl;r<=v_addresses->nch;c++) {
			if (mh->element[r][c]!=novalue) v_addresses->element[r][c]=v_addresses->element[r][c]+n_horizontal;
		}

	}

	free_doublevector(v_horizontal);
	free_doublevector(v_vertical);
	printf ("function get_doublevector_for_lines n_lines=%ld n_horizontal=%ld \n",v->nh,n_horizontal);
	return v;
}








int write_raster_from_doublevector_v2(char *filename, DOUBLEVECTOR *v, long n0, T_INIT *UVref, LONGMATRIX *indices, DOUBLEMATRIX *Mref, short float_type, short map_format) {
	/*!
	 * \author Emanuele Cordano
	 *
	 * \date March 2009
	 *

	 * \param filename (char *) name of the file
	 * \param v - (DOUBLECTOR *) - vector to be mapped
	 * \param n0 - (long) starting point of v (i-th element of v is put in the pixels whose index (in indices matrix) is i+n0)
	 * \param UVref - (T_INIT *) T_IIT struct containg ewres and nwres information
	 * \param indices  (LONGVECTOR *) - matrix of indices
     * \param Mref - (DOUBLEMATRIX *) - reference raster
     * \param float_type - float type for printed values (floating = 0, recommended)
     * \param map_format - map asccii format (FLUIDTURLE, GRASSASCI or ESRIASCII) (GRASSASCI=2 recommended)

	 *

	 * \brief it writes a map from a vector using get_doublematrix_from_doublevecto and write_map
	 */

	DOUBLEMATRIX *M;
	long r,c;

	for(r=indices->nrl;r<=indices->nrh;r++){
		for(c=indices->ncl;c<=indices->nch;c++){
			indices->element[r][c]=indices->element[r][c]-n0;
		}
	}
	M=get_doublematrix_from_doublevector(v,indices,Mref,UVref->V->element[2]);
	write_map(filename,float_type,map_format,M,UVref);

	for(r=indices->nrl;r<=indices->nrh;r++){
		for(c=indices->ncl;c<=indices->nch;c++){
			indices->element[r][c]=indices->element[r][c]+n0;
		}
	}
	free_doublematrix(M);

	return 0;

}

int write_line_map(char *filename, DOUBLEVECTOR *v, long n_horizontal, T_INIT *UVref, LONGMATRIX *indices_h, LONGMATRIX *indices_v,DOUBLEMATRIX *Mref_h,DOUBLEMATRIX *Mref_v, short float_type, short map_format) {
/*!
  * \author Emanuele Cordano
	 *
	 * \date March 2009
	 *

	 * \param filename (char *) name of the file
	 * \param v - (DOUBLECTOR *) - vector to be mapped
	 * \param nhorizontal - (long) - number of horizontal lines
	 * \param UVref - (T_INIT *) T_IIT struct containg ewres and nwres information
	 * \param indices_h  (LONGVECTOR *) - matrix of indices for horizontal lines
	 * \param indices_v  (LONGVECTOR *) - matrix of indices for vertical lines
     * \param Mref_h - (DOUBLEMATRIX *) - reference raster for horizotal lines
     * \param Mref_v - (DOUBLEMATRIX *) - reference raster for vertical lines
     * \param float_type - float type for printed values (floating = 0, recommended)
     * \param map_format - map asccii format (FLUIDTURLE, GRASSASCI or ESRIASCII) (GRASSASCI=2 recommended)

	 *
 *
 */

	char *filename_h=join_strings(filename,"_horizontal_lines");
	char *filename_v=join_strings(filename,"_vertical_lines");
	DOUBLEVECTOR *v_vertical, *v_horizontal;

	long i;
	v_horizontal=new_doublevector(n_horizontal);
	v_vertical=new_doublevector(v->nh-n_horizontal);
	for (i=v_horizontal->nl;i<=v_horizontal->nh;i++) {
		v_horizontal->element[i]=v->element[i];
	}

	for (i=v_vertical->nl;i<=v_vertical->nh;i++) {
		v_vertical->element[i]=v->element[i+v_horizontal->nh]; /* modified by Emanuele Cordano replace v_vertical->nh with v_horizontal->nh with 9 September 2009 */
	}


	write_raster_from_doublevector_v2(filename_h,v_horizontal,0,UVref,indices_h,Mref_h,float_type,map_format);
	write_raster_from_doublevector_v2(filename_v,v_vertical,n_horizontal,UVref,indices_v,Mref_v,float_type,map_format);

	free_doublevector(v_horizontal);
	free_doublevector(v_vertical);

	//int write_line_map(char *filename, DOUBLEVECTOR *v, long n_horizontal, T_INIT *UVref, LONGMATRIX *indices_h, LONGMATRIX *indices_v,DOUBLEMATRIX *Mref_h,DOUBLEMATRIX *Mref_v,float_type,map_format);
	free(filename_h);
	free(filename_v);
	return 0;


}

int write_lines(char *filename,long n_lines, long n_horizontal, T_INIT *UVref, LONGMATRIX *indices_h, LONGMATRIX *indices_v,DOUBLEMATRIX *Mref_h,DOUBLEMATRIX *Mref_v, short float_type, short map_format) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date April 2009
	 *
	 * \param filename (char *) name of the file
	 * \param n_lines - (long) - number of lines
	 * \param nhorizontal - (long) - number of horizontal lines
	 * \param UVref - (T_INIT *) T_IIT struct containg ewres and nwres information
	 * \param indices_h  (LONGVECTOR *) - matrix of indices for horizontal lines
	 * \param indices_v  (LONGVECTOR *) - matrix of indices for vertical lines
     * \param Mref_h - (DOUBLEMATRIX *) - reference raster for horizotal lines
     * \param Mref_v - (DOUBLEMATRIX *) - reference raster for vertical lines
     * \param float_type - float type for printed values (floating = 0, recommended)
     * \param map_format - map asccii format (FLUIDTURLE, GRASSASCI or ESRIASCII) (GRASSASCI=2 recommended)

	 *
	 * \brief write the addresse map for lines
	 *
	 */

	long i;
	DOUBLEVECTOR *v;
	v=new_doublevector(n_lines);
	for (i=v->nl;i<=v->nh;i++) {
		v->element[i]=(double)i;
	}
	write_line_map(filename,v,n_horizontal,UVref,indices_h,indices_v,Mref_h,Mref_v,float_type,map_format);
	free_doublevector(v);
	return 0;
}

int write_cell(char *filename, long n_cells, T_INIT *UVref, LONGMATRIX *indices, DOUBLEMATRIX *Mref, short float_type, short map_format) {
	/*!
	 * \author Emanuele Cordano
	 * \date April 2008
	 *
	 * \param filename (char *) name of the file
	 * \param ncells - (long) number of cells
	 * \param UVref - (T_INIT *) T_IIT struct containg ewres and nwres information
	 * \param indices  (LONGVECTOR *) - matrix of indices
     * \param Mref - (DOUBLEMATRIX *) - reference raster
     * \param float_type - float type for printed values (floating = 0, recommended)
     * \param map_format - map asccii format (FLUIDTURLE, GRASSASCI or ESRIASCII) (GRASSASCI=2 recommended)
	 *
	 *
	 */
	long i;
	DOUBLEVECTOR *v;
	v=new_doublevector(n_cells);
	for (i=v->nl;i<=v->nh;i++){
		v->element[i]=(double)i;
	}
	write_raster_from_doublevector_v2(filename,v,0,UVref,indices,Mref,float_type,map_format);
	free_doublevector(v);
	return 0;
}

int zeros(DOUBLEVECTOR *v) {
/*!
 * \author Emanuele Cordanp
 *
 * \date 19 April 2009
 *
 */
long i;
	for(i=v->nl;i<=v->nh;i++){
		v->element[i]=0.0;
	}
return 0;
}

int check_matrices(DOUBLEMATRIX *M1,LONGMATRIX *L1, double dnull, long lnull,short print) {
/*!
 * \author Emanule Cordano
 * \date 26 April 2009
 *
 */
	long r,c;
	int s;
	if ((M1->nrh!=L1->nrh) || (M1->nch!=L1->nch)) printf("Error:: in check_matrices [%ld,%ld] and M [%ld,%ld] has different sizes! \n",L1->nrh,L1->nch,M1->nrh,M1->nch);
	s=0;
	for (r=M1->nrl;r<=M1->nrh;r++) {
		for (c=M1->ncl;c<=M1->nch;c++) {
			if ((M1->element[r][c]==dnull) && (L1->element[r][c]==lnull)) {
			} else if ((M1->element[r][c]==dnull) || (L1->element[r][c]==lnull)) {
				printf ("Error Function check_matrices: dmatrix and lomatrix [%ld,%ld] have different values %lf and %ld!! \n",r,c,M1->element[r][c],L1->element[r][c]);
				printf("Neighnourg ponts: dmatrix[%ld][%ld]=%lf and dmatrix[%ld][%ld]=%lf\n",r-1,c,M1->element[r-1][c],r,c-1,M1->element[r][c-1]);
				return -1;
			}
		}
	}
	if (print==1) printf("Function check matrix Exit %d\n",s);

	return s;
}

double get_value_from_doublevector(long k,DOUBLEVECTOR *v) {
	/*!
	 *
	 *
	 *\author Emanuele Cordano
	 *\date 14 October 2009
	 *
	 *\param k (long) - element requested
	 *\param v (DOUBLECTOR *)
	 *
	 *\return the k-th element of the vector v, an error value (-9999) in case k exceeds the size of v
	 */
	double error_v=-9999.0;
	if ((k<v->nl) || (k>v->nh)) {
		printf("Warning in function in get_value_from_doublevector: value %ld execeeds the size of the vector %ld %ld, function returns %lf ! \n",k,v->nl,v->nh,error_v);
		return error_v;
	} else {
		return v->element[k];
	}
	printf("Error in function in get_value_from_doublevector: value %ld (between %ld %ld) is not read correctly, function returns %lf ! \n",k,v->nl,v->nh,error_v);
	return error_v;
}

int write_doublevector_in_a_ascii_ft_format_single_file(char *filename, DOUBLEVECTOR *v) {
	/*
	 * \author Emanuele Cordano
	 * \date 18 March 2010
	 *
	 *
	 */
	char *function_name="write_doublevector_in_ft_format";
	FILE *fd;
	long l;
//	if (!v->name) v->name="missing_name";
	fd=t_fopen(filename,"w");
	fprintf(fd,"index{1} \n");
	fprintf(fd,"\n");
	fprintf(fd,"\n");
	fprintf(fd,"1: double array STANDARD_NAME  {\n");
	for (l=v->nl;l<v->nh;l++) {
		fprintf(fd,"%le, \n",v->element[l]);
	}
	fprintf(fd,"%le} \n",v->element[v->nh]);

	t_fclose(fd);

	return 0;

}

DOUBLEVECTOR *read_doublevector_from_a_ft_format_single_file (char *filename,short print) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date 18 March 2010
	 *
	 *
	 */
	DOUBLEVECTOR *vect=NULL;
	FILE *fd;
	char *function_name="read_doublevector_from_a_ft_format_single_file";
	fd=fopen(filename,"r");

	if (fd==NULL) {
		printf("Waring in %s: velocity file is missing",function_name);

		return NULL;
	}
	return NULL;
	short index=read_index(fd,print);
	vect=read_doublearray(fd,print);
	fclose(fd);
	return vect;
}
