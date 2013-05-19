/*!
 * \file b_readgrid.c
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
#include "t_utilities.h"
#include "get_filenames.h"
#include "rw_maps.h"
#include "linear_span.h"
#include "geometry.h"
#include "geometry_utilities.h"
#include "read_command_line.h"
#include "additional_read_functions.h"
#include "geometry_io.h"
#include "geometry_attribute.h"
#include "geometry_freememory.h"


#include "g_raster2plvector.h"
#include "bigcells2.h"
#include "geometry2.h"
#include "keywords_file_b.h"
#include "b_solver.h"
#include "b_volumes.h"
#include "b_utilities.h"
#include "b_readgrid.h"

#define NULL_ELEVATION -9999
#define INTEGER_NULL -99
#define POLYGON_SUFFIX "_polygons.txt"
#define LINE_SUFFIX   "_lines.txt"
#define CONNECTION_SUFFIX "_connections.txt"
#define no_PRINT 0

LINEVECTOR *read_linevector (char *filename, short print) {
	/*!
	 * \author Emanuele Cordano
	 * \date May 2009
	 *
	 *\param filename (char *) - name of the file where to read line_information
	 *\param print (short)
	 *
	 *
	 *\brief It creates and reads a linevector from the following options:
	 *\brief index{1}
	 *  FILE CONTAINIG NECESSARY INFORMATION FOR LINES
			x    y    line_index    lenght2d    x_P1    y_P1    x_P2   y_P2

			1: double matrix lines information  {55550,8}
	 */

	FILE *fd;
	DOUBLEMATRIX *ldata;
	LINEVECTOR *lines;
	POINT *P1,*P2;
	long index,j,ja;

	int ix=1; /*! x coordinate of the middle point */
	int iy=ix+1; /*! y coordinate of the midddele point */
	int iline_index=iy+1; /*! iline_index of the line */
	int ilength2d=iline_index+1; /*! lenght of the line */
	int iP1_x=ilength2d+1; /*! P1 x coordinate */
	int iP1_y=iP1_x+1; /*! P1 y coordinate */
	int iP2_x=iP1_y+1; /*! P2 x coordinate */
	int iP2_y=iP2_x+1; /*! P2 y coordinate */
	int indata=iP2_y; /* number of data */

	fd=t_fopen(filename,"r");
	index=(long)read_index(fd,no_PRINT);
	ldata=read_doublematrix(fd,"a",no_PRINT);
	if (ldata->nch!=indata) printf("Warning in read_linevector: inconstancy on number of columns (data) : %ld %d \n",ldata->nch,indata);
	//if (ldata->nrh!=index)  printf("Error in read_linevector: inconstancy on number of lines : %ld %ld",ldata->nrh,index);
	lines=new_linevector(ldata->nrh);
	for (j=lines->nl;j<=lines->nh;j++) {
		ja=(long)ldata->element[j][iline_index];
		if (j!=ja) printf("Error in read_linevector (line %ld of %ld): inconstancy on line index : %ld %ld \n",j,lines->nh,j,ja);
		P1=new_point(j-1,ldata->element[j][iP1_x],ldata->element[j][iP1_y],NULL_ELEVATION);
		P2=new_point(j,ldata->element[j][iP2_x],ldata->element[j][iP2_y],NULL_ELEVATION);
		lines->element[j]=new_line_from_points(j,P1,P2);
		free_point(P1);
		free_point(P2);

	}

	free_doublematrix(ldata);
	if (print==1) printf("Function read_linevector (number of lines %ld) was successfully executed!! \n",lines->nh);

	return lines;

}

POLYGON *read_polygon(FILE *fd,short print) {
	/*!
	 * \autor Emanuele Cordano
	 * \date May 2009
	 *
	 * \param (FILE *) - file pointer
	 * \param (short) - print
	 *
	 */

	int ix=1; /* x coordinate of the centroid  */
	int iy=ix+1; /* y coordinate of the centroid */
	int ipolygon_index=iy+1; /* index of the polygon */
	int iarea2d=ipolygon_index+1; /* area of the polygon  */
	int n_data=iarea2d;
	long i;
	DOUBLEVECTOR *v_data;
	POLYGON *po;

	v_data=read_doublearray(fd,print);

	if (v_data->nh<=n_data) printf ("Error in read_polygon there no sufficient data !!\n");

	po=(POLYGON *)malloc(sizeof(POLYGON));
	if (!po) t_error("Polygon in read_polygon struct was not allocated");

	po->area2D=v_data->element[iarea2d];
	po->index=v_data->element[ipolygon_index];

	po->centroid=new_point(po->index,v_data->element[ix],v_data->element[iy],NULL_ELEVATION);

	po->edge_indices=new_longvector(v_data->nh-n_data);
	for(i=po->edge_indices->nl;i<=po->edge_indices->nh;i++) {
		po->edge_indices->element[i]=v_data->element[i+n_data];
	}

	free_doublevector(v_data);

	return po;
}

POLYGONVECTOR *read_polygonvector(char *filename,short print) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \author May 2008
	 *
	 *\param (char*) - name of filename
	 *\param (short) -
	 *
	 *\brief It creates and reads a polygonvector from the following options:
	 */
	POLYGONVECTOR *polygons;
	FILE *fd;
	long j,n_po;
	fd=t_fopen(filename,"r");

	n_po=(long)read_index(fd,no_PRINT);
	polygons=new_polygonvector(n_po);

	for (j=polygons->nl;j<=polygons->nh;j++) {
		polygons->element[j]=read_polygon(fd,print);
		if (polygons->element[j]->index!=j) printf ("Error in read_polygonvector (polygon %ld) inconstancy: %ld %ld \n",j,j,polygons->element[j]->index);
	}

	t_fclose(fd);
	if (print==1) printf("Function read_polygonvector (number of polygons %ld) was successfully executed!!",polygons->nh);

	return polygons;

}

polygon_connection_attributes *read_connections(FILE *fd,short print) {
	/*!
	 * \author Emanuele Cordano
	 * \date May 2009
	 *
	 * \param fd - (FILE *) file pointer
	 * \param print - (short)
	 *
	 *
	 */

	polygon_connection_attributes *pca;
	DOUBLEVECTOR *v_data;
	long j;

	v_data=read_doublearray(fd,no_PRINT);
	int s=(v_data->nh-1)%2;
	if (s!=0) printf("Error in read_connections (index %ld) odd number of elements in the vector after the first one which is the polygon index",v_data->nh);

	long l=(v_data->nh-1)/2;


	pca=(polygon_connection_attributes *)malloc((sizeof(polygon_connection_attributes)));
	if (!pca) printf("Error: polygon_connection_attributes was not allocated at %ld polygon",(long)v_data->element[1]);

	pca->connections=new_longvector(l);
	pca->d_connections=new_doublevector(l);

    for (j=pca->connections->nl;j<=pca->connections->nh;j++) {
    	pca->connections->element[j]=(long)(v_data->element[j*2]);
    	pca->d_connections->element[j]=v_data->element[j*2+1];

    }

    free_doublevector(v_data);

	return pca;


}


polygon_connection_attribute_array *read_connection_attribute_array(char *filename,short print) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \author May 2008
	 *
	 *\param (char*) - name of filename
	 *\param (short) -
	 *
	 *\brief It creates and reads a polygonvector from the following options:
	 */
	polygon_connection_attribute_array *pca;
	FILE *fd;
	long j,n_po;

	fd=t_fopen(filename,"r");
	n_po=(long)read_index(fd,no_PRINT);
	pca=new_connection_attributes(n_po);

	for (j=pca->nl;j<=pca->nh;j++) {
		pca->element[j]=read_connections(fd,print);

	//	if (pca->element[j]->index!=j) printf ("Error in read_connection_attributes (polygon %ld) inconstancy: %ld %ld \n",j,j,pca->element[j]->index);
	}

	t_fclose(fd);
	if (print==1) printf("Function read_connection_attributes (number of polygons %ld) was successfully executed!!",pca->nh);

	return pca;

}


GRID *read_grid(char *keyname,short print) {
	/*!
	 *\author Emanuele Cordano
	 *\date May 2008
	 *
	 */
	GRID *grid;
	long j,c;

	grid=(GRID *)malloc(sizeof(GRID));
	if (!grid) t_error("Grid in read_grid was not allocated");

	grid->file_resume_lines=copy_string(join_strings(keyname,LINE_SUFFIX));
	grid->file_resume_polygons=copy_string(join_strings(keyname,POLYGON_SUFFIX));
	grid->file_resume_connections=copy_string(join_strings(keyname,CONNECTION_SUFFIX));

	grid->lines=read_linevector(grid->file_resume_lines,print);
	grid->polygons=read_polygonvector(grid->file_resume_polygons,print);
	grid->links=read_connection_attribute_array(grid->file_resume_connections,print);

	/* check the boudaries */

	if (grid->polygons->nh!=grid->links->nh) printf ("Error in read_grid inconstancy betenn numbers of polygons %ld and of connection attributes %ld \n",grid->lines->nh,grid->links->nh);


	//j=link->nl;
	//while ((cf==0) || (df==0)) {
	for (j=grid->links->nl;j<=grid->links->nh;j++) {
		if (grid->links->element[j]->connections->nh!=grid->polygons->element[j]->edge_indices->nh)  printf ("Error in read_grid inconstancy between numbers of edge %ld and of connections  %ld at the polygon  %ld \n",grid->polygons->element[j]->edge_indices->nh,grid->links->element[j]->connections->nh,j);

		for (c=grid->links->element[j]->connections->nl;c<=grid->links->element[j]->connections->nh;c++) {
			if (grid->links->element[j]->connections->element[c]<0) grid->boundary_indicator=(long)grid->links->element[j]->connections->element[c];
		}
	}
	//grid->boundary_indicator=-10;
	if (print==1) printf("\n Function read_grid was successfully executed (number of polygons: %ld - number of lines: %ld - boundary indactor : %ld) \n",grid->polygons->nh,grid->lines->nh,grid->boundary_indicator);


	//}



//	char *file_resume_lines_fine=join_strings(resume_filenames,"_lines_fine.txt");
//	char *file_resume_polygons_fine=join_strings(resume_filenames,"_polygons_fine.txt");
//	char *file_resume_connections_fine=join_strings(resume_filenames,"_connections_fine.txt");
//	char *file_resume_c_polygon=join_strings(resume_filenames,"_c_polygon.txt");
//	char *file_resume_c_line=join_strings(resume_filenames,"_c_line.txt");

	return grid;
}

SQUARE_GRID *read_square_grid(char *keyname,DOUBLEMATRIX *DTM,long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index),DOUBLEVECTOR *V,int (*check_novalues)(double x, DOUBLEVECTOR *V),short print) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date May 2009
	 *
	 *\param keyname (char *) - keyname of the SQUARE_GRID
	 *\param DTM     (DOUBLEMATRIX *) - Digital Terrain Model
	 *\param long (*index_pixel_from_a_bin)(long r, long c,LONGVECTOR *s_index) - function for the space filling curve
	 *\param V - (DOUBLEVECTOR *) - vector containig novalue information
	 *\param int (*check_novalues)(double x, DOUBLEVECTOR *V) - function identintifyng novalue information
	 *\param print (short)
	 *
	 *\return a SQUARE GRID whose data are written in files with FLUIdTurtle formalism
	 *
	 *
	 *
	 */

	SQUARE_GRID *sq;
	long count,r,c;

//	sq=(SQUARE_GRID *)malloc(sizeof(SQUARE_GRID));
//	if (!sq) t_error("Square Grid sq was not allocated");

//	sq->indices_pixel=m_indices_from_mask(DTM,0,0,INTEGER_NULL,0,index_pixel_from_a_bin,V,(*check_novalues));
//	sq->indices_vertex=m_indices_from_mask(DTM,-1,-1,INTEGER_NULL,0,(*index_pixel_from_a_bin),V,(*check_novalues));
//	sq->indices_horizontal_lines=m_indices_from_mask(DTM,-1,0,INTEGER_NULL,0,(*index_pixel_from_a_bin),V,(*check_novalues));

	sq=(SQUARE_GRID *)malloc(sizeof(SQUARE_GRID));
	if (!sq) t_error("Square Grid sq was not allocated");
	sq->indices_pixel=m_indices_from_mask(DTM,0,0,INTEGER_NULL,0,index_pixel_from_a_bin,V,(*check_novalues));
	sq->indices_vertex=m_indices_from_mask(DTM,-1,-1,INTEGER_NULL,0,(*index_pixel_from_a_bin),V,(*check_novalues));
	sq->indices_horizontal_lines=m_indices_from_mask(DTM,-1,0,INTEGER_NULL,0,(*index_pixel_from_a_bin),V,(*check_novalues));
	count=0;
	for(r=sq->indices_horizontal_lines->nrl;r<=sq->indices_horizontal_lines->nrh;r++){
		for (c=sq->indices_horizontal_lines->ncl;c<=sq->indices_horizontal_lines->nch;c++){
	    	if (sq->indices_horizontal_lines->element[r][c]!=INTEGER_NULL) count++;
		}
	}
	sq->nhorizontal_lines=count;
	sq->novalue=INTEGER_NULL;
	sq->indices_vertical_lines=m_indices_from_mask(DTM,0,-1,INTEGER_NULL,count,(*index_pixel_from_a_bin),V,(*check_novalues));


    printf("Warning: index matrix of SQUARE_GRID (already created)  %s was not checked!!\n",keyname);
    sq->novalue=sq->indices_pixel->element[1][1]; /*! first element of sqp->inices_pixel is assumed as a no-value */
	sq->grid=read_grid(keyname,print);
    if (print==1) printf("Function read_square_grid (%s) was correctly executed!!! \n",keyname);
	return sq;

}

LONGMATRIX_VECTOR *read_fine_indices(char *filename,short print) {
	/*!
	 * \author Emanuele Cordano
	 * date May 2009
	 *
	 * \param (char*) - name of filename
	 *\param (short) -
	 */
	LONGMATRIX_VECTOR *lm;
	FILE *fd;
	long j,n_po;
	fd=t_fopen(filename,"r");

	n_po=(long)read_index(fd,no_PRINT);
	lm=new_longmatrix_vector(n_po);

	for (j=lm->nl;j<=lm->nh;j++) {
		if (print==1) printf("Function read_fine_indices element %ld of %ld is being defined \n",j,lm->nh);
		lm->element[j]=read_longmatrix(fd,"a",no_PRINT);
		if (print==1) printf("Function read_fine_indices element %ld of %ld is read [%ld,%ld] f.e. %ld \n",j,lm->nh,lm->element[j]->nrh,lm->element[j]->nrh,lm->element[j]->element[1][1]);
		if (!lm->element[j]) t_error("lm->element[j] in read_fine_indices was not defined");

	//		if (polygons->element[j]->index!=j) printf ("Error in read_polygonvector (polygon %ld) inconstancy: %ld %ld \n",j,j,polygons->element[j]->index);
	}

	t_fclose(fd);

	if (print==1) printf("Function read_fine_indices was successfully executed! \n");

	return lm;
}


LONGBIN *read_line_indices(char *filename,short print) {
	/*
	 *
	 * \author Emanuele Cordano
	 * date May 2009
	 *
	 * \param (char*) - name of filename
	 *\param (short) -
	 *
	 */
	LONGVECTOR **lv;
	LONGBIN *lb;
	LONGVECTOR *vi;
	FILE *fd;
	long j,n_l,c;
	fd=t_fopen(filename,"r");

	n_l=(long)read_index(fd,print);
//	printf("n_l=%ld",n_l);

	lv=(LONGVECTOR **)malloc((size_t)(n_l*sizeof(LONGVECTOR *)));
//	stop_execution();
	for (j=0;j<=n_l-1;j++) {
		lv[j]=read_longarray(fd,print);
	}
	vi=new_longvector(n_l);
	for (j=vi->nl;j<=vi->nh;j++) {
		vi->element[j]=lv[j-1]->nh;
	}
	lb=new_longbin(vi);
	for (j=lb->index->nl;j<=lb->index->nh;j++) {
		for(c=1;c<=lb->index->element[j];c++) {
			lb->element[j][c]=lv[j-1]->element[c];
		}
	}

	free_longvector(vi);
	for (j=0;j<=n_l-1;j++) {
			free_longvector(lv[j]);
	}
	free(lv);
	t_fclose(fd);

	if (print==1) printf("Function read_line_indices was successfully executed! \n");




	return lb;


}



DOUBLESQUARE_GRID *read_doublesquare_grid (DOUBLERASTER_MAP *draster, char *keyname,long (*index_pixel_from_a_bin_coarse)(long r, long c,LONGVECTOR *s_index),long (*index_pixel_from_a_bin_fine)(long r, long c,LONGVECTOR *s_index),short print){
	/*
	 *
	 * \author Emanuele Cordano
	 * \date May 2009
	 *
	 *\param - DOUBLERASTER_MAP *draster
	 *\param - (char *) root neme of textfiles containing the struct information.
	 *\param long (*index_pixel_from_a_bin_coarse)(long r, long c,LONGVECTOR *s_index) - equation of the filling curve for the pixel of a coarse grid
	 *\param long (*index_pixel_from_a_bin_fine)(long r, long c,LONGVECTOR *s_index),long d_coarse,long d_fine,short print) - equation of the filling curve for the pixel of a fine grid

	 *\param (short) print
	 *
	 */

	DOUBLESQUARE_GRID *dsq;

	char *keyname_coarse=join_strings(keyname,"__coarse");
	char *keyname_fine=join_strings(keyname,"__fine");

	char *filename_c_polygon=join_strings(keyname,"_c_polygon.txt");
	char *filename_c_line=join_strings(keyname,"_c_line.txt");

	dsq=(DOUBLESQUARE_GRID *)malloc(sizeof(DOUBLESQUARE_GRID));
	if (!dsq) t_error("Double Square Grid dsq in read_doublesquare_grid was not allocated");

	dsq->big=read_square_grid(keyname_coarse,draster->coarse->layer[draster->coarse->reference_index_map],index_pixel_from_a_bin_coarse,draster->coarse->UV->V,draster->coarse->check_novalues,print);
	dsq->fine=read_square_grid(keyname_fine,draster->fine->layer[draster->fine->reference_index_map],index_pixel_from_a_bin_fine,draster->fine->UV->V,draster->fine->check_novalues,print);

	dsq->file_resume_c_line=copy_string(filename_c_line);
	dsq->file_resume_c_polygon=copy_string(filename_c_polygon);

	dsq->small_content_line=read_line_indices(dsq->file_resume_c_line,print);
	dsq->small_content_polygon=read_fine_indices(dsq->file_resume_c_polygon,print);

	if (dsq->small_content_line->index->nh!=dsq->big->grid->lines->nh) {
		printf("Error in read_doublesquare_grid inconstancy with number of (big) lines %ld and %ld \n",dsq->small_content_line->index->nh,dsq->big->grid->lines->nh);
		stop_execution();
	}
	if (dsq->small_content_polygon->nh!=dsq->big->grid->polygons->nh) {
		printf("Error in read_doublesquare_grid inconstancy with number of (big) lines %ld and %ld \n",dsq->small_content_polygon->nh,dsq->big->grid->polygons->nh);
		stop_execution();
	}

	if (print==1) printf("Function read_doublesquare_grid was successfully executed! \n ");

	return dsq;
}



