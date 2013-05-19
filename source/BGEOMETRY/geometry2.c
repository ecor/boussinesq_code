
#include "turtle.h"
#include "rw_maps.h"
#include "geometry.h"
#include "geometry_utilities.h"
#include "geometry_attribute.h"
#include "geometry_freememory.h"
#include "geometry_io.h"
#include "g_raster2plvector.h"
#include "sorting.h"
#include "linear_span.h"
#include "bigcells2.h"
#include "geometry2.h"

POINT *new_point_from_point(POINT *point) {
	/*!
	 *
	 * \param (POINT *) point
	 *
	 * \author Emanuele Cordano
	 * \date April 2010
	 *
	 */
	POINT *p_new;

	p_new=new_point(point->index,point->x,point->y,point->z);

	return p_new;
}

LINE *new_line_from_line(LINE *line) {
	/*!
	 *
	 * \param (LINE *) line
	 * \author Emanuele Cordano
	 *
	 * \date April 2010
	 *
	 *
	 */
	LINE *l_new;
//	POINT *p_start,*p_end;


	l_new=new_line_from_points(line->index,line->begin,line->end);
//	free_line(l_new);
//	free_line(line);
//	printf("freed lines");
//	stop_execution();

	return l_new;
}


LINEVECTOR *new_linevector_from_linevector(LINEVECTOR *lines) {
	/*!
	 *
	 * \param (LINEVECTOR *) line
	 *
	 * \author Emanuele Cordano
	 *
	 * \date April 2010
	 */
	LINEVECTOR *lines_new;
	long j;

	lines_new=new_linevector(lines->nh);

	for (j=lines->nl;j<=lines->nh;j++){
		lines_new->element[j]=new_line_from_line(lines->element[j]);

	}

	return lines_new;

}

POLYGON *new_polygon_from_polygon(POLYGON *polygon) {
	/*!
	 * \param (POLYGON *) polygon;
	 *
	 * \author Emanuele Cordano
	 *
	 */
	POLYGON *pol_new;
	char *function_name="new_polygon_from_polygon";
	long j;

	pol_new=(POLYGON *)malloc(sizeof(POLYGON));
	if (!pol_new) printf("Error in %s , polygon was not allocated \n",function_name);

	pol_new->index=polygon->index;
	pol_new->area2D=polygon->area2D;

	pol_new->centroid=new_point_from_point(polygon->centroid);
	//linee

	pol_new->edge_indices=new_longvector(polygon->edge_indices->nh);
	for (j=pol_new->edge_indices->nl;j<=pol_new->edge_indices->nh;j++) {
		pol_new->edge_indices->element[j]=polygon->edge_indices->element[j];
	}

	return pol_new;

}

POLYGONVECTOR *new_polygonvector_from_polygonvector(POLYGONVECTOR *polygons) {
	/*!
	 *
	 * \param (POLYGONVECTOR *) - polygons
	 *
	 * \author Emanuele Cordano
	 * \date April 2010
	 *
	 */

	POLYGONVECTOR *polygons_new;
	long j;

	polygons_new=new_polygonvector(polygons->nh);

	for (j=polygons->nl;j<=polygons->nh;j++) {
		polygons_new->element[j]=new_polygon_from_polygon(polygons->element[j]);
	}

	return polygons_new;

}

polygon_connection_attributes *new_connection_from_connection(polygon_connection_attributes *pc){
	/*
	 *
	 *
	 * \param (polygon_connection_attributes *) - polygon_connection_attributes
	 *
	 *\author Emanuele Cordano
	 *\date April 2010

	 */
	polygon_connection_attributes *pc_new;
	char *function_name="new_connection_from_connection";
	long j;

	pc_new=(polygon_connection_attributes *)malloc((sizeof(polygon_connection_attributes)));
	if (!pc_new) printf("Error in %s : polygon_connection_attributes was not allocated \n",function_name);

	pc_new->connections=new_longvector(pc->connections->nh);
	pc_new->d_connections=new_doublevector(pc->d_connections->nh);

	for(j=pc_new->connections->nl;j<=pc_new->connections->nh;j++) {
		pc_new->connections->element[j]=pc->connections->element[j];
	}

	for(j=pc_new->d_connections->nl;j<=pc_new->d_connections->nh;j++) {
		pc_new->d_connections->element[j]=pc->d_connections->element[j];
	}


	return pc_new;
}


polygon_connection_attribute_array *new_connection_array_from_connection_array(polygon_connection_attribute_array *pca) {
	/*!
	 *
	 * \param (polygon_connection_attribute_array *) - polygon_connection_attribute_array
	 *
	 *\author Emanuele Cordano
	 *\date April 2010
	 *
	 *
	 */
	long j;
	polygon_connection_attribute_array *pca_new;

	pca_new=new_connection_attributes(pca->nh);

	for(j=pca_new->nl;j<=pca_new->nh;j++) {
		pca_new->element[j]=new_connection_from_connection(pca->element[j]);
	}

	return pca_new;

}


GRID *new_grid_from_grid(GRID *grid) {
	/*!
	 *
	 * \author Emanuele Cordano
	 *
	 * \date April 2010
	 */

	GRID *grid_new;
	char *function_name="new_grid_from_grid";

	grid_new=(GRID *)malloc(sizeof(GRID));
	if (!grid_new) printf("Error in %s, grid was not allocated",function_name);

	grid_new->lines=new_linevector_from_linevector(grid->lines);

	grid_new->polygons=new_polygonvector_from_polygonvector(grid->polygons);
//	free_polygonvector(grid_new->polygons);
//	free_polygonvector(grid->polygons);
//	printf("ba!!");
//	stop_execution();

	grid_new->links=new_connection_array_from_connection_array(grid->links);
//	free_polygon_connection_attribute_array(grid->links);
//	free_polygon_connection_attribute_array(grid_new->links);
//	printf("ba2!!");
//	stop_execution();
	grid_new->boundary_indicator=grid->boundary_indicator;

	grid_new->file_resume_lines=copy_string(grid->file_resume_lines);
	grid_new->file_resume_polygons=copy_string(grid->file_resume_polygons);
	grid_new->file_resume_connections=copy_string(grid->file_resume_connections);
//	free(grid_new->file_resume_lines);
//	free(grid->file_resume_lines);
//	printf("ba24!!");
//	stop_execution();


	return grid_new;

}


/*! The following functions transforms a DOUBLESQUARE_GRID into a DOUBLEGRID */

LONGBIN *new_longbin_from_doublematrix_array(LONGMATRIX_VECTOR *lmv) {
	/*!
	 *
	 * \param (LONGMATRIX_VECTOR * ) - lmv
	 * \param (long) - novalue
	 * \author Emanuele Cordano
	 *
	 */
	LONGBIN *lb;
	LONGVECTOR *index;
	long j,r,c,cnt;
	long initialization_value=-99;
	char *function_name="new_longbin_from_doublematrix_array";

	index=new_longvector(lmv->nh);

	for(j=lmv->nl;j<=lmv->nh;j++) {
		cnt=0;
		for(r=lmv->element[j]->nrl;r<=lmv->element[j]->nrh;r++) {
			for(c=lmv->element[j]->ncl;c<=lmv->element[j]->nch;c++) {
				if ((lmv->element[j]->element[r][c]>0) && (lmv->element[j]->element[r][c]==lmv->element[j]->element[r][c])) cnt++;
			}
		}
		index->element[j]=cnt;
	}

	lb=new_longbin(index);
	for (j=lb->index->nl;j<=lb->index->nh;j++) {
		for (r=NL;r<=index->element[j];r++) {
			lb->element[j][r]=initialization_value;
		}
	}

	for (j=lb->index->nl;j<=lb->index->nh;j++) {
		cnt=0;
		for(r=lmv->element[j]->nrl;r<=lmv->element[j]->nrh;r++) {
			for(c=lmv->element[j]->ncl;c<=lmv->element[j]->nch;c++) {
				if ((lmv->element[j]->element[r][c]>0) && (lmv->element[j]->element[r][c]==lmv->element[j]->element[r][c])) {
					cnt++;
					// andar avanti !!!
					if (cnt<=lb->index->element[j]) {
					//	printf("Note n %s:filling %ld %ld on %ld res=%ld \n",function_name,j,cnt,lb->index->element[j],lb->element[j][cnt]);
						lb->element[j][cnt]=lmv->element[j]->element[r][c];
					//	printf("Note n %s:filling %ld %ld on %ld res=%ld \n",function_name,j,cnt,lb->index->element[j],lb->element[j][cnt]);
					} else {
						printf("Error in %s  : %ld at j=%ld exceeds number of elements %ld \n",function_name,cnt,j,lb->index->element[j]);
					}
				}
			}
		}

	}
	/* verify of the longbin */
	for (j=lb->index->nl;j<=lb->index->nh;j++) {
		for (r=NL;r<=index->element[j];r++) {
			if (lb->element[j][r]==initialization_value) printf("Error in %s: lb[%ld][%ld] was not set correctly \n",function_name,r,j);
		}
	}
	free_longvector(index);

	return lb;

}

LONGBIN *new_longbin_from_longbin(LONGBIN *lb) {
	/*!
	 *
	 * \param (LONGBIN *) - lb
	 *
	 * \author Emanuele Cordano
	 *
	 * \date April 2010
	 */
	//LONGVECTOR *index;
	LONGBIN *lb_new;
	long r,j;

	lb_new=new_longbin(lb->index);

	for(r=lb_new->index->nl;r<=lb_new->index->nh;r++) {
		for (j=NL;j<=lb_new->index->element[r];j++) {
			lb_new->element[r][j]=lb->element[r][j];
		}
	}

	return lb_new;
}


DOUBLE_GRID *new_double_grid_from_doublesquare_grid(DOUBLESQUARE_GRID *dsq) {
	/*!
	 *
	 *
	 * \author Emanuele Cordano
	 *
	 * \date April 2010
	 *
	 */
	 DOUBLE_GRID *dgrid;
	 char *function_name="new_double_grid_from_doublesquare_grid";

	 dgrid=(DOUBLE_GRID *)malloc(sizeof(DOUBLE_GRID));
	 if (!dgrid) printf("Error in %s, double_grid was not allocated \n",function_name);

	 dgrid->coarse=new_grid_from_grid(dsq->big->grid);
	 dgrid->fine=new_grid_from_grid(dsq->fine->grid);


	 dgrid->small_line_content=new_longbin_from_longbin_cleaning_novalues(dsq->small_content_line);
	 dgrid->small_polygon_content=new_longbin_from_doublematrix_array(dsq->small_content_polygon);

	 dgrid->novalue=dsq->fine->novalue;

	 return dgrid;

}


void free_DOUBLE_GRID(DOUBLE_GRID *dgrid) {
	/*!
	 *
	 * \param
	 *
	 * \author Emanuele Cordano
	 * \date April 2010
	 *
	 *
	 */
	free_grid(dgrid->coarse);
	free_grid(dgrid->fine);
	free_longbin(dgrid->small_line_content);
	free_longbin(dgrid->small_polygon_content);

	free(dgrid);
}


LONGBIN *new_longbin_from_longbin_cleaning_novalues(LONGBIN *lb) {
	/*!
	 *
	 * \param (LONGBIN *) - lb
	 *
	 * \author Emanuele Cordano
	 *
	 * \date April 2010
	 */
	LONGVECTOR *index;
	LONGBIN *lb_new;
	long r,j,cnt=0;
	long novalue=-99;
	index=new_longvector(lb->index->nh);
	char *function_name="new_longbin_from_longbin_cleaning_novalues";


	for(r=lb->index->nl;r<=lb->index->nh;r++) {
		cnt=lb->index->element[r];
		index->element[r]=cnt;
		for (j=NL;j<=lb->index->element[r];j++) {
		//	if (lb->index->element[])
		//	if (lb->index->element[])
			if (lb->element[r][j]<0) cnt--;
		}
		if (cnt<NL) cnt=NL;
		index->element[r]=cnt;

//		printf("%ld,%ld \d",lb->index->element[r],index->element[r]);
//		stop_execution();

	}

	lb_new=new_longbin(index);

	for(r=lb_new->index->nl;r<=lb_new->index->nh;r++) {
		for (j=NL;j<=lb_new->index->element[r];j++) {
			lb_new->element[r][j]=novalue;
		}
	}

	for(r=lb_new->index->nl;r<=lb_new->index->nh;r++) {
		lb_new->element[r][NL]=lb->element[r][NL];

		if ((lb_new->index->element[r]>NL)) {

			cnt=NL-1;
			for (j=NL;j<=lb->index->element[r];j++) {
				if (lb->element[r][j]>0) {
					cnt++;
					if (cnt>lb_new->index->element[r]) {
						printf("Error in %s counter %ld exceeds bin size %ld at row %ld (%ld) function returns NULL\n",function_name,cnt,lb_new->index->element[r],lb->index->element[r],r);
						return NULL;
					}
					lb_new->element[r][cnt]=lb->element[r][j];

				}

			}
		}
	}

	for(r=lb_new->index->nl;r<=lb_new->index->nh;r++) {
		for (j=NL;j<=lb_new->index->element[r];j++) {
//			printf("val=%ld val=%ld r=%ld j=%ld \n",lb_new->element[r][j],lb->element[r][j],r,j);
//			stop_execution();
			if ((lb_new->element[r][j]<=0) && (lb_new->index->element[r]!=NL)) {
				printf("Error in %s, negative (null) value at r=%ld j=%ld of   (%ld %ld) the created longbin function returns NULL\n",function_name,r,j,lb_new->element[r][j],lb->element[r][j]);
				stop_execution();
			//	return NULL;
			}

		}
	}

	free_longvector(index);

	return lb_new;
}



long bubble_sort_eleveation(long *cell_index, long nh, DOUBLEVECTOR *elevation) {
	/*!
	 *
	 * \author Emanuele Cordano
	 *
	 * \date April 2001
	 *
	 *\param -(long *) cell_index :  index of the cells
	 *\param -(long )  number of the cells in the rows
	 *\param (DOUBLEVECTOR *) - map of elevation
	 *
	 *\brief sort the cell_indexes according with the increase of elevation.
	 *
	 */

	char *function_name="bubble_sort_eleveation";

	long i,mk,lk;

	int swapped=1;

	long n=nh;

	for (i=0;i<nh-1;i++) {		if ((cell_index[i]<elevation->nl) || (cell_index[i]>elevation->nh)) {
			printf("Error in %s , cell_index %ld exceeds elevation %ld %ld function return -1 \n",function_name,cell_index[i],elevation->nl,elevation->nh);
	//		stop_execution();
			return -1;
		}
	}


//	for (i=la->nl;i<=la->nh;i++) {
//		la->element[i]=i;
//	}


	do {
		swapped=0;
		n=n-1;

		for(i=0;i<=n-1;i++){

//			if (elevation->element[cell_index[i]]>cell_index[i+1]){
			if (elevation->element[cell_index[i]]<elevation->element[cell_index[i+1]]){
				/*! swap  cell_index[i] and cell_index[i+1]*/
				mk=cell_index[i];
				cell_index[i]=cell_index[i+1];
				cell_index[i+1]=mk;
				/*! swap  cell_index[i] and cell_index[i+1]*/
		//		lk=la->element[i];
		//		la->element[i]=la->element[i+1];
		//		la->element[i+1]=lk;
				swapped=1;
			}
		}
//		printf ("n=%ld nh=%ld swapped=%d \n",n,nh,swapped);

		if (n<1 && swapped==1) {
			printf("Error in %s ; cells were not correctly sorted!! \n",function_name);
		//	print_longvector_elements(v,print);
			swapped=1;
		}
	} while (swapped==1);
//	stop_execution();
	return n;

}


int bubble_sort_elevation_in_longbin(LONGBIN *lb,DOUBLEVECTOR *elevation) {
	/*!
	 *
	 * \author Emanuele Cordano
	 * \date April 2010
	 *
	 */

	char *function_name="bubble_sort_eleveation_in_longbin";
	long s;
	long r;

	for (r=lb->index->nl;r<=lb->index->nh;r++) {
		if (lb->index->element[r]>NL) s=bubble_sort_eleveation(&(lb->element[r][NL]),lb->index->element[r],elevation);
	}

	return 0;

}
