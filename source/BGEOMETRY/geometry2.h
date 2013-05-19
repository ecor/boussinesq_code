

typedef struct {


	char *name;
	long novalue;
	GRID *coarse;
	GRID *fine;

	LONGBIN *small_polygon_content;
	LONGBIN *small_line_content;

//	char *file_resume_c_polygon;
//	char *file_resume_c_line;


} DOUBLE_GRID ;

POINT *new_point_from_point(POINT *point);

LINE *new_line_from_line(LINE *line);
LINEVECTOR *new_linevector_from_linevector(LINEVECTOR *lines);
POLYGON *new_polygon_from_polygon(POLYGON *polygon);
POLYGONVECTOR *new_polygonvector_from_polygonvector(POLYGONVECTOR *polygons);
polygon_connection_attributes *new_connection_from_connection(polygon_connection_attributes *pc);
polygon_connection_attribute_array *new_connection_array_from_connection_array(polygon_connection_attribute_array *pca);
GRID *new_grid_from_grid(GRID *grid);
LONGBIN *new_longbin_from_doublematrix_array(LONGMATRIX_VECTOR *lmv);
LONGBIN *new_longbin_from_longbin(LONGBIN *lb);
DOUBLE_GRID *new_double_grid_from_doublesquare_grid(DOUBLESQUARE_GRID *dsq);
void free_DOUBLE_GRID(DOUBLE_GRID *dgrid);
LONGBIN *new_longbin_from_doublematrix_array(LONGMATRIX_VECTOR *lmv);
LONGBIN *new_longbin_from_longbin_cleaning_novalues(LONGBIN *lb);
long bubble_sort_eleveation(long *cell_index, long nh, DOUBLEVECTOR *elevation);
int bubble_sort_elevation_in_longbin(LONGBIN *lb,DOUBLEVECTOR *elevation);
