/*!
 * \file keywords_file_b.h
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
#define A_FLAG 2
#define MAX_ELEVATION_VALUE 1e+10

#define I_MORPHO_ELEVATION_COARSE 1 /*!< Digital terrain model at the coarse grid resolution */
#define I_MORPHO_ELEVATION_FINE  I_MORPHO_ELEVATION_COARSE+1 /*!< Digital terrain model at the fine grid resolution */
#define I_MORPHO_ELEVATION_TERRAINSURFACE_FINE  I_MORPHO_ELEVATION_FINE+1 /*!<  Digital terrain model (of surface) at the fine grid resolution */
#define I_INITCOND_WATERSURFACE_ELEVATION  I_MORPHO_ELEVATION_TERRAINSURFACE_FINE+1 /*!<map of initial water surface  (coarse grid),*/
#define I_INITCOND_WATERDEPTH I_INITCOND_WATERSURFACE_ELEVATION+1 /*!<map of initial water depth (optional) (coarse grid)*/
#define I_POROSITY_FINE I_INITCOND_WATERDEPTH+1/*!< map of porosity (fine grid),*/
#define I_OUTLET_COEFFICIENT_FINE I_POROSITY_FINE+1 /*!<map of the coefficient of the rating curve in the outlet q_discharge=C*h_sup^p */
#define I_OUTLET_COEFFICIENT_SURF_FINE I_OUTLET_COEFFICIENT_FINE+1  /*!<map of the coefficient of the rating curve in the outlet (surface flow)  q_discharge=C*h_sup^m */
#define I_RUNOFF_COEFFICIENT_COARSE I_OUTLET_COEFFICIENT_SURF_FINE+1/*!<map of the coefficient of the dissipation runoff coefficient  C*(velocity)^p ,*/
#define I_SOURCEMAPSERIES_COARSE I_RUNOFF_COEFFICIENT_COARSE+1     /*!<series of maps for sources at different time instant (suffix)*/
#define I_DIRICHLETMAPSERIES_COARSE I_SOURCEMAPSERIES_COARSE+1 /*!< series of maps at different for Dirichlet nodes time instant (suffix) */
#define I_PARAM_FT_FILE I_DIRICHLETMAPSERIES_COARSE+1 /*!< Fluidturtle ascii files with parameters */
#define I_TIMES_FT_FILE I_PARAM_FT_FILE+1 /*!< Fluidturtle ascii files with time information for source/rainfalls (name with extension) */
#define I_DIRICHLETTIMES_FT_FILE I_TIMES_FT_FILE+1 /*!< Fluidturtle ascii files with time information for Dirichlet nodes (name with extension) */
#define I_INITIAL_VELOCITY_FT_ARRAY I_DIRICHLETTIMES_FT_FILE+1 /*<!fluidturtle array for surface velocity,*/
#define O_RESUME_FILES I_INITIAL_VELOCITY_FT_ARRAY+1
#define O_COARSEMAP_INDEX_CELLS O_RESUME_FILES+1 /*!< Output map containing indicesc of coarse cells */
#define O_COARSEMAP_INDEX_LINES O_COARSEMAP_INDEX_CELLS+1 /*!< Output map containing indicesc of coarse lines */
#define O_FINEMAP_INDEX_CELLS O_COARSEMAP_INDEX_LINES+1 /*!< Output map containing indicesc of fine cells */
#define O_FINEMAP_INDEX_LINES O_FINEMAP_INDEX_CELLS+1 /*!< Output map containing indicesc of fine lines */
#define O_COARSEMAP_WATERSURFACE_ELEVATION O_FINEMAP_INDEX_LINES+1 /*!< Output map containing water surface elevation */
#define O_COARSEMAP_WATERMASS_ERROR O_COARSEMAP_WATERSURFACE_ELEVATION+1 /*!< Output map containing water surface mass error*/
#define NFILES O_COARSEMAP_WATERMASS_ERROR


/* number of layer for coarse map*/

#define DTM_MASK 1
#define H_MASK DTM_MASK+1
#define V_MASK H_MASK+1
#define N_MAPS  V_MASK

//#define N_MAPS H_MASK
//#define V_MASK H_MASK+1
//#define N_MAPS  V_MASK

/* number of gridded variables */
#define GRIDDED_VARIABLES_START 0                /*! initialization value  */
#define BOTTOM_ELEVATION_FINE GRIDDED_VARIABLES_START   /*! map of bottom elevation (fine grid)  */
#define BOTTOM_ELEVATION_COARSE BOTTOM_ELEVATION_FINE+1 /*! map of bottom elevation (coarse grid)  */
#define BOTTOM_ELEVATION_FLINES BOTTOM_ELEVATION_COARSE+1
#define POROSITY_FINE BOTTOM_ELEVATION_FLINES+1 /*! map of porosity (fine grid) */
#define WATER_SURFACE_ELEVATION POROSITY_FINE+1
#define WATER_MASS_ERROR WATER_SURFACE_ELEVATION+1
#define N_GRIDDED_VARIABLES WATER_MASS_ERROR
