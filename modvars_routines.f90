!
! Copyright 2013 RBINS-MUMM
!
! Licensed under the EUPL, Version 1.1 or - as soon they will be approved by
! the European Commission - subsequent versions of the EUPL (the "Licence");
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!
! http://ec.europa.eu/idabc/eupl
!
! Unless required by the applicable law or agreed to in writing, software
! distributed under the Licence is distributed on an "AS IS" basis,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and
! limitations under the Licence.goto-li


MODULE modvars_routines
!************************************************************************
!
! *modvars_routines* Utility routines for model variables and file attributes
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.x
!
! $Date: 2022-01-24 09:33:32 +0100 (Mo, 24 Jan 2022) $
!
! $Revision: 1456 $
!
! Description -
!
! Reference -
!
! Routines - inquire_var, inquire_varid, set_modfiles_atts, set_modfiles_name,
!            set_modvars_atts, set_modvars_init_flags_2d,
!            set_modvars_init_flags_3d
!
!************************************************************************
!
USE iopars
USE syspars


IMPLICIT NONE

CONTAINS

!========================================================================

SUBROUTINE inquire_var(varid,f90_name,standard_name,long_name,units,node,&
                     & vector_name,format,kind_type,nrank,global_dims,&
                     & local_dims,halo_dims,isub,varatts)
!************************************************************************
!
! *inquire_var* Obtain information about a specific model variable
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.12
!
! Description -
!
! Module calls - inquire_biovar, inquire_partvar, inquire_sedvar
!
!************************************************************************
!
USE datatypes
USE gridpars
USE modids
USE nestgrids
USE obconds
USE paralpars
USE physpars
USE structures
USE switches
USE tide
USE wavepars
USE biovars_routines, ONLY: inquire_biovar
USE partvars_routines, ONLY: inquire_partvar
USE sedvars_routines, ONLY: inquire_sedvar

!
!*Arguments
!
CHARACTER (LEN=lenname), INTENT(OUT), OPTIONAL :: f90_name
CHARACTER (LEN=lendesc), INTENT(OUT), OPTIONAL :: long_name, standard_name, &
                                                & vector_name
CHARACTER (LEN=lenunit), INTENT(OUT), OPTIONAL :: units
CHARACTER (LEN=lenformat), INTENT(OUT), OPTIONAL :: format
CHARACTER (LEN=lennode), INTENT(OUT), OPTIONAL :: node
INTEGER, INTENT(IN) :: varid
INTEGER, INTENT(IN), OPTIONAL :: isub
INTEGER, INTENT(OUT), OPTIONAL :: kind_type, nrank
TYPE (VariableAtts), INTENT(OUT), OPTIONAL :: varatts
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(4) :: halo_dims
INTEGER, INTENT(OUT), OPTIONAL, DIMENSION(MaxVarDims) :: global_dims, local_dims

!
! Name           Type    Purpose
!------------------------------------------------------------------------------
!*varid*         INTEGER Variable id
!*f90_name*      CHAR    FORTRAN name
!*standard_name* CHAR    netCDF CF compliant standard name
!*long_name*     CHAR    Long descriptive name
!*units*         CHAR    Variable unit
!*node*          CHAR    Nodal type of model array
!*vector_name*   CHAR    Associated vector name
!*format*        CHAR    FORTRAN I/O format
!*kind_type*     INTEGER FORTRAN type of model variable
!*nrank*         INTEGER Rank of model array
!*global_dims*   INTEGER Global shape of model array (excluding halos)
!*local_dims*    INTEGER Local shape of model array (excluding halos)
!*halo_dims*     INTEGER Halo size of model array (WESN directions)
!*isub*          INTEGER Sub-index number in case ivarid is a multiple-variable
!                        key id
!*varatts*       DERIVED Attributes of model variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lenname) :: fname
CHARACTER (LEN=lendesc) :: lname, sname, vname
CHARACTER (LEN=lenunit) :: unit
CHARACTER (LEN=lenformat) :: formatx 
CHARACTER (LEN=lennode) :: cnode
CHARACTER (LEN=12) :: csub
INTEGER :: dtype, nodim
INTEGER,  DIMENSION(4) :: nhdims
INTEGER,  DIMENSION(MaxVarDims) :: ngdims, nldims
TYPE (VariableAtts) :: varattsx

!
!1. Non-physical parameters
!--------------------------
!

SELECT CASE (varid)
   CASE (MaxModArids+1:MaxModArids+MaxSedArids)
      CALL inquire_sedvar(varid,fname,sname,lname,unit,cnode,vname,formatx,&
                        & dtype,nodim,ngdims,nldims,nhdims,varattsx)
   CASE (MaxModArids+MaxSedArids+1:MaxModArids+MaxSedArids+MaxBioArids)
      CALL inquire_biovar(varid,fname,sname,lname,unit,cnode,vname,formatx,&
                        & dtype,nodim,ngdims,nldims,nhdims,varattsx)
   CASE (MaxModArids+MaxSedArids+MaxBioArids+1:MaxTotArids)      
      CALL inquire_partvar(varid,fname,sname,lname,unit,cnode,vname,formatx,&
                         & dtype,nodim,ngdims,nldims,nhdims,varattsx)
END SELECT

IF (varid.GT.MaxModArids.AND.varid.LE.MaxTotArids) THEN
 
   IF (PRESENT(f90_name)) THEN
      IF (PRESENT(isub)) THEN
         WRITE (csub,'(I12)') isub
         csub = ADJUSTL(csub)
         f90_name = TRIM(fname)//'_'//TRIM(csub)
      ELSE
         f90_name = fname
      ENDIF
   ENDIF

   IF (PRESENT(standard_name)) standard_name = sname
   IF (PRESENT(long_name)) long_name = lname
   IF (PRESENT(units)) units = unit
   IF (PRESENT(format)) format = formatx
   IF (PRESENT(node)) node = cnode
   IF (PRESENT(vector_name)) vector_name = vname
   IF (PRESENT(kind_type)) kind_type = dtype
   IF (PRESENT(nrank)) nrank = nodim
   IF (PRESENT(global_dims)) global_dims = ngdims
   IF (PRESENT(halo_dims)) halo_dims = nhdims
   IF (PRESENT(local_dims)) local_dims = nldims
   IF (PRESENT(varatts)) varatts = varattsx

   RETURN

ENDIF

!
!2. Initialise
!-------------
!

unit = '_'; cnode = 'C'; vname = ''; formatx = ''
dtype = float_type; nodim = 0
ngdims = 0; nhdims = 0; nldims = 0

!
!3. Get model array info
!-----------------------
!

SELECT CASE (varid)

!
!3.1 Model grid
!--------------
!

CASE (iarr_alphatc_fld)
   sname = 'factor_for_flooding_drying_scheme'
   fname = 'alphatc_fld'
   lname = 'Drying/wetting factor at C-nodes'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_alphatu_fld)
   sname = 'factor_for_flooding_drying_scheme_at_u_nodes' 
   fname = 'alphatu_fld'
   lname = 'Drying/wetting factor at U-nodes'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_alphatv_fld)
   sname = 'factor_for_flooding_drying_scheme_at_v_nodes'
   fname = 'alphatv_fld'
   lname = 'Drying/wetting factor at V-nodes'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_coriolatu)
   sname = 'coriolis_frequency_at_u_nodes'
   fname = 'coriolatu'
   lname = 'Coriolis frequency at U-nodes'
   unit = 'radian s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,1,0/)
CASE (iarr_coriolatv)
   sname = 'coriolis_frequency_at_v_nodes'
   fname = 'coriolatv'
   lname = 'Coriolis frequency at V-nodes'
   unit = 'radian s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,1/)
CASE (iarr_delxatc)
   sname = 'magnitude_of_derivative_of_position_wrt_x_coordinate_index'
   fname = 'delxatc'
   lname = 'X-spacing at C-nodes'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delxatu)
   sname = 'magnitude_of_derivative_of_position_wrt_x_coordinate_index_'//&
         & 'at_u_nodes' 
   fname = 'delxatu'
   lname = 'X-spacing at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delxatuv)
   sname = 'magnitude_of_derivative_of_position_wrt_x_coordinate_index_'//&
         & 'at_uv_nodes'
   fname = 'delxatuv'
   lname = 'X-spacing at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delxatv)
   sname = 'magnitude_of_derivative_of_position_wrt_x_coordinate_index_'//&
         & 'at_v_nodes' 
   fname = 'delxatv'
   lname = 'X-spacing at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delyatc)
   sname = 'magnitude_of_derivative_of_position_wrt_y_coordinate_index'
   fname = 'delyatc'
   lname = 'Y-spacing at C-nodes'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delyatu)
   sname = 'magnitude_of_derivative_of_position_wrt_y_coordinate_index_'//&
         & 'at_u_nodes' 
   fname = 'delyatu'
   lname = 'Y-spacing at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delyatuv)
   sname = 'magnitude_of_derivative_of_position_wrt_y_coordinate_index_'//&
         & 'at_uv_nodes' 
   fname = 'delyatuv'
   lname = 'Y-spacing at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_delyatv)
   sname = 'magnitude_of_derivative_of_position_wrt_y_coordinate_index_'//&
         & 'at_v_nodes' 
   fname = 'delyatv'
   lname = 'Y-spacing at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   nhdims = nhalo
CASE (iarr_delzatc)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number'
   fname = 'delzatc'
   lname = 'Vertical grid spacing at C-nodes'
   unit = 'm'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_delzatu)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_u_nodes' 
   fname = 'delzatu'
   lname = 'Vertical grid spacing at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,1,0,0/)
CASE (iarr_delzatuv)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_uv_nodes' 
   fname = 'delzatuv'
   lname = 'Vertical grid spacing at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,1,0,1/)
CASE (iarr_delzatuw)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_uw_nodes' 
   fname = 'delzatuw'
   lname = 'Vertical grid spacing at UW-nodes'
   unit = 'm'
   cnode = 'UW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/0,1,0,0/)
CASE (iarr_delzatv)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_v_nodes' 
   fname = 'delzatv'
   lname = 'Vertical grid spacing at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,1,1/)
CASE (iarr_delzatvw)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_vw_nodes' 
   fname = 'delzatvw'
   lname = 'Vertical grid spacing at VW-nodes'
   unit = 'm'
   cnode = 'VW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/0,0,0,1/)
CASE (iarr_delzatw)
   sname = 'magnitude_of_derivative_of_position_wrt_model_level_number_'//&
         & 'at_w_nodes' 
   fname = 'delzatw'
   lname = 'Vertical grid spacing at W-nodes'
   unit = 'm'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_dryfac)
   sname = 'dry_area_fraction'
   fname = 'dryfac'
   lname = 'Dry area fraction'
   unit = '1'
   cnode = ''
   nodim = 0
CASE (iarr_gaccatc)
   sname = 'acceleration_at_c_nodes_due_to_gravity' 
   fname = 'gaccatc'
   lname = 'Acceleration of gravity at C-nodes'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gaccatu)
   sname = 'acceleration_at_u_nodes_due_to_gravity' 
   fname = 'gaccatu'
   lname = 'Acceleration of gravity at U-nodes'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_gaccatv)
   sname = 'acceleration_at_v_nodes_due_to_gravity' 
   fname = 'gaccatv'
   lname = 'Acceleration of gravity at V-nodes'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,0,1/)
CASE (iarr_gangleatc)
   sname = 'angle_of_model_grid_wrt_reference_grid'
   fname = 'gangleatc'
   lname = 'Grid angle at C-nodes'
   unit = 'radian'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_garea)
   sname = 'horizontal_cell_area'
   fname = 'garea'
   lname = 'Horizontal grid cell area'
   unit = 'm2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_gdelxglb)
   sname = 'rectangular_grid_spacings_in_x_direction'
   fname = 'gdelxglb'
   lname = 'Global grid spacings in X-direction for rectangular grid' 
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   nodim = 1
   ngdims(1) = nc
   nhdims(1:2) = 1
CASE (iarr_gdelyglb)
   sname = 'rectangular_grid_spacings_in_y_direction'
   fname = 'gdelyglb'
   lname = 'Global grid spacings in Y-direction for rectangular grid' 
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   nodim = 1
   ngdims(1) = nr
   nhdims(1:2) = 1
CASE (iarr_gscoordatc)
   sname = 'ocean_sigma_coordinate_at_c_nodes'
   fname = 'gscoordatc'
   lname = 'Sigma coordinates at C-nodes'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_gscoordatu)
   sname = 'ocean_sigma_coordinate_at_u_nodes' 
   fname = 'gscoordatu'
   lname = 'Sigma coordinates at U-nodes'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_gscoordatuvw)
   sname = 'local_ocean_sigma_coordinate_at_uvw_nodes' 
   fname = 'gscoordatuvw'
   lname = 'Local sigma coordinates at W-nodes'
   cnode = 'UVW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gscoordatuw)
   sname = 'ocean_sigma_coordinate_at_uw_nodes' 
   fname = 'gscoordatuw'
   lname = 'Sigma coordinates at UW-nodes'
   cnode = 'UW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gscoordatv)
   sname = 'ocean_sigma_coordinate_at_v_nodes' 
   fname = 'gscoordatv'
   lname = 'Sigma coordinates at V-nodes'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_gscoordatvw)
   sname = 'ocean_sigma_coordinate_at_vw_nodes' 
   fname = 'gscoordatvw'
   lname = 'Sigma coordinates at VW-nodes'
   cnode = 'VW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gscoordatw)
   sname = 'ocean_sigma_coordinate_at_w_nodes' 
   fname = 'gscoordatw'
   lname = 'Sigma coordinates at W-nodes'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gscoordglb)
   sname = 'global_ocean_sigma_coordinate_at_w_nodes'
   fname = 'gscoordglb'
   lname = 'Global sigma coordinates at W-nodes'
   cnode = 'W'
   nodim = 3
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = 1
CASE (iarr_gsigcoordatc)
   sname = 'vertical_ocean_sigma_coordinate'
   fname = 'gsigcoordatc'
   lname = 'Sigma coordinates at centered nodes on uniform grid'
   cnode = 'C'
   nodim = 1
   ngdims(1) = nz
   CASE (iarr_gsigcoordatw)
   sname = 'vertical_ocean_sigma_coordinate'
   fname = 'gsigcoordatw'
   lname = 'Sigma coordinates at W-nodes on uniform grid'
   cnode = 'W'
   nodim = 1
   ngdims(1) = nz+1
CASE (iarr_gxcoord)
   sname = 'local_x_coordinate_at_uv_nodes' 
   fname = 'gxcoord'
   lname = 'Local X-coordinates at UV-nodes'
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gxcoordglb)
   sname = 'global_x_coordinate_at_uv_nodes' 
   fname = 'gxcoordglb'
   lname = 'Global X-coordinates at UV-nodes'
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   cnode = 'UV'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gxlon)
   sname = 'longitude'
   fname = 'gxlon'
   lname = 'Longitude at C-nodes'
   unit = 'degree_east'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_gycoord)
   sname = 'local_y_coordinate_at_uv_nodes' 
   fname = 'gycoord'
   lname = 'Local Y-coordinates at UV-nodes'
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gycoordglb)
   sname = 'global_y_coordinate_at_uv_nodes' 
   fname = 'gycoordglb'
   lname = 'Global Y-coordinates at UV-nodes'
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   cnode = 'UV'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_gylat)
   sname = 'latitude' 
   fname = 'gylat'
   lname = 'Latitude at C-nodes'
   unit = 'degree_north'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_indexobu)
   sname = 'index_at_u_open_boundaries_for_local_to_global_mapping'
   fname = 'indexobu'
   lname = 'Index at U-open_boundaries for local to global mapping'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_indexobuprocs)
   sname = 'index_at_u_open_boundaries_for_local_to_global_mapping_per_process'
   fname = 'indexobuprocs'
   lname = 'Index at U-open boundaries for local to global mapping per process'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nprocs/)
CASE (iarr_indexobv)
   sname = 'index_at_v_open_boundaries_for_local_to_global_mapping'
   fname = 'indexobv'
   lname = 'Index at V-open_boundaries for local to global mapping'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_indexobvprocs)
   sname = 'index_at_v_open_boundaries_for_local_to_global_mapping_per_process'
   fname = 'indexobvprocs'
   lname = 'Index at V-open boundaries for local to global mapping per process'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nprocs/)
CASE (iarr_indexobx)
   sname = 'index_at_x_open_boundaries_for_local_to_global_mapping'
   fname = 'indexobx'
   lname = 'Index at X-open_boundaries for local to global mapping'
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1) = nobx
CASE (iarr_indexobxprocs)
   sname = 'index_at_x_open_boundaries_for_local_to_global_mapping_per_process'
   fname = 'indexobxprocs'
   lname = 'Index at V-open boundaries for local to global mapping per process'
   dtype = int_type
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,nprocs/)
CASE (iarr_indexoby)
   sname = 'index_at_y_open_boundaries_for_local_to_global_mapping'
   fname = 'indexoby'
   lname = 'Index at Y-open_boundaries for local to global mapping'
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1) = noby
CASE (iarr_indexobyprocs)
   sname = 'index_at_y_open_boundaries_for_local_to_global_mapping_per_process'
   fname = 'indexobyprocs'
   lname = 'Index at Y-open boundaries for local to global mapping per process'
   dtype = int_type
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,nprocs/)
CASE (iarr_iobu)
   sname = 'global_x_index_at_u_open_boundaries'
   fname = 'iobu'
   lname = 'Global X-index of U-open boundaries' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1:1) = nobu
CASE (iarr_iobuloc)
   sname = 'local_x_index_at_u_open_boundaries'
   fname = 'iobuloc'
   lname = 'Local X-index of U-open boundaries' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1:1) = nobuloc_ext
CASE (iarr_iobv)
   sname = 'global_x_index_at_v_open_boundaries'
   fname = 'iobv'
   lname = 'Global X-index of V-open boundaries' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1:1) = nobv
CASE (iarr_iobvloc)
   sname = 'local_x_index_at_v_open_boundaries'
   fname = 'iobvloc'
   lname = 'Local X-index of V-open boundaries' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1:1) = nobvloc_ext
CASE (iarr_iobx)
   sname = 'global_x_index_at_x_open_boundaries'
   fname = 'iobx'
   lname = 'Global X-index of X-open boundaries' 
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1:1) = nobx
CASE (iarr_iobxloc)
   sname = 'local_x_index_at_x_open_boundaries'
   fname = 'iobxloc'
   lname = 'Local X-index of X-open boundaries' 
   dtype = int_type
   cnode = 'X'
   nodim = 1
   nldims(1:1) = nobxloc_ext
CASE (iarr_ioby)
   sname = 'global_x_index_at_y_open_boundaries'
   fname = 'ioby'
   lname = 'Global X-index of Y-open boundaries' 
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1:1) = noby
CASE (iarr_iobyloc)
   sname = 'local_x_index_at_y_open_boundaries'
   fname = 'iobyloc'
   lname = 'Local X-index of Y-open boundaries' 
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   nldims(1:1) = nobyloc_ext
CASE (iarr_jobu)
   sname = 'global_y_index_at_u_open_boundaries'
   fname = 'jobu'
   lname = 'Global Y-index of U-open boundaries' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1:1) = nobu
CASE (iarr_jobuloc)
   sname = 'local_y_index_at_u_open_boundaries'
   fname = 'jobuloc'
   lname = 'Local Y-index of U-open boundaries' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1:1) = nobuloc_ext
CASE (iarr_jobv)
   sname = 'global_y_index_at_v_open_boundaries'
   fname = 'jobv'
   lname = 'Global Y-index of V-open boundaries' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1:1) = nobv
CASE (iarr_jobvloc)
   sname = 'local_y_index_at_v_open_boundaries'
   fname = 'jobvloc'
   lname = 'Local Y-index of V-open boundaries' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1:1) = nobvloc_ext
CASE (iarr_jobx)
   sname = 'global_y_index_at_x_open_boundaries'
   fname = 'jobx'
   lname = 'Global Y-index of X-open boundaries' 
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1:1) = nobx
CASE (iarr_jobxloc)
   sname = 'local_y_index_at_x_open_boundaries'
   fname = 'jobxloc'
   lname = 'Local Y-index of X-open boundaries' 
   dtype = int_type
   cnode = 'X'
   nodim = 1
   nldims(1:1) = nobxloc_ext
CASE (iarr_joby)
   sname = 'global_y_index_at_y_open_boundaries'
   fname = 'joby'
   lname = 'Global Y-index of Y-open boundaries' 
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1:1) = noby
CASE (iarr_jobyloc)
   sname = 'local_y_index_at_y_open_boundaries'
   fname = 'jobyloc'
   lname = 'Local Y-index of Y-open boundaries' 
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   nldims(1:1) = nobyloc_ext
CASE (iarr_maskatc_int)
   sname = 'Local dry_mask_at_c_nodes'
   fname = 'maskatc_int'
   lname = 'Local dry mask at C-nodes'
   unit = 'log'
   dtype = log_type
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 0
CASE (iarr_mgvars_nc1procs)
   sname = 'global_x_index_of_local_sw_corner_per_process'
   fname = 'mgvars_nc1procs'
   lname = 'Global X-index of local SW corner per process' 
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nomglevels/)
CASE (iarr_mgvars_nc2procs)
   sname = 'global_x_index_of_local_se_corner_per_process'
   fname = 'mgvars_nc2procs'
   lname = 'Global X-index of local SE corner per process' 
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nomglevels/)
CASE (iarr_mgvars_nr1procs)
   sname = 'global_y_index_of_local_sw_corner_per_process'
   fname = 'mgvars_nr1procs'
   lname = 'Global Y-index of local SW corner per process' 
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nomglevels/)
CASE (iarr_mgvars_nr2procs)
   sname = 'global_y_index_of_local_se_corner_per_process'
   fname = 'mgvars_nr2procs'
   lname = 'Global Y-index of local SE corner per process' 
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nomglevels/)
CASE (iarr_ncprocs)
   sname = 'x_local_grid_dimension_per_process'
   fname = 'ncprocs'
   lname = 'Local grid dimension in X-direction per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nobuprocs)
   sname = 'number_of_local_interior_u_open_boundaries_per_process' 
   fname = 'nobuprocs'
   lname = 'Number of local interior U-open boundaries'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nobvprocs)
   sname = 'number_of_local_interior_v_open_boundaries_per_process' 
   fname = 'nobvprocs'
   lname = 'Number of local interior V-open boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nobxprocs)
   sname = 'number_of_local_interior_x_open_boundaries_per_process'
   fname = 'nobxprocs'
   lname = 'Number of local interior X-open boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nobyprocs)
   sname = 'number_of_local_interior_y_open_boundaries_per_process' 
   fname = 'nobyprocs'
   lname = 'Number of local interior Y-open boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nodeatc)
   sname = 'grid_pointer_at_c_nodes'
   fname = 'nodeatc'
   lname = 'Grid pointer at C-nodes'
   dtype = int_type
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_nodeatu)
   sname = 'grid_pointer_at_u_nodes'
   fname = 'nodeatu'
   lname = 'Grid pointer at U-nodes'
   dtype = int_type
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_nodeatuv)
   sname = 'grid_pointer_at_uv_nodes'
   fname = 'nodeatuv'
   lname = 'Grid pointer at UV-nodes'
   dtype = int_type
   cnode = 'UV'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_nodeatuw)
   sname = 'grid_pointer_at_uw_nodes'
   fname = 'nodeatuw'
   lname = 'Grid pointer at UW-nodes'
   dtype = int_type
   cnode = 'UW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_nodeatv)
   sname = 'grid_pointer_at_v_nodes'
   fname = 'nodeatv'
   lname = 'Grid pointer at V-nodes'
   dtype = int_type
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_nodeatvw)
   sname = 'grid_pointer_at_vw_nodes'
   fname = 'nodeatvw'
   lname = 'Grid pointer at VW-nodes'
   dtype = int_type
   cnode = 'VW'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_node2du)
   sname = '2d_grid_pointer_at_u_nodes'
   fname = 'node2du'
   lname = '2-D grid pointer at U-nodes'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_node2duv)
   sname = '2d_grid_pointer_at_uv_nodes'
   fname = 'node2duv'
   lname = '2-D grid pointer at UV-nodes'
   dtype = int_type
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_node2dv)
   sname = '2d_grid_pointer_at_v_nodes'
   fname = 'node2dv'
   lname = '2-D grid pointer at V-nodes'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_nosbuprocs)
   sname = 'number_of_local_interior_u_open_sea_boundaries_per_process' 
   fname = 'nosbuprocs'
   lname = 'Number of local interior U-open sea boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nosbvprocs)
   sname = 'number_of_local_interior_v_open_sea_boundaries_per_process' 
   fname = 'nosbvprocs'
   lname = 'Number of local interior V-open sea boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nrprocs)
   sname = 'y_local_grid_dimension_per_process'
   fname = 'nrprocs'
   lname = 'Local grid dimension in Y-direction'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nrvbuprocs)
   sname = 'number_of_local_interior_u_river_boundaries_per_process' 
   fname = 'nrvbuprocs'
   lname = 'Number of local interior U-river boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nrvbvprocs)
   sname = 'number_of_local_interior_v_river_boundaries_per_process' 
   fname = 'nrvbvprocs'
   lname = 'Number of local interior V-river boundaries per process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_rlxobcatu)
   sname = 'relaxation_factor_for_momentum_advection_at_u_nodes'
   fname = 'rlxobcatu'
   lname = 'Relaxation factor for momentum advection at U-nodes'
   cnode ='U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_rlxobcatv)
   sname = 'relaxation_factor_for_momentum_advection_at_v_nodes'
   fname = 'rlxobcatv'
   lname = 'Relaxation factor for momentum advection at V-nodes'
   cnode ='V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_seapoint)
   sname = 'Local land_mask_at_c_nodes'
   fname = 'seapoint'
   lname = 'Local land mask'
   unit = 'log'
   dtype = log_type
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_seapointglb)
   sname = 'global_land_mask_at_c_nodes'
   fname = 'seapointglb'
   lname = 'Global land mask'
   unit = 'log'
   dtype = log_type
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_soutobv)
   sname = 'south_open_boundary_at_v_node'
   fname = 'soutobv'
   lname = 'South open boundary at V-node'
   unit = 'log'
   dtype = log_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_soutoby)
   sname = 'south_open_boundary_at_y_node'
   fname = 'soutoby'
   lname = 'South open boundary at Y-node'
   unit = 'log'
   dtype = log_type
   cnode = 'UV'
   nodim = 1
   ngdims(1) = noby
CASE (iarr_westobu)
   sname = 'west_open_boundary_at_u_node'
   fname = 'westobu'
   lname = 'West open boundary at U-node'
   unit = 'log'
   dtype = log_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_westobx)
   sname = 'west_open_boundary_at_x_node'
   fname = 'westobx'
   lname = 'West open boundary at X-node'
   unit = 'log'
   dtype = log_type
   cnode = 'UV'
   nodim = 1
   ngdims(1) = nobx

!
!3.2 Depths
!----------
!

CASE (iarr_depmeanatc)
   sname = 'sea_floor_depth_below_mean_sea_level'
   fname = 'depmeanatc'
   lname = 'Mean water depth'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_depmeanatu)
   sname = 'sea_floor_depth_below_mean_sea_level_at_u_nodes'
   fname = 'depmeanatu'
   lname = 'Mean water depth at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,0,0/)
CASE (iarr_depmeanatuv)
   sname = 'sea_floor_depth_below_mean_sea_level_at_uv_nodes'
   fname = 'depmeanatuv'
   lname = 'Mean water depth at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,1/)
CASE (iarr_depmeanatv)
   sname = 'sea_floor_depth_below_mean_sea_level_at_v_nodes'
   fname = 'depmeanatv'
   lname = 'Mean water depth at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,1/)
CASE (iarr_depmeanglb)
   sname = 'global_sea_floor_depth_below_mean_sea_level'
   fname = 'depmeanglb'
   lname = 'Global mean water depth'
   nhdims = 1
   unit = 'm'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_deptotatc)
   sname = 'sea_floor_depth_below_sea_surface'
   fname = 'deptotatc'
   lname = 'Total water depth'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_deptotatc_err)
   sname = 'error_in_sea_floor_depth_below_sea_surface_due_to_mass_'//&
         & 'conservation_violation'
   fname = 'deptotatc_err'
   lname = 'Total water depth error'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_deptotatc_old)
   sname = 'sea_floor_depth_below_sea_surface_at_old_time_level'
   fname = 'deptotatc_old'
   lname = 'Total water depth at old time level'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_deptotatu)
   sname = 'sea_floor_depth_below_sea_surface_at_u_nodes'
   fname = 'deptotatu'
   lname = 'Total water depth at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)
CASE (iarr_deptotatu_old)
   sname = 'sea_floor_depth_below_sea_surface_at_u_nodes_at_old_time_level'
   fname = 'deptotatu_old'
   lname = 'Total water depth at U-nodes and old time level'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,0,0/)
CASE (iarr_deptotatuv)
   sname = 'sea_floor_depth_below_sea_surface_at_uv_nodes'
   fname = 'deptotatuv'
   lname = 'Total water depth at UV-nodes'
   unit = 'm'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,1/)
CASE (iarr_deptotatv)
   sname = 'sea_floor_depth_below_sea_surface_at_v_nodes'
   fname = 'deptotatv'
   lname = 'Total water depth at v-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,nhalo,nhalo/)
CASE (iarr_deptotatv_old)
   sname = 'sea_floor_depth_below_sea_surface_at_v_nodes_at_old_time_level'
   fname = 'deptotatv_old'
   lname = 'Total water depth at V-nodes and old time level'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,1/)
CASE (iarr_dzeta)
   sname = 'sea_surface_elevation_difference_between_new_and_previous_'//&
         & 'iteration_level'
   fname = 'dzeta'
   lname = 'Surface elevation difference at new and previous iteration level'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_zeta)
   sname = 'sea_surface_height_above_geoid'
   fname = 'zeta'
   lname = 'Surface elevation'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1
CASE (iarr_zeta_old)
   sname = 'sea_surface_height_above_geoid_at_old_time_level'
   fname = 'zeta_old'
   lname = 'Surface elevation at old time level'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = 1

!
!3.3 Currents
!------------
!

CASE (iarr_hdvelmag)
   sname = 'integral_wrt_depth_of_horizontal_sea_water_velocity'
   fname = 'hdvelmag'
   lname = 'Magnitude of depth-integrated horizontal current'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hmvelmag)
   sname = 'depth_mean_magnitude_of_horizontal_sea_water_velocity'
   fname = 'hmvelmag'
   lname = 'Magnitude of depth-mean horizontal current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hmvelmag_hadv_cfl)
   sname = 'CFL_number_of_depth_mean_horizontal_current'
   fname = 'hmvelmag_hadv_cfl'
   lname = 'CFL number depth mean horizontal current'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hvelmag)
   sname = 'magnitude_of_horizontal_sea_water_velocity'
   fname = 'hvelmag'
   lname = 'Magnitude of horizontal current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_hvelmag_hadv_cfl)
   sname = 'CFL_number_horizontal_sea_water_velocity'
   fname = 'hvelmag_hadv_cfl'
   lname = 'CFL number 3-D horizontal current'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_p2dbcgradatu)
   sname = 'integral_wrt_depth_of_x_derivative_of_kinematic_baroclinic_'//&
         & 'ocean_pressure'
   fname = 'p2dbcgradatu'
   lname = 'X-component of depth-integrated kinematic baroclinic '//&
         & 'pressure gradient'
   vname = 'Depth-integrated kinematic baroclinic pressure gradient'
   unit = 'm2 s-2'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_p2dbcgradatv)
   sname = 'integral_wrt_depth_of_y_derivative_of_kinematic_baroclinic_'//&
         & 'ocean_pressure'
   fname = 'p2dbcgradatv'
   lname = 'Y-component of depth-integrated kinematic baroclinic '//&
         & 'pressure gradient'
   vname = 'Depth-integrated kinematic baroclinic pressure gradient'
   unit = 'm2 s-2'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_p3dbcgradatu)
   sname = 'x_derivative_of_kinematic_baroclinic_ocean_pressure'
   fname = 'p3dbcgradatu'
   lname = 'X-component of kinematic baroclinic pressure gradient'
   vname = 'Kinematic baroclonic pressure gradient'
   unit = 'm s-2'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_p3dbcgradatv)
   sname = 'y_derivative_of_kinematic_baroclinic_ocean_pressure'
   fname = 'p3dbcgradatv'
   lname = 'Y-component of kinematic baroclinic pressure gradient'
   vname = 'Kinematic baroclonic pressure gradient'
   unit = 'm s-2'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_udevint)
   sname = 'baroclinic_term_in_U_momentum_equation' 
   fname = 'udevint'
   lname = 'Baroclinic term in U-momentum equation'
   vname = 'Baroclinic term in 2-D momentum equations'
   unit = 'm2 s-2'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_udfvel)
   sname = 'integral_wrt_depth_of_x_sea_water_filtered_velocity'
   fname = 'udfvel'
   lname = 'X-component of filtered depth-integrated current'
   vname = 'Filtered depth-integrated current'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_udvel)
   sname = 'integral_wrt_depth_of_x_sea_water_velocity'
   fname = 'udvel'
   lname = 'X-component of depth-integrated current'
   vname = 'Depth-integrated current'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_udvel_old)
   sname = 'integral_wrt_depth_of_x_sea_water_velocity_at_old_time_level'
   fname = 'udvel_old'
   lname = 'X-component of depth-integrated current at old time level'
   vname = 'Depth-integrated current at old time level'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_ufvel)
   sname = 'x_sea_water_filtered_velocity'
   fname = 'ufvel'
   lname = 'X-component of filtered current'
   vname = 'Filtered current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/nhalo-1,nhalo,0,0/)
CASE (iarr_umpred)
   sname = 'depth_mean_predicted_x_sea_water_velocity'
   fname = 'umpred'
   lname = 'X-component of depth-mean predicted current'
   vname = 'Depth-mean predicted current'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,0,0/)
CASE (iarr_umvel)
   sname = 'depth_mean_x_sea_water_velocity'
   fname = 'umvel'
   lname = 'X-component of depth-mean current'
   vname = 'Depth-mean current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)
CASE (iarr_umvel_hadv_cfl)
   sname = 'CFL_number_of_depth_mean_x_sea_water_velocity'
   fname = 'umvel_hadv_cfl'
   lname = 'CFL number X-component of depth-mean current'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)
CASE (iarr_umvel_old)
   sname = 'depth_mean_x_sea_water_velocity_at_old_time_level'
   fname = 'umvel_old'
   lname = 'X-component of depth-mean current at old time level'
   vname = 'Depth-mean current'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)
CASE (iarr_uvel)
   sname = 'x_sea_water_velocity'
   fname = 'uvel'
   lname = 'X-component of current'
   vname = 'Current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_uvel_hadv_cfl)
   sname = 'CFL_number_of_x_sea_water_velocity'
   fname = 'uvel_hadv_cfl'
   lname = 'CFL number X-component of 3-D current'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_uvel_old)
   sname = 'x_sea_water_velocity_at_old_time_level'
   fname = 'uvel_old'
   lname = 'X-component of current at old time level'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_vdevint)
   sname = 'baroclinic_term_in_V_momentum_equation' 
   fname = 'vdevint'
   lname = 'Baroclinic term in V-momentum equation'
   vname = 'Baroclinic term in 2-D momentum equations'
   unit = 'm2 s-2'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vdfvel)
   sname = 'integral_wrt_depth_of_y_sea_water_filtered_velocity'
   fname = 'vdfvel'
   lname = 'Y-component of filtered depth-integrated current'
   vname = 'Filtered depth-integrated current'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,0,1/)
CASE (iarr_vdvel)
   sname = 'integral_wrt_depth_of_y_sea_water_velocity'
   fname = 'vdvel'
   lname = 'Y-component of depth-integrated current'
   vname = 'Depth-integrated current'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_vdvel_old)
   sname = 'integral_wrt_depth_of_y_sea_water_velocity_at_old_time_level'
   fname = 'vdvel_old'
   lname = 'Y-component of depth-integrated current at old time level'
   vname = 'Depth-integrated current at old time level'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
CASE (iarr_vel2d)
   sname = 'integral_wrt_depth_of_sea_water_velocity_at_nested_boundaries'  
   fname = 'vel2d'
   lname = 'Depth-integrated current at nested boundaries'
   unit = 'm2 s-1'
   cnode = ''
   nodim = 2
CASE (iarr_vel3d)
   sname = 'baroclinic_sea_water_velocity_at_nested_boundaries' 
   fname = 'vel3d'
   lname = 'Baroclinic current at nested boundaries'
   unit = 'm s-1'
   cnode = ''
   nodim = 3
CASE (iarr_vfvel)
   sname = 'y_sea_water_filtered_velocity'
   fname = 'vfvel'
   lname = 'Y-component of filtered current'
   vname = 'Filtered current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,nhalo-1,nhalo/)
CASE (iarr_vmpred)
   sname = 'depth_mean_predicted_x_sea_water_velocity'
   fname = 'vmpred'
   lname = 'Y-component of depth-mean predicted current'
   vname = 'Depth-mean predicted current'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,1/)
CASE (iarr_vmvel)
   sname = 'depth_mean_y_sea_water_velocity'
   fname = 'vmvel'
   lname = 'Y-component of depth-mean current'
   vname = 'Depth-mean current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,nhalo,nhalo/)
CASE (iarr_vmvel_hadv_cfl)
   sname = 'CFL_number_of_depth_mean_y_sea_water_velocity'
   fname = 'vmvel_hadv_cfl'
   lname = 'CFL number Y-component of depth-mean current'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,1/)   
CASE (iarr_vmvel_old)
   sname = 'depth_mean_y_sea_water_velocity_at_old_time_level' 
   fname = 'vmvel_old'
   lname = 'Y-component of depth-mean current at old time level'
   vname = 'Depth-mean current at old time level'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,1,nhalo,nhalo/)
CASE (iarr_vvel)
   sname = 'y_sea_water_velocity'
   fname = 'vvel'
   lname = 'Y-component of current'
   vname = 'Current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_vvel_hadv_cfl)
   sname = 'CFL_number_of_y_sea_water_velocity'
   fname = 'vvel_hadv_cfl'
   lname = 'CFL number Y-component of 3-D current'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo   
CASE (iarr_vvel_old)
   sname = 'y_sea_water_velocity_at_old_time_level'
   fname = 'vvel_old'
   lname = 'Y-component of current at old time level'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_wphys)
   sname = 'physical_upward_sea_water_velocity'
   fname = 'wphys'
   lname = 'Physical vertical velocity'
   vname = 'Current'
   unit = 'm s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_wvel)
   sname = 'transformed_upward_sea_water_velocity'
   fname = 'wvel'
   lname = 'Transformed vertical velocity'
   vname = 'Current'
   unit = 'm s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wvel_vadv_cfl)
   sname = 'CFL_number_transformed_upward_sea_water_velocity'
   fname = 'wvel_vadv_cfl'
   lname = 'CFL number transformed vertical velocity'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
   
!
!3.4 Density
!-----------
!

CASE (iarr_beta_sal)
   sname = 'salinity_coefficient_for_sea_water_expansion'
   fname = 'beta_sal'
   lname = 'Salinity expansion coefficient'
   unit = 'PSU-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_beta_temp)
   sname = 'temperature_coefficient_for_sea_water_expansion'
   fname = 'beta_temp'
   lname = 'Temperature expansion coefficient'
   unit = 'degC-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_dens)
   sname = 'sea_water_density'
   fname = 'dens'
   lname = 'Mass density'
   unit = 'kg m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 1
CASE (iarr_sal)
   sname = 'sea_water_salinity'
   fname = 'sal'
   lname = 'Salinity'
   unit = 'PSU'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_temp)
   sname = 'sea_water_temperature'
   fname = 'temp'
   lname = 'Temperature'
   unit = 'degC'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo

!
!3.5 Diffusion coefficients
!--------------------------
!

CASE (iarr_hdifcoef2datc)
   sname = 'depth_mean_ocean_momentum_diffusivity'
   fname = 'hdifcoef2datc'
   lname = 'Depth-mean diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef2d_mom)
   sname = 'depth_mean_ocean_momentum_diffusivity'
   fname = 'hdifcoef2d_mom'
   lname = 'Depth-mean momentum diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef2d_scal)
   sname = 'depth_mean_ocean_scalar_diffusivity'
   fname = 'hdifcoef2d_scal'
   lname = 'Depth-mean scalar diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef2datuv)
   sname = 'integral_wrt_depth_of_ocean_diffusivity_at_uv_nodes'
   fname = 'hdifcoef2datuv'
   lname = 'Depth-integrated diffusion coefficient at UV-nodes'
   unit = 'm3 s-1'
   cnode = 'UV'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,1/)
CASE (iarr_hdifcoef3datc)
   sname = 'ocean_horizontal_diffusivity'
   fname = 'hdifcoef3datc'
   lname = 'Horizontal 3-D diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef3d_mom)
   sname = 'ocean_horizontal_momentum_diffusivity'
   fname = 'hdifcoef3d_mom'
   lname = 'Horizontal 3-D momentum diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef3d_scal)
   sname = 'ocean_horizontal_scalar_diffusivity'
   fname = 'hdifcoef3d_mom'
   lname = 'Horizontal 3-D scalar diffusion coefficient at C-nodes'
   unit = 'm2 s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,1,0/)
CASE (iarr_hdifcoef3datu)
   sname = 'ocean_horizontal_diffusivity_at_u_nodes'
   fname = 'hdifcoef3datu'
   lname = 'Horizontal 3-D diffusion coefficient at U-nodes'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,1,0,0/)
CASE (iarr_hdifcoef3datuv)
   sname = 'ocean_horizontal_diffusivity_at_uv_nodes' 
   fname = 'hdifcoef3datuv'
   lname = 'Horizontal 3-D diffusion coefficient at UV-nodes'
   unit = 'm2 s-1'
   cnode = 'UV'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,1,0,1/)
CASE (iarr_hdifcoef3datv)
   sname = 'ocean_horizontal_diffusivity_at_v_nodes'
   fname = 'hdifcoef3datv'
   lname = 'Horizontal 3-D diffusion coefficient at V-nodes'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,0,1/)
CASE (iarr_hdifcoef3datw)
   sname = 'ocean_horizontal_diffusivity_at_w_nodes'
   fname = 'hdifcoef3datv'
   lname = 'Horizontal 3-D diffusion coefficient at W-nodes'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 0
CASE (iarr_kinvisc)
   sname = 'sea_water_kinematic_viscosity'
   fname = 'kinvisc'
   lname = 'Kinematic viscosity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/1,0,1,0/)
CASE (iarr_mom_vdif_cfl)
   sname = 'CFL_number_of_vertical_momentum_dfiffusion'
   fname = 'mom_vdif_cfl'
   lname = 'CFL number for vertical momentum diffusion'
   cnode = 'C'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo
CASE (iarr_scal_vdif_cfl)
   sname = 'CFL_number_of_vertical_scalar_diffusion'
   fname = 'scal_vdif_cfl'
   lname = 'CFL number for vertical scalar diffusion'
   cnode = 'C'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = nhalo   
CASE (iarr_vdifcoefmom)
   sname = 'sea_water_turbulent_viscosity'
   fname = 'vdifcoefmom'
   lname = 'Eddy viscosity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/1,0,1,0/)
CASE (iarr_vdifcoefscal)
   sname = 'sea_water_turbulent_diffusivity'
   fname = 'vdifcoefscal'
   lname = 'Eddy diffusivity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_vdifcoefscal_norot)
   sname = 'non_rotated_component_sea_water_turbulent_diffusivity'
   fname = 'vdifcoefscal_norot'
   lname = 'Non-rotated component of eddy diffusivity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_vdifcoefscal_rot)
   sname = 'rotated_component_sea_water_turbulent_diffusivity'
   fname = 'vdifcoefscal_rot'
   lname = 'Rotated component of eddy diffusivity'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_vdifcoeftke)
   sname = 'vertical_diffusion_coefficient_for_turbulence_energy'
   fname = 'vdifcoeftke'
   lname = 'Vertical diffusion coefficient for turbulence energy'
   unit = 'm2 s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_xslopeatu_geo)
   sname = 'isolevel_slopes_in_x_direction_with_respect_sigma_surfaces_at_u_'//&
         & 'nodes'
   fname = 'xslopeatu_geo'
   lname = 'Isolevel slopes in the X-direction at U-nodes'
   cnode  = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,1,0,0/)
CASE (iarr_xslopeatu_ziso)
   sname = 'isoneutral_slopes_in_x_direction_at_u_nodes'
   fname = 'xslopeatu_iso'
   lname = 'Isolevel slopes in the X-direction at U-nodes'
   cnode  = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,1,0,0/)
CASE (iarr_xslopeatw_geo)
   sname = 'isoneutral_slopes_in_x_direction_with_respect_sigma_surfaces_at_'//&
         & 'w_nodes'
   fname = 'xslopeatw_geo'
   lname = 'Isolevel slopes in the X-direction at W-nodes'
   cnode  = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = 0
CASE (iarr_yslopeatv_geo)
   sname = 'isolevel_slopes_in_y_direction_with_respect_sigma_surfaces_at_'//&
         & 'v_nodes'
   fname = 'yslopeatv_geo'
   lname = 'Isolevel slopes in the Y-direction at V-nodes'
   cnode  = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,0,1/)
CASE (iarr_yslopeatv_ziso)
   sname = 'isoneutral_slopes_in_Y_direction_at_v_nodes'
   fname = 'yslopeatv_ziso'
   lname = 'Isoneutral slopes in the Y-direction at V-nodes'
   cnode  = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,1,1/)
CASE (iarr_yslopeatw_geo)
   sname = 'isoneutral_slopes_in_y_direction_with_respect_sigma_surfaces_at_'//&
         & 'w_nodes'
   fname = 'yslopeatw_geo'
   lname = 'Isolevel slopes in the Y-direction at W-nodes'
   cnode  = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_2D_hdif_cfl)
   sname = 'CFL_number_of_horizontal_dfiffusion'
   fname = '2D_hdif_cfl'
   lname = 'CFL number for horizontal diffusion'
   cnode = 'C'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

!
!3.6 Turbulence
!--------------
!

CASE (iarr_buofreq2)
   sname = 'square_of_buoyancy_frequency_in_water'
   fname = 'buofreq2'
   lname = 'Squared buoyancy frequency'
   unit = 's-2'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_dissip)
   sname = 'kinetic_energy_dissipation_in_sea_water'
   fname = 'dissip'
   lname = 'Dissipation of turbulence kinetic energy'
   unit = 'W kg-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_ricnum)
   sname = 'richardson_number_in_sea_water'
   fname = 'ricnum'
   lname = 'Richardson number'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_shearfreq2)
   sname = 'square_of_shear_frequency_in_sea_water'
   fname = 'shearfreq2'
   lname = 'Squared shear frequency'
   unit = 's-2'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
CASE (iarr_tke)
   sname = 'turbulent_kinetic_energy_in_sea_water'
   fname = 'tke'
   lname = 'Turbulent kinetic energy'
   unit = 'J kg-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_tke_old)
   sname = 'turbulent_kinetic_energy_in_sea_water_at_old_time_level'
   fname = 'tke_old'
   lname = 'Turbulent kinetic energy at old time level'
   unit = 'J kg-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_tkezl)
   sname = 'product_of_turbulent_kinetic_energy_and_mixing_length'
   fname = 'tkezl'
   lname = 'Turbulent kinetic energy times mixing length'
   unit = 'J m kg-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo
CASE (iarr_zlmix)
   sname = 'mixing_length_in_sea_water'
   fname = 'zlmix'
   lname = 'Mixing length'
   unit = 'm'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = nhalo

!
!3.7 Tides
!---------
!

CASE (iarr_astro_earth)
   sname = 'elasticity_factor_for_earth_tides'
   fname = 'astro_earth'
   lname = 'Elasticity factor for Earth tides'
   cnode = ''
   nodim = 1
   ngdims(1) = MaxAstroTides
   CASE (iarr_fnodal_anal)
   sname = 'nodal_factors_of_tidal_constituents_used_for_harmonic_analysis'
   fname = 'fnodal_anal'
   lname = 'Nodal factors of tidal constituents used for harmonic analysis'
   cnode = ''
   nodim = 1
   ngdims(1) = nofreqsanal
CASE (iarr_fnodal_astro)
   sname = 'nodal_factor_of_tidal_force'
   fname = 'fnodal_astro'
   lname = 'Nodal factor of tidal force'
   cnode = ''
   nodim = 1
   ngdims(1) = nconastro
CASE (iarr_fnodal_obc)
   sname = 'nodal_factors_of_tidal_constituents_at_open_boundaries'
   fname = 'fnodal_obc'
   lname = 'Nodal factors of tidal constituents at open boundaries'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_fxastro)
   sname = 'x_tidal_force'
   fname = 'fxastro'
   lname = 'X-component of tidal force'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_fyastro)
   sname = 'y_tidal_force'
   fname = 'fyastro'
   lname = 'Y-component of tidal force'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_index_astro)
   sname = 'key_ids_of_tidal_force_constituents'
   fname = 'index_astro'
   lname = 'Key ids of tidal force constituents'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nconastro
CASE (iarr_index_obc)
   sname = 'key_ids_of_tidal_force_constituents_at_open_boundaries'
   fname = 'index_obc'
   lname = 'Key ids of tidal force constituents at open boundaries'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_ispec_tides)
   sname = 'tidal_species' 
   fname = 'ispec_tides'
   lname = 'Tidal species'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = MaxConstituents
CASE (iarr_phase_anal)
   sname = 'astronomical_phases_for_harmonic_analysis'
   fname = 'phase_anal'
   lname = 'Astronomical phases for harmonin analysis'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nofreqsanal
CASE (iarr_phase_astro)
   sname = 'astronomical_phases_of_tidal_force'
   fname = 'phase_astro'
   lname = 'Astronomical phases of tidal force'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconastro
CASE (iarr_phase_obc)
   sname = 'astronomical_phases_at_open_boundaries'
   fname = 'phase_obc'
   lname = 'Astronomical phases at open boundaries'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_tidal_spectrum)
   sname = 'tidal_frequency' 
   fname = 'tidal_spectrum'
   lname = 'Tidal frequencies'
   unit = 'radian s-1'
   cnode = ''
   nodim = 1
   ngdims(1) = MaxConstituents

!
!3.8 Meteo forcing
!-----------------
!

CASE (iarr_airtemp)
   sname = 'air_temperature'
   fname = 'airtemp'
   lname = 'Air temperature'
   unit = 'degC'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_atmpres)
   sname = 'air_pressure_at_sea_level'
   fname = 'atmpres'
   lname = 'Atmospheric pressure'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_cloud_cover)
   sname = 'cloud_area_fraction'
   fname = 'cloud_cover'
   lname = 'Cloud cover'
   unit = '1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_evapminprec)
   sname = 'evaporation_minus_precipitation_flux'
   fname = 'evapminprec'
   lname = 'Evaporation minus precipitation rate'
   unit = 'kg m-2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_evaporation)
   sname = 'evaporation_flux'
   fname = 'evaporation'
   lname = 'Evaporation rate'
   unit = 'kg m-2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_meteodata)
   sname = 'meteo_forcing_data'
   fname = 'meteodata'
   lname = 'Meteo forcing data'
   nodim = 1
CASE (iarr_precipitation)
   sname = 'precipitation_flux' 
   fname = 'precipitation'
   lname = 'Precipitation rate'
   unit = 'kg m-2 s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_relhum)
   sname = 'relative_humidity'
   fname = 'relhum'
   lname = 'Relative humidity'
   unit = '1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_sst)
   sname = 'sea_surface_temperature'
   fname = 'sst'
   lname = 'Sea surface temperature'
   unit = 'degC'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_uwindatc)
   sname = 'x_wind' 
   fname = 'uwindatc'
   lname = 'X-component of surface wind'
   vname = 'Surface wind'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_uwindatc_old)
   sname = 'x_wind_at_old_time_level'
   fname = 'uwindatc_old'
   lname = 'X-component of surface wind at old time level'
   unit = 'm s-1'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vappres_air)
   sname = 'saturated_vapour_pressure' 
   fname = 'vappres_air'
   lname = 'Saturated vapour essure'
   unit = 'mbar'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vwindatc)
   sname = 'y_wind' 
   fname = 'vwindatc'
   lname = 'Y-component of surface wind'
   vname = 'Surface wind'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_vwindatc_old)
   sname = 'Y_wind_at_old_time_level'
   fname = 'vwindatc_old'
   lname = 'Y-component of surface wind at old time level'
   unit = 'm s-1'
   nodim = 2
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_windatc)
   sname = 'wind_speed'
   fname = 'windatc'
   lname = 'Wind speed'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)

!
!3.9 Waves
!---------
!

CASE (iarr_gxcoordglbwav)
   sname = 'global_x_coordinates_wave_grid'
   fname = 'gxcoordglbwav'
   lname = 'Global X-coordinates on the wave model grid'
   unit = MERGE('m          ','degree_east',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/ncwav,nrwav/)
CASE (iarr_gycoordglbwav)
   sname = 'global_y_coordinates_wave_grid'
   fname = 'gycoordglbwav'
   lname = 'Global Y-coordinates on the wave model grid'
   unit = MERGE('m           ','degree_north',iopt_grid_sph.EQ.0)
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/ncwav,nrwav/)
CASE (iarr_hbwdissipmag)
   sname = 'magnitude_of_bottom_wave_dissipation'
   fname = 'hbwdissipmag'
   lname = 'Magnitude of bottom wave dissipation'
   unit = 'm s-2'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_hmbwdissipmag)
   sname = 'magnitude_of_depth_mean_bottom_wave_dissipation'
   fname = 'hmbwdissipmag'
   lname = 'Magnitude of depth-mean bottom wave dissipation'
   unit = 'm s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hmswdissipmag)
   sname = 'magnitude_of_depth_mean_surface_wave_dissipation'
   fname = 'hmswdissipmag'
   lname = 'Magnitude of depth-mean surface wave dissipation'
   unit = 'm s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hswdissipmag)
   sname = 'magnitude_of_surface_wave_dissipation'
   fname = 'hswdissipmag'
   lname = 'Magnitude of surface wave dissipation'
   unit = 'm s-2'
   nodim = 2
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_maskglbwav)
   sname = 'global_land_mask_wave_grid'
   fname = 'maskglbwav'
   lname = 'Global land mask on the wave model grid'
   unit = 'log'
   dtype = log_type
   nodim = 2
   ngdims(1:2) = (/ncwav,nrwav/)
CASE (iarr_ubwdissipatc)
   sname = 'x_bed_wave_dissipation_force_at_c_nodes'
   fname = 'ubwdissipatc'
   lname = 'X-component of bed wave dissipation force at C-nodes'
   unit = 'm s-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,0,0/)
CASE (iarr_umbwdissipatc)
   sname = 'integral_wrt_depth_of_x_bed_wave_dissipation_force_at_c_nodes'
   fname = 'umbwdissipatc'
   lname = 'X-component of depth-integrated bed wave dissipation force at '//&
           'C-nodes'
   unit = 'm s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,0/)
CASE (iarr_umswdissipatc)
   sname = 'integral_wrt_depth_of_x_surface_wave_dissipation_force_at_c_nodes'
   fname = 'umswdissipatc'
   lname = 'X-component of depth-integrated surface wave dissipation force '//&
         & 'at C-nodes'
   unit = 'm s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,0/)
CASE (iarr_uswdissipatc)
   sname = 'x_surface_wave_dissipation_force_at_c_nodes'
   fname = 'uswdissipatc'
   lname = 'X-component of surface wave dissipation force at C-nodes'
   unit = 'm s-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/1,0,0,0/)
CASE (iarr_vbwdissipatc)
   sname = 'y_bed_wave_dissipation_force_at_c_nodes'
   fname = 'vbwdissipatc'
   lname = 'Y-component of bed wave dissipation force at C-nodes'
   unit = 'm s-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,1,0/)
CASE (iarr_vmbwdissipatc)
   sname = 'integral_wrt_depth_of_y_bed_wave_dissipation_force_at_c_nodes'
   fname = 'vmbwdissipatc'
   lname = 'Y-component of depth-integrated bed wave dissipation force at '//&
         & 'C-nodes'
   unit = 'm2 s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,0/)
CASE (iarr_vmswdissipatc)
   sname = 'integral_wrt_depth_of_y_surface_wave_dissipation_force_at_c_nodes'
   fname = 'vmswdissipatc'
   lname = 'Y-component of depth-integrated bed wave dissipation force at '//&
         & 'C-nodes'
   unit = 'm2 s-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,1,0/)
CASE (iarr_vswdissipatc)
   sname = 'y_surface_wave_dissipation_force_at_c_nodes'
   fname = 'vswdissipatc'
   lname = 'Y-component of surface wave dissipation force at C-nodes'
   unit = 'm s-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims = (/0,0,1,0/)
CASE (iarr_wavedir)
   sname = 'surface_wave_direction_wrt_reference_grid'
   fname = 'wavedir'
   lname = 'Wave direction'
   unit = 'degrees'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_waveexcurs)
   sname = 'near_bed_orbital_wave_excursion'
   fname = 'waveexcurs'
   lname = 'Near bed orbital wave excursion'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wavefreq)
   sname = 'wave_frequency'
   fname = 'wavefreq'
   lname = 'Peak wave frequency'
   unit = 'radian s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_waveheight)
   sname = 'wave_height'
   fname = 'waveheight'
   lname = 'Significant wave height'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wavenum)
   sname = 'wave_number'
   fname = 'wavenum'
   lname = 'Wave number'
   unit = 'm -1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_waveperiod)
   sname = 'wave_period'
   fname = 'waveperiod'
   lname = 'Peak wave period'
   unit = 's'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wavepres)
   sname = 'wave_induced_pressure'
   fname = 'wavepres'
   lname = 'Wave induced pressure'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_wavevel)
   sname = 'near_bed_wave_orbital_velocity'
   fname = 'wavevel'
   lname = 'Near bed wave orbital velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
  
!
!3.10 Stokes velocities
!----------------------
!

CASE (iarr_hmstokesmag)
   sname = 'depth_mean_magnitude_of_stokes_velocity'
   fname = 'hmstokesmag'
   lname = 'Magnitude of depth-mean Stokes current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hmveltotmag)
   sname = 'depth_mean_magnitude_of_total_horizontal_sea_water_velocity'
   fname = 'hmveltotmag'
   lname = 'Magnitude of total depth-mean current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_hstokesmag)
   sname = 'magnitude_of_horizontal_stokes_velocity'
   fname = 'hstokesmag'
   lname = 'Magnitude of Stokes current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_hveltotmag)
   sname = 'magnitude_of_total_horizontal_velocity'
   fname = 'hveltotmag'
   lname = 'Magnitude of total current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_stokessource2du)
   sname = 'stokes_source_terms_2d_u_momentum_equation'
   fname ='stokessource2du'
   lname = 'Stokes source term in the U-equation'
   unit = 'm2 s2'
   cnode  = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_stokessource2dv)
   sname = 'stokes_source_terms_2d_v_momentum_equation'
   fname ='stokessource2dv'
   lname = 'Stokes source term in the V-equation'
   unit = 'm2 s2'
   cnode  = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_udstokesatu)
   sname = 'integral_wrt_depth_of_x_sea_water_stokes_velocity_at_u_nodes'
   fname = 'udstokesatu'
   lname = 'X-component of depth-integrated Stokes current at U-nodes'
   vname = 'Depth-integrated Stokes current'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_umstokesatc)
   sname = 'depth_mean_x_sea_water_stokes_velocity_at_c_nodes'
   fname = 'umstokesatc'
   lname = 'X-component of depth-mean Stokes current at C-nodes'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   cnode = 'C'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
   nhdims = (/nhalo,nhalo,1,0/)
CASE (iarr_umstokesatu)
   sname = 'depth_mean_x_sea_water_stokes_velocity_at_u_nodes'
   fname = 'umstokesatu'
   lname = 'X-component of depth-mean Stokes current at U-nodes'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/nhalo,nhalo,1,0/)
CASE (iarr_umveltot)
   sname = 'depth_mean_x_sea_water_total_velocity'
   fname = 'umveltot'
   lname = 'X-component of depth-mean total current'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ustokesatc)
   sname = 'x_sea_water_stokes_velocity_at_c_nodes'
   fname = 'ustokesatc'
   lname = 'X-component of Stokes current at C-nodes'
   vname = 'Stokes current'
   unit = 'm s-1'
   cnode = 'C'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims =  (/nhalo,nhalo,1,0/)
CASE (iarr_ustokesatu)
   sname = 'x_sea_water_stokes_velocity_at_u_nodes'
   fname = 'ustokesatu'
   lname = 'X-component of Stokes current at U-nodes'
   vname = 'Stokes current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims =  (/nhalo,nhalo,1,0/)
CASE (iarr_uveltot)
   sname = 'x_sea_water_total_velocity'
   fname = 'uveltot'
   lname = 'X-component of total current'
   vname = 'Total current'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_vdstokesatv)
   sname = 'integral_wrt_depth_of_y_sea_water_stokes_velocity_at_v_nodes'
   fname = 'vdstokesatv'
   lname = 'Y-component of depth-integrated Stokes current at V-nodes'
   vname = 'Depth-integrated Stokes current'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,0,0,1/)
CASE (iarr_vmstokesatc)
   sname = 'depth_mean_y_sea_water_stokes_velocity'
   fname = 'vmstokesatc'
   lname = 'Y-component of depth-mean Stokes current'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = nhalo
   nhdims = (/1,0,nhalo,nhalo/)
CASE (iarr_vmstokesatv)
   sname = 'depth_mean_y_sea_water_stokes_velocity_at_v_nodes'
   fname = 'vmstokesatv'
   lname = 'Y-component of depth-mean Stokes current at V-nodes'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,nhalo,nhalo/)
CASE (iarr_vmveltot)
   sname = 'depth_mean_y_sea_water_total_velocity'
   fname = 'vmveltot'
   lname = 'Y-component of depth-mean total current'
   vname = 'Depth-mean Stokes current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vstokesatc)
   sname = 'y_sea_water_stokes_velocity_at_c_nodes'
   fname = 'vstokesatc'
   lname = 'Y-component of Stokes current at C-nodes'
   vname = 'Stokes current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims =  (/1,0,nhalo,nhalo/)
CASE (iarr_vstokesatv)
   sname = 'y_sea_water_stokes_velocity_at_v_nodes'
   fname = 'vstokesatv'
   lname = 'Y-component of Stokes current at V-nodes'
   vname = 'Stokes current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
   nhdims =  (/1,0,nhalo,nhalo/)
CASE (iarr_vveltot)
   sname = 'y_sea_water_total_velocity'
   fname = 'vveltot'
   lname = 'Y-component of total current'
   vname = 'Total current'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_wstokesatw)
   sname = 'transformed_upward_sea_water_stokes_velocity'
   fname = 'wstokesatw'
   lname = 'Transformed vertical Stokes velocity'
   vname = 'Stokes current'
   unit = 'm s-1'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)
   nhdims = (/1,0,1,0/)

!
!3.11 Optical arrays
!-------------------
!

CASE (iarr_optattcoef2)
   sname = 'volume_attenuation_coefficient_of_short_wave_radiative_flux_'//&
         & 'in_sea_water'
   fname = 'optattcoef2'
   lname = 'Inverse optical attenuation depth for short waves'
   unit = 'm-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qrad)
   sname = 'surface_downwelling_spherical_irradiance_in_sea_water'
   fname = 'qrad'
   lname = 'Surface solar irradiance'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_radiance)
   sname = 'downwelling_sperical_irradiance_in_sea_water'
   fname = 'radiance'
   lname = 'Solar irradiance'
   unit = 'W m-2'
   cnode = 'W'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz+1/)
   ngdims(1:3) = (/nc,nr,nz+1/)

!
!3.12 Bottom/surface fluxes
!--------------------------
!

CASE (iarr_bdragcoefatc)
   sname = 'drag_coefficient_at_sea_floor_at_c_nodes'
   fname = 'bdragcoefatc'
   lname = 'Bottom drag coefficient'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_bfricatu)
   sname = 'friction_velocity_term_at_sea_floor_at_u_nodes'
   fname = 'bfricatu'
   lname = 'Bottom friction velocity term at U-nodes'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricatv)
   sname = 'friction_velocity_term_at_sea_floor_at_v_nodes'
   fname = 'bfricatv'
   lname = 'Bottom friction velocity term at V-nodes'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricvel)
   sname = 'friction_velocity_at_sea_floor'
   fname = 'bfricvel'
   lname = 'Bottom friction velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricvel_max)
   sname = 'maximum_current_wave_friction_velocity_at_sea_floor'
   fname = 'bfricvel_max'
   lname = 'Maximum current-wave bottom friction velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bfricvel_wav)
   sname = 'wave_friction_velocity_at_sea_floor'
   fname = 'bfricvel_wav'
   lname = 'Bottom wave friction'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatc)
   sname = 'stress_at_sea_floor'
   fname = 'bstresatc'
   lname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatc_max)
   sname = 'maximum_current_wave_stress_at_sea_floor'
   fname = 'bstresatc_max'
   lname = 'Maximum current-wave bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatc_wav)
   sname = 'maximum_wave_stress_at_sea_floor'
   fname = 'bstresatc_wav'
   lname = 'Maximum wave bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatu)
   sname = 'stress_at_sea_floor'
   fname = 'bstresatu'
   lname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_bstresatv)
   sname = 'stress_at_sea_floor'
   fname = 'bstresatv'
   lname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_cds)
   sname = 'drag_coefficient_at_sea_surface'
   fname = 'cds'
   lname = 'Surface drag coefficient'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_ces)
   sname = 'exchange_coefficient_for_latent_heat_flux' 
   fname = 'ces'
   lname = 'Exchange coefficient for latent heat flux'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_chs)
   sname = 'exchange_coefficient_for_sensible_heat_flux'
   fname = 'chs'
   lname = 'Exchange coefficient for sensible heat flux'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_fwave)
   sname = 'wave_friction_factor'
   fname = 'fwave'
   lname = 'Wave friction factor'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qlatent)
   sname = 'surface_upward_latent_heat_flux'
   fname = 'qlatent'
   lname = 'Latent heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qlwave)
   sname = 'surface_net_upward_longwave_flux'
   fname = 'qlwave'
   lname = 'Long-wave heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qnonsol)
   sname = 'surface_upward_non_solar_heat_flux_in_sea_water'
   fname = 'qnonsol'
   lname = 'Non-solar heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qsensible)
   sname = 'surface_upward_sensible_heat_flux'
   fname = 'qsensible'
   lname = 'Sensible heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_qtot)
   sname = 'surface_downward_heat_flux_in_sea_water'
   fname = 'qtot'
   lname = 'Total downward surface heat flux'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_sfricatc)
   sname = 'surface_firction_velocity'
   fname = 'sfricatc'
   lname = 'Surface friction velocity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ssalflux)
!  ? + units : PSU m s-1 instead of kg m-2 s -1!!!!! 
   sname = 'virtual_salt_flux_into_sea_water' 
   fname = 'ssalflux'
   lname = 'Salinity flux'
   unit = 'PSU m s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_sstresatc)
   sname = 'surface_downward_stress'
   fname = 'sstresatc'
   lname = 'Surface stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_surdat1d)
   sname = '1d_forcing_surface_data'
   fname = 'surdat1d'
   lname = '1-D surface forcing data'
   nodim = 1
CASE (iarr_ubstresatc)
   sname = 'upward_x_stress_at_sea_floor'
   fname = 'ubstresatc'
   lname = 'X-component of bottom stress at C-nodes'
   vname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ubstresatu)
   sname = 'upward_x_stress_at_sea_floor'
   fname = 'ubstresatu'
   lname = 'X-component of bottom stress at U-nodes'
   vname = 'Bottom stress'
   unit = 'Pa'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/0,1,0,0/)
CASE (iarr_usstresatc)
   sname = 'surface_downward_x_stress'
   fname = 'usstresatc'
   lname = 'X-component of surface stress'
   vname = 'Surface stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,0/)
CASE (iarr_usstresatu)
   sname = 'surface_downward_x_stress_at_u_nodes'
   fname = 'usstresatu'
   lname = 'X-component of surface stress at U-nodes'
   vname = 'Surface stress'
   unit = 'Pa'
   cnode = 'U'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vbstresatc)
   sname = 'upward_y_stress_at_sea_floor'
   fname = 'vbstresatc'
   lname = 'Y-component of bottom stress at C-nodes'
   vname = 'Bottom stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vbstresatv)
   sname = 'upward_y_stress_at_sea_floor'
   fname = 'vbstresatv'
   lname = 'Y-component of bottom stress'
   vname = 'Bottom stres'
   unit = 'Pa'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims =(/0,0,0,1/)
CASE (iarr_vsstresatc)
   sname = 'surface_downward_x_stress'
   fname = 'vsstresatc'
   lname = 'X-component of surface stress'
   vname = 'Surface stress'
   unit = 'Pa'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,0,0/)
CASE (iarr_vsstresatv)
   sname = 'surface_downward_y_stress_at_v_nodes'
   fname = 'vsstresatv'
   lname = 'Y-component of surface stress at V-nodes'
   vname = 'Surface stress'
   unit = 'Pa'
   cnode = 'V'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_wavedata)
   sname = 'wave_forcing_data'
   fname = 'wavedata'
   lname = 'Wave forcing data'
   nodim = 1   
CASE (iarr_wavethickatc)
   sname = 'wave_boundary_layer_thickness'
   fname = 'wavethickatc'
   lname = 'Wave boundary layer thickness '
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_zaroughatc)
   sname = 'sea_floor_apparent_roughness_length' 
   fname = 'zaroughatc'
   lname = 'Bottom apparent roughness length'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)
CASE (iarr_zroughatc)
   sname = 'sea_floor_roughness_length' 
   fname = 'zroughatc'
   lname = 'Bottom roughness length'
   unit = 'm'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
   nhdims = (/1,0,1,0/)

!
!3.13 Open boundary forcing
!--------------------------
!

CASE (iarr_floutobu)
   sname = 'last_output_time_at_u_open_boundaries'
   fname = 'floutobu'
   lname = 'Last output time at U-open boundaries'
   unit = 's'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nz/)
CASE (iarr_floutobv)
   sname = 'last_output_time_at_v_open_boundaries'
   fname = 'floutobv'
   lname = 'Last output time at V-open boundaries'
   unit = 's'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nz/)
CASE (iarr_gxslope)
   sname = 'x_derivative_of_kinematic_ocean_pressure'
   fname = 'gxslope'
   lname = 'X-component of kinematic pressure gradient'
   unit = 'm s-2'
   cnode = ''
   nodim = 0
CASE (iarr_gxslope_amp)
   sname = 'amplitude_of_x_derivative_of_kinematic_ocean_pressure'
   fname = 'gxslope_amp'
   lname = 'Amplitude of X-component of kinematic pressure gradient'
   unit = 'm s-2'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_gxslope_pha)
   sname = 'phase_of_x_derivative_of_kinematic_ocean_pressure'
   fname = 'gxslope_pha'
   lname = 'Phase of X-component of kinematic pressure gradient'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_gyslope)
   sname = 'y_derivative_of_kinematic_ocean_pressure'
   fname = 'gyslope'
   lname = 'Y-component of kinematic pressure gradient'
   unit = 'm s-2'
   cnode = ''
   nodim = 0
CASE (iarr_gyslope_amp)
   sname = 'Amplitude_of_y_derivative_of_kinematic_ocean_pressure'
   fname = 'gyslope_amp'
   lname = 'Amplitude of Y-component of kinematic pressure gradient'
   unit = 'm s-2'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_gyslope_pha)
   sname = 'phase_of_y_derivative_of_kinematic_ocean_pressure'
   fname = 'gyslope_pha'
   lname = 'Phase of Y-component of kinematic pressure gradient'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_iloczobu)
   sname = 'position_of_elevation_points_with_respect_to_u_open_boundaries'
   fname = 'iloczobu'
   lname = 'Position of elevation points with respect to U-open boundaries'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_iloczobv)
   sname = 'position_of_elevation_points_with_respect_to_v_open_boundaries'
   fname = 'iloczobv'
   lname = 'Position of elevation points with respect to V-open boundaries'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_indexprof)
   sname = 'mapping_array_of_profile_numbers'
   fname = 'indexprof'
   lname = 'Map between of profile numbers'
   dtype = int_type
   cnode = ''
   nodim = 3
CASE (iarr_index2dobuv)
   sname = 'mapping_array_of_data_points_to_normal_open_boundary_locations'
   fname = 'index2dobuv'
   lname = 'Map of data points to normal open boundary locations'  
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nobu+nobv,maxdatafiles(io_2uvobc,1)-1/)
CASE (iarr_index2dobxy)
   sname = 'mapping_array_of_data_points_to_tangential_open_boundary_locations'
   fname = 'index2dobxy'
   lname = 'Map of data points to tangential open boundary locations'  
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nobx+noby,maxdatafiles(io_2xyobc,1)-1/)
CASE (iarr_iobc2dtype)
   sname = 'type_2d_open_boundary_data_input'
   fname = 'iobc2dtype'
   lname = 'Type 2-D open boundary data input'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = maxdatafiles(io_2uvobc,1)-1
CASE (iarr_iprofobu)
   sname = 'profile_number_at_u_open_boundaries'
   fname = 'iprofobu'
   lname = 'Profile number at U-open boundaries'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,1/)
CASE (iarr_iprofobv)
   sname = 'profile_number_at_v_open_boundaries'
   fname = 'iprofobv'
   lname = 'Profile number at V-open boundaries'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,1/)
CASE (iarr_iprofobx)
   sname = 'profile_number_at_x_open_boundaries'
   fname = 'iprofobx'
   lname = 'Profile number at X-open boundaries'
   dtype = int_type
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,1/)
CASE (iarr_iprofoby)
   sname = 'profile_number_at_y_open_boundaries'
   fname = 'iprofoby'
   lname = 'Profile number at Y-open boundaries'
   dtype = int_type
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,1/)
CASE (iarr_iqsecobu)
   sname = 'start_end_locations_at_u_node_discharge_boundaries'
   fname = 'iqsecobu'
   lname = 'Start/end locations at U-node discharge boundaries'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nqsecobu,2/)
CASE (iarr_itypobu)
   sname = 'type_of_open_boundary_condition_at_u_open_boundaries'
   fname = 'itypobu'
   lname = 'Type of open boundary condition at U-open boundaries'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_itypobv)
   sname = 'type_of_open_boundary_condition_at_v_open_boundaries'
   fname = 'itypobv'
   lname = 'Type of open boundary condition at V-open boundaries'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_itypobx)
   sname = 'type_of_open_boundary_condition_at_x_open_boundaries'
   fname = 'itypobx'
   lname = 'Type of open boundary condition at X-open boundaries'
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1) = nobx
CASE (iarr_itypoby)
   sname = 'type_of_open_boundary_condition_at_y_open_boundaries'
   fname = 'itypoby'
   lname = 'Type of open boundary condition at Y-open boundaries'
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1) = noby
CASE (iarr_ityp2dobu)
   sname = 'type_of_2d_open_boundary_condition_at_u_nodes'
   fname = 'ityp2dobu'
   lname = 'Type of 2-D open boundary condition at U-nodes'
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = nobu
CASE (iarr_ityp2dobv)
   sname = 'type_of_2d_open_boundary_condition_at_v_nodes'
   fname = 'ityp2dobv'
   lname = 'Type of 2-D open boundary condition at V-nodes'
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_ityp2dobx)
   sname = 'type_of_2d_open_boundary_condition_at_x_nodes'
   fname = 'ityp2dobx'
   lname = 'Type of 2-D open boundary condition at X-nodes'
   dtype = int_type
   cnode = 'X'
   nodim = 1
   ngdims(1) = nobx
CASE (iarr_ityp2doby)
   sname = 'type_of_2d_open_boundary_condition_at_y_nodes'
   fname = 'ityp2doby'
   lname = 'Type of 2-D open boundary condition at Y-nodes'
   dtype = int_type
   cnode = 'Y'
   nodim = 1
   ngdims(1) = nobv
CASE (iarr_jqsecobv)
   sname = 'start_end_locations_at_v_node_discharge_boundaries'
   fname = 'jqsecobv'
   lname = 'Start/end locations at V-node discharge boundaries'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nqsecobv,2/)
CASE (iarr_noprofsd)
   sname = 'number_of_profiles_per_input_file'
   fname = 'noprofsd'
   lname = 'Number of profiles per input file'
   dtype = int_type
   cnode = ''
   nodim = 1
CASE (iarr_no2dobuv)
   sname = 'number_of_input_data_per_input_file_at_normal_open_boundaries'
   fname = 'no2dobuv'
   lname = 'Number of input data per input file'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = maxdatafiles(io_2uvobc,1)-1
CASE (iarr_no2dobxy)
   sname = 'number_of_input_data_per_input_file_at_tangential_open_boundaries'
   fname = 'no2dobxy'
   lname = 'Number of input data per input file'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = maxdatafiles(io_2xyobc,1)-1
CASE (iarr_obcdatuv2d)
   sname = '2d_input_data_at_normal_open_boundaries'
   fname = 'obcdatuv2d'
   lname = '2-D input data at normal open boundaries'
   cnode = ''
   nodim = 2
CASE (iarr_obcdatxy2d)
   sname = '2d_input_data_at_tangential_open_boundaries'
   fname = 'obcdatxy2d'
   lname = '2-D input data at tangential open boundaries'
   unit = 'm s-1'
   cnode = ''
   nodim = 2
CASE (iarr_obcsalatu)
   sname = 'storage_array_for_salinity_at_u_open_boundaries'
   fname = 'obcsalatu'
   lname = 'Storage array for salinity at U-open boundaries'
   unit = 'PSU'
   cnode = 'U'
   nodim = 4
   ngdims(1:4) = (/nobu,nz,3,1/)
CASE (iarr_obcsalatv)
   sname = 'storage_array_for_salinity_at_v_open_boundaries'
   fname = 'obcsalatv'
   lname = 'Storage array for salinity at V-open boundaries'
   unit = 'PSU'
   cnode = 'V'
   nodim = 4
   ngdims(1:4) = (/nobv,nz,3,1/)
CASE (iarr_obctmpatu)
   sname = 'storage_array_for_temperature_at_u_open_boundaries'
   fname = 'obctmpatu'
   lname = 'Storage array for temperature at U-open boundaries'
   unit = 'degC'
   cnode = 'U'
   nodim = 4
   ngdims(1:4) = (/nobu,nz,3,1/)
CASE (iarr_obctmpatv)
   sname = 'storage_array_for_temperature_at_v_open_boundaries'
   fname = 'obctmpatv'
   lname = 'Storage array for temperature at V-open boundaries'
   unit = 'degC'
   cnode = 'V'
   nodim = 4
   ngdims(1:4) = (/nobv,nz,3,1/)
CASE (iarr_obc2uvatu)
   sname = 'storage_array_for_2d_mode_at_u_open_boundaries'
   fname = 'obc2uvatu'
   lname = 'Storage array for 2-D mode at U-open boundaries'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,2/)
CASE (iarr_obc2uvatu_old)
   sname = 'storage_array_for_2d_mode_at_u_open_boundaries_at_old_time_level'
   fname = 'obc2uvatu_old'
   lname = 'Storage array for 2-D mode at U-open boundaries and old time level'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,2/)
CASE (iarr_obc2uvatv)
   sname = 'storage_array_for_2d_mode_at_v_open_boundaries'
   fname = 'obc2uvatv'
   lname = 'Storage array for 2-D mode at V-open boundaries'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,2/)
CASE (iarr_obc2uvatv_old)
   sname = 'storage_array_for_2d_mode_at_v_open_boundaries_at_old_time_level'
   fname = 'obc2uvatv_old'
   lname = 'Storage array for 2-D mode at V-open boundaries and old time level'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,2/)
CASE (iarr_obc3uvatu)
   sname = 'storage_array_for_3d_mode_at_u_open_boundaries'
   fname = 'obc3uvatu'
   lname = 'Storage array for 3-D mode at U-open boundaries'
   unit = 'm s-1'
   cnode = 'U'
   nodim = 3
   ngdims(1:3) = (/nobu,nz,2/)
CASE (iarr_obc3uvatv)
   sname = 'storage_array_for_3d_mode_at_v_open_boundaries'
   fname = 'obc3uvatv'
   lname = 'Storage array for 3-D mode at V-open boundaries'
   unit = 'm s-1'
   cnode = 'V'
   nodim = 3
   ngdims(1:3) = (/nobv,nz,2/)
CASE (iarr_prof3dsal)
   sname = 'sea_water_salinity_at_open_boundaries'
   fname = 'prof3dsal'
   lname = 'Salinity profiles at open boundaries'
   unit = 'PSU'
   cnode = ''
   nodim = 2
CASE (iarr_prof3dtmp)
   sname = 'sea_water_temperature_at_open_boundaries'
   fname = 'prof3dtmp'
   lname = 'Temperature profiles at open boundaries'
   unit = 'degC'
   cnode = ''
   nodim = 2
CASE (iarr_prof3dveluv)
   sname = 'baroclinic_sea_water_velocity_at_normal_open_boundaries'
   fname = 'prof3dveluv'
   lname = 'Baroclinic current profiles at normal open boundaries'
   unit = 'm s-1'
   cnode = ''
   nodim = 2
CASE (iarr_prof3dvelxy)
   sname = 'baroclinic_sea_water_velocity_at_tangential_open_boundaries'
   fname = 'prof3dvelxy'
   lname = 'Baroclinic current profiles at tangential open boundaries'
   unit = 'm s-1'
   cnode = ''
   nodim = 2
CASE (iarr_return_time)
   sname = 'return_time_at_open_boundaries'
   fname = 'return_time'
   lname = 'Vertical profile of return times for TH-open boundary condition'
   unit = 's'
   cnode = ''
   nodim = 1
   ngdims(1) = nz
CASE (iarr_udatobu)
   sname = 'integral_wrt_depth_of_x_sea_water_velocity_at_u_open_boundaries'
   fname = 'udatobu'
   lname = 'X-component of depth-integrated current or discharge at '//&
         & 'U-normal boundaries'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,2/)
CASE (iarr_udatobu_amp)
   sname = 'amplitude_of_integral_wrt_depth_of_x_sea_water_velocity_at_'//&
         & 'u_open_boundaries'
   fname = 'udatobu_amp'
   lname = 'Amplitude of X-component of depth-integrated current or '//&
         & 'discharge at U-open boundaries'
   unit = 'm2 s-1'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nconobc/)
CASE (iarr_udatobu_pha)
   sname = 'phase_of_integral_wrt_depth_of_x_sea_water_velocity_at_'//&
         & 'u_open_boundaries'
   fname = 'udatobu_pha'
   lname = 'Phase of X-component of depth-integrated current or discharge '//&
         & 'at U-open boundaries'
   unit = 'radian'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nconobc/)
CASE (iarr_udatoby)
   sname = 'integral_wrt_depth_of_x_sea_water_velocity_wrt_depth_at_'//&
         & 'y_open_boundaries'
   fname = 'udatoby'
   lname = 'X-component of depth-integrated current at U-normal boundaries'
   unit = 'm2 s-1'
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,2/)
CASE (iarr_udatoby_amp)
   sname = 'amplitude_of_integral_wrt_depth_of_x_sea_water_velocity_at_'//&
         & 'y_open_boundaries'
   fname = 'udatoby_amp'
   lname = 'Amplitude of X-component of depth-integrated current at Y-open '//&
         & 'boundaries'
   unit = 'm2 s-1'
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,nconobc/)
CASE (iarr_udatoby_pha)
   sname = 'phase_of_integral_wrt_depth_of_x_sea_water_velocity_at_'//&
         & 'y_open_boundaries'
   fname = 'udatoby_pha'
   lname = 'Phase of X-component of depth-integrated current at Y-open '//&
         & 'boundaries'
   unit = 'radian'
   cnode = 'Y'
   nodim = 2
   ngdims(1:2) = (/noby,nconobc/)
CASE (iarr_vdatobv)
   sname = 'integral_wrt_depth_of_y_sea_water_velocity_at_v_open_boundaries' 
   fname = 'vdatobv'
   lname = 'Y-component of depth-integrated current or discharge at V-open '//&
         & 'boundaries'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,2/)
CASE (iarr_vdatobv_amp)
   sname = 'amplitude_of_integral_wrt_depth_of_y_sea_water_velocity_at_v_'//&
         & 'open_boundaries'
   fname = 'vdatobv_amp'
   lname = 'Amplitude of Y-component of depth-integrated current or '//&
         & 'discharge at V-open boundaries'
   unit = 'm2 s-1'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nconobc/)
CASE (iarr_vdatobv_pha)
   sname = 'phase_of_integral_wrt_depth_of_y_sea_water_velocity_at_v_'//&
         & 'open_boundaries'
   fname = 'vdatobv_pha'
   lname = 'Phase of Y-component of depth-integrated current or discharge '//&
         & 'at V-open boundaries'
   unit = 'radian'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nconobc/)
CASE (iarr_vdatobx)
   sname = 'integral_wrt_depth_of_y_sea_water_velocity_at_x_open_boundaries'
   fname = 'vdatobx'
   lname = 'Y-component of depth-integrated current at X-normal boundaries'
   unit = 'm2 s-1'
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,2/)
CASE (iarr_vdatobx_amp)
   sname = 'amplitude_of_integral_wrt_depth_of_y_sea_water_velocity_at_x_'//&
         & 'open_boundaries'
   fname = 'vdatobx_amp'
   lname = 'Amplitude of Y-component of depth-integrated current at X-open '//&
         & 'boundaries'
   unit = 'm2 s-1'
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,nconobc/)
CASE (iarr_vdatobx_pha)
   sname = 'phase_of_integral_wrt_depth_of_y_sea_water_velocity_at_x_'//&
         & 'open_boundaries'
   fname = 'vdatobx_pha'
   lname = 'Phase of Y-component of depth-integrated current at X-open '//&
         & 'boundaries'
   unit = 'radian'
   cnode = 'X'
   nodim = 2
   ngdims(1:2) = (/nobx,nconobc/)
CASE (iarr_zdatobu)
   sname = 'sea_surface_height_above_sea_level_at_u_open_boundaries'
   fname = 'zdatobu'
   lname = 'Surface elevation at U-open boundaries'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,2/)
CASE (iarr_zdatobu_amp)
   sname = 'sea_surface_height_amplitude_at_u_open_boundaries'
   fname = 'zdatobu_amp'
   lname = 'Amplitude of surface elevation at U-open boundaries'
   unit = 'm'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) =(/nobu,nconobc/)
CASE (iarr_zdatobu_pha)
   sname = 'sea_surface_height_phase_at_u_open_boundaries'
   fname = 'zdatobu_pha'
   lname = 'Phase of surface elevation at U-open boundaries'
   unit = 'radian'
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/nobu,nconobc/)
CASE (iarr_zdatobv)
   sname = 'sea_surface_height_at_v_open_boundaries' 
   fname = 'zdatobv'
   lname = 'Surface elevation at V-open boundaries'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,2/)
CASE (iarr_zdatobv_amp)
   sname = 'sea_surface_height_amplitude_at_v_open_boundaries'
   fname = 'zdatobv_amp'
   lname = 'Amplitude of surface elevation at V-open boundaries'
   unit = 'm'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nconobc/)
CASE (iarr_zdatobv_pha)
   sname = 'sea_surface_height_phase_at_v_open_boundaries'
   fname = 'zdatobv_pha'
   lname = 'Phase of surface elevation at V-open boundaries'
   unit = 'radian'
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/nobv,nconobc/)
CASE (iarr_zeta_amp)
   sname = 'sea_surface_height_amplitude'
   fname = 'zeta_amp'
   lname = 'Amplitude of surface elevation'
   unit = 'm'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc
CASE (iarr_zeta_pha)
   sname = 'sea_surface_height_phase'
   fname = 'zeta_pha'
   lname = 'Phase of surface elevation'
   unit = 'radian'
   cnode = ''
   nodim = 1
   ngdims(1) = nconobc

!
!3.14 Structure module
!---------------------
!

CASE (iarr_idry)
   sname = 'global_x_index_at_dry_cells'
   fname = 'idry'
   lname = 'Global X-index of dry cells' 
   dtype = int_type
   nodim = 1
   ngdims(1) = numdry
CASE (iarr_indexwbaru)
   sname = 'index_at_u_node_weirs_barriers_for_local_to_global_mapping'
   fname ='indexwbaru'
   lname = 'Index at U-node weirs and barriers for local to global mapping'
   dtype = int_type
   cnode ='U'
   nodim = 1
   nldims(1) = numwbaruloc
CASE (iarr_indexwbaruprocs)
   sname = 'index_at_u_node_weirs_barriers_for_local_to_global_mapping_'//&
         & 'per_process'
   fname = 'indexwbaruprocs'
   lname = 'Index at U-node weirs and barriers for local to global '//&
         & 'mapping per process'
   dtype = int_type
   cnode = 'U'
   nodim = 2
   ngdims(1:2) = (/numwbaru,nprocs/)
CASE (iarr_indexwbarv)
   sname = 'index_at_v_node_weirs__barriers_for_local_to_global_mapping'
   fname ='indexwbarv'
   lname = 'Index at V-node weirs and barriers for local to global mapping'
   dtype = int_type
   cnode ='V'
   nodim = 1
   nldims(1) = numwbarvloc
CASE (iarr_indexwbarvprocs)
   sname = 'index_at_v_node_weirs_barriers_for_local_to_global_mapping_'//&
         & 'per_process'
   fname = 'indexwbarvprocs'
   lname = 'Index at V-node weirs and barriers for local to global mapping '//&
         & 'per process'
   dtype = int_type
   cnode = 'V'
   nodim = 2
   ngdims(1:2) = (/numwbarv,nprocs/)
CASE (iarr_ithinu)
   sname = 'global_x_index_at_thin_dams_at_u_nodes'
   fname = 'ithinu'
   lname = 'Global X-index of U-thin dams' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = numthinu
CASE (iarr_ithinuloc)
   sname = 'local_x_index_at_thin_dams_at_u_nodes'
   fname = 'ithinuloc'
   lname = 'Local X-index of U-thin dams' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1) = numthinuloc
CASE (iarr_ithinv)
   sname = 'global_x_index_at_thin_dams_at_v_nodes'
   fname = 'ithinv'
   lname = 'Global X-index of V-thin dams' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = numthinv
CASE (iarr_ithinvloc)
   sname = 'local_x_index_at_thin_dams_at_v_nodes'
   fname = 'ithinvloc'
   lname = 'Local X-index of V-thin dams' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1) = numthinvloc
CASE (iarr_iwbaru)
   sname = 'global_x_index_at_weirs_barriers_at_u_nodes'
   fname = 'iwbaru'
   lname = 'Global X-index of U-node weirs/barriers' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_iwbaruloc)
   sname = 'local_x_index_at_weirs_barriers_at_u_nodes'
   fname = 'iwbaruloc'
   lname = 'Local X-index of U-node weirs/barriers' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1) = numwbaruloc
CASE (iarr_iwbarv)
   sname = 'global_x_index_at_weirs_barriers_at_v_nodes'
   fname = 'iwbarv'
   lname = 'Global X-index of V-node weirs/barriers' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_iwbarvloc)
   sname = 'local_x_index_at_weirs_barriers_at_v_nodes'
   fname = 'iwbarvloc'
   lname = 'Local X-index of V-node weirs/barriers' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1) = numwbarvloc
CASE (iarr_jdry)
   sname = 'global_y_index_at_dry_cells'
   fname = 'jdry'
   lname = 'Global Y-index of dry cells' 
   dtype = int_type
   nodim = 1
   ngdims(1) = numdry
CASE (iarr_jthinu)
   sname = 'global_y_index_at_thin_dams_at_u_nodes'
   fname = 'jthinu'
   lname = 'Global Y-index of U-thin dams' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = numthinu
CASE (iarr_jthinuloc)
   sname = 'local_y_index_at_thin_dams_at_u_nodes'
   fname = 'jthinuloc'
   lname = 'Local Y-index of U-thin dams' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1) = numthinuloc
CASE (iarr_jthinv)
   sname = 'global_y_index_at_thin_dams_at_v_nodes'
   fname = 'jthinv'
   lname = 'Global Y-index of V-thin dams' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = numthinv
CASE (iarr_jthinvloc)
   sname = 'local_y_index_at_thin_dams_at_v_nodes'
   fname = 'jthinvloc'
   lname = 'Local Y-index of V-thin dams' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1) = numthinvloc
CASE (iarr_jwbaru)
   sname = 'global_y_index_at_weirs_barriers_at_u_nodes'
   fname = 'jwbaru'
   lname = 'Global Y-index of U-weirs/barriers' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_jwbaruloc)
   sname = 'local_y_index_at_weirs_barriers_at_u_nodes'
   fname = 'jwbaruloc'
   lname = 'Global Y-index of U-node weirs/barriers' 
   dtype = int_type
   cnode = 'U'
   nodim = 1
   nldims(1) = numwbaruloc
CASE (iarr_jwbarv)
   sname = 'global_y_index_at_weirs_barriers_at_v_nodes'
   fname = 'jwbarv'
   lname = 'Global Y-index of V-node weirs/barriers' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_jwbarvloc)
   sname = 'local_y_index_at_weirs_barriers_at_v_nodes'
   fname = 'jwbarvloc'
   lname = 'Local Y-index of V-node weirs/barriers' 
   dtype = int_type
   cnode = 'V'
   nodim = 1
   nldims(1) = numwbarvloc
CASE (iarr_nowbaruprocs)
   sname = 'number_of_local_u_weirs_barriers_per_process' 
   fname = 'nowbaruprocs'
   lname = 'Number of local U-node weirs/barriers'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_nowbarvprocs)
   sname = 'number_of_local_v_weirs_barriers_per_process' 
   fname = 'nowbarvprocs'
   lname = 'Number of local V-node weirs/barriers'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_oricoefu)
   sname = 'discharge_coefficient_at_u_node_orifices'
   fname = 'oricoefu'
   lname = 'Discharge coefficient for orifices at U-nodes'
   unit = 'm1/2 s-1'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_oricoefv)
   sname = 'discharge_coefficient_at_v_node_orifices'
   fname = 'oricoefv'
   lname = 'Discharge coefficient for orifices at V-nodes'
   unit = 'm1/2 s-1'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_oriheightu)
   sname = 'heigth_of_u_node_orifices_wtr_sea_bed'
   fname = 'oriheightu'
   lname = 'Orifice height at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_oriheightv)
   sname = 'heigth_of_v_node_orifices_wtr_sea_bed'
   fname = 'oriheightv'
   lname = 'Orifice height at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv   
CASE (iarr_oriwidthu)
   sname = 'width_of_u_node_orifices'
   fname = 'oriwidthu'
   lname = 'Orifice width at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_oriwidthv)
   sname = 'width_of_u_node_orifices'
   fname = 'oriwidthv'
   lname = 'Orifice width at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_wbarcoefu)
   sname = 'discharge_coefficient_at_u_node_weirs'
   fname = 'wbarcoefu'
   lname = 'Discharge coefficient for weirs/barriers at U-nodes'
   unit = 'm1/2 s-1'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_wbarcoefv)
   sname = 'discharge_coefficient_at_v_node_weirs'
   fname = 'wbarcoefv'
   lname = 'Discharge coefficient for weirs/barriers at V-nodes'
   unit = 'm1/2 s-1'
   cnode = 'V'
   nodim = 1
   nldims(1) = numwbarv
CASE (iarr_wbarcrestu)
   sname = 'crest_heigth_wrt_sea_bed_at_u_node_weirs'
   fname = 'wbarcrestu'
   lname = 'Height of weir crest at U-nodes'
   unit = 'm'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_wbarcrestv)
   sname = 'crest_heigth_wrt_sea_bed_at_v_node_weirs'
   fname = 'wbarcrestv'
   lname = 'Height of weir crest at V-nodes'
   unit = 'm'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_wbarmodlu)
   sname = 'modular_limit_at_u_node_weirs'
   fname = 'wbarmodlu'
   lname = 'Modular limit at U-node weirs'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_wbarmodlv)
   sname = 'modular_limit_at_v_node_weirs'
   fname = 'wbarmodlv'
   lname = 'Modular limit at V-node weirs'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_wbarelossu)
   sname= 'energy_loss_at_u_node_weirs_barriers'
   fname = 'wbarelossu'
   lname = 'Energy loss sink term at U-node weirs/barriers'
   unit = 's-1'
   cnode = 'U'
   nodim = 1
   ngdims(1) = numwbaru
CASE (iarr_wbarelossv)
   sname= 'energy_loss_at_v_node_weirs_barriers'
   fname = 'wbarelossv'
   lname = 'Energy loss sink term at V-node weirs/barriers'
   unit = 's-1'
   cnode = 'V'
   nodim = 1
   ngdims(1) = numwbarv
CASE (iarr_mpvcov)
   sname= 'mpv_partial_surface_coverage'
   fname = 'mpvcov'
   lname = 'Fractional areal coverage of MPV structures'
   unit = '-'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc-1,nr-1/)
!
!3.15 Discharges
!---------------
!

CASE (iarr_disarea)
   sname = 'discharge_area'
   fname = 'disarea'
   lname = 'Discharge area'
   unit = 'm2'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_disdir)
   sname = 'discharge_direction'
   fname = 'disdir'
   lname = 'Discharge direction'
   unit = 'radian'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_disflag)
   sname = 'discharge_flag'
   fname = 'disflag'
   lname = 'Discharge flag'
   unit = 'log'
   dtype = log_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_dissal)
   sname = 'salinity_discharge'
   fname = 'dissal'
   lname = 'Salinity discharge'
   unit = 'PSU'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_disspeed)
   sname = 'discharge_speed'
   fname = 'disspeed'
   lname = 'Discharge speed'
   unit = 'm s-1'
   nodim = 1
   nldims(1) = numdisloc
CASE (iarr_distmp)
   sname = 'temperature_discharge'
   fname = 'distmp'
   lname = 'Temperature discharge'
   unit = 'degC'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_disvol)
   sname = 'volume_discharge_rate'
   fname = 'disvol'
   lname = 'Volume discharge rate'
   unit = 'm3 s-1'
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_idis)
   sname = 'global_x_index_at_discharge_location'
   fname = 'idis'
   lname = 'Global X-index of discharge location'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_idisloc)
   sname = 'local_x_index_at_discharge_locations'
   fname = 'idisloc'
   lname = 'Local X-index of discharge locations'
   dtype = int_type
   nodim = 1
   nldims(1) = numdisloc_ext
CASE (iarr_indexdisloc)
   sname = 'index_at_discharge_locations_for_local_to_global_mapping'
   fname = 'indexdisloc'
   lname = 'Index at discharge locations for local to global mapping'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_indexdisprocs)
   sname = 'index_at_discharge_locations_for_local_to_global_mapping_'//&
         & 'per_process'
   fname = 'indexdisprocs'
   lname = 'Index at discharge locations for local to global mapping '//&
         & 'per process'
   dtype = int_type
   nodim = 2
   ngdims(1:2) = (/numdis,nprocs/)
CASE (iarr_jdis)
   sname = 'global_y_index_at_discharge_locations'
   fname = 'jdis'
   lname = 'Global Y-index of discharge locations'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_jdisloc)
   sname = 'local_y_index_at_discharge_locations'
   fname = 'jdisloc'
   lname = 'Local Y-index of discharge locations'
   dtype = int_type
   nodim = 1
   nldims(1) = numdisloc_ext
CASE (iarr_kdis)
   sname = 'global_z_index_at_discharge_locations'
   fname = 'kdis'
   lname = 'Vertical grid index of discharge locations'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_kdistype)
   sname = 'type_of_vertical_discharge_location'
   fname = 'kdistype'
   lname = 'Type of vertical discharge location'
   dtype = int_type
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_nodisprocs)
   sname = 'number_of_local_discharge_locations_per_process' 
   fname = 'nodisprocs'
   lname = 'Number of local discharge locations perv process'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs
CASE (iarr_xdiscoord)
   sname = 'x_coordinate_of_discharge_locations'
   fname = 'xdiscoord'
   lname = 'X-coordinates of discharge locations'
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_ydiscoord)
   sname = 'y_coordinate_of_discharge_locations'
   fname = 'ydiscoord'
   lname = 'Y-coordinates of discharge locations'
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   nodim = 1
   ngdims(1) = numdis
CASE (iarr_zdiscoord)
   sname = 'z_coordinate_of_discharge_locations'
   fname = 'zdiscoord'
   lname = 'Vertical coordinates of discharge locations'
   unit = 'm'
   nodim = 1
   ngdims(1) = numdis

!
!3.16 Parameters for parallel mode
!---------------------------------
!

CASE (iarr_idprocs)
   sname = 'process_ids'
   fname = 'idprocs'
   lname = 'Process ids'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nprocs

!
!3.17 Energy equation, enstrophy, vorticity
!------------------------------------------
!

CASE (iarr_edens0d)
   sname = 'baroclinic_energy_in_sea_water'
   fname = 'edens0d'
   lname = 'Domain integrated baroclinic energy'
   unit = 'J'
   cnode = ''
   nodim = 0
CASE (iarr_edens2d)
   sname = 'integral_wrt_depth_of_baroclinic_energy__in_sea_water'
   fname = 'edens2d'
   lname = 'Vertically integrated baroclinic energy'
   unit = 'J m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_edens3d)
   sname = 'baroclinic_energy__in_sea_water'
   fname = 'edens3d'
   lname = 'Baroclinic energy'
   unit = 'J m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_edissip0d)
   sname = 'energy_dissipation_in_sea_water'
   fname = 'edissip0d'
   lname = 'Domain integrated energy dissipation'
   unit = 'W'
   cnode = ''
   nodim = 0
CASE (iarr_edissip2d)
   sname = 'integral_wrt_depth_of_energy_dissipation_in_sea_water'
   fname = 'edissip2d'
   lname = 'Vertically integrated energy dissipation'
   unit = 'W m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_edissip3d)
   sname = 'energy_dissipation_in_sea_water'
   fname = 'edissip3d'
   lname = 'Energy dissipation'
   unit = 'W m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_eflux2du)
   sname = 'integral_wrt_depth_of_x_energy_flux'
   fname = 'eflux2du'
   lname = 'X-component of depth-integrated energy flux'
   vname = 'Depth-integrated energy flux'
   unit = 'W m-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_eflux2dv)
   sname = 'integral_wrt_depth_of_y_energy_flux'
   fname = 'eflux2dv'
   lname = 'Y-component of depth-integrated energy flux'
   vname = 'Depth-integrated energy flux'
   unit = 'W m-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_eflux3du)
   sname = 'x_energy_flux' 
   fname = 'eflux3du'
   lname = 'X-component of energy flux'
   vname = 'Energy flux'
   unit = 'W m-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_eflux3dv)
   sname = 'y_energy_flux'
   fname = 'eflux3dv'
   lname = 'Y-component of energy flux'
   vname = 'Energy flux'
   unit = 'W m-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_eflux3dw)
   sname = 'upward_energy_flux'
   fname = 'eflux3dw'
   lname = 'Vertical component of energy flux'
   vname = 'Energy flux'
   unit = 'W m-2'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ekin0d)
   sname = 'kinetic_energy_in_sea_water'
   fname = 'ekin0d'
   lname = 'Domain integrated kinetic energy'
   unit = 'J'
   cnode = ''
   nodim = 0
CASE (iarr_ekin2d)
   sname = 'integral_wrt_depth_of_kinetic_energy_in_sea_water'
   fname = 'ekin2d'
   lname = 'Vertically integrated kinetic energy'
   unit = 'J m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ekin3d)
   sname = 'kinetic_energy_in_sea_water'
   fname = 'ekin3d'
   lname = 'Kinetic energy'
   unit = 'J m-3'
   cnode = ''
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_enstr0d)
   sname = 'enstrophy_in_sea_water'
   fname = 'enstr0d'
   lname = 'Domain integrated enstrophy'
   unit = 'm3 s-2'
   cnode = ''
   nodim = 0
CASE (iarr_epot0d)
   sname = 'potential_energy_in_sea_water'
   fname = 'epot0d'
   lname = 'Domain integrated potential energy'
   unit = 'J'
   cnode = ''
   nodim = 0
CASE (iarr_epot2d)
   sname = 'integral_wrt_depth_of_potential_energy_in_sea_water'
   fname = 'epot2d'
   lname = 'Vertically integrated potential energy'
   unit = 'J m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_etot0d)
   sname = 'total_energy_in_sea_water'
   fname = 'etot0d'
   lname = 'Domain integrated total energy'
   unit = 'J'
   cnode = ''
   nodim = 0
CASE (iarr_etot2d)
   sname = 'integral_wrt_depth_of_total_energy_in_sea_water'
   fname = 'etot2d'
   lname = 'Vertically integrated total energy'
   unit = 'J m-2'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_etot3d)
   sname = 'total_energy_density_in_sea_water'
   fname = 'etot3d'
   lname = 'Total energy'
   unit = 'J m-3'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_vortic2d)
   sname = 'integral_wrt_depth_of_sea_water_relative_vorticity'
   fname = 'vortic2d'
   lname = 'Vertically integrated vorticity'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_vortic3d)
   sname = 'sea_water_relative_vorticity'
   fname = 'vortic3d'
   lname = 'Vorticity'
   unit = 's-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

!
!3.18 Nesting
!------------
!

CASE (iarr_inst2dtype)
   sname = 'type_2d_open_boundary_data_output'
   fname = 'inst2dtype'
   lname = 'Type 2-D open boundary data output'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstatc)
   sname = 'index_of_first_local_nest_c_node_point_in_full_array'
   fname = 'lbhnstatc'
   lname = 'Index of first local nest C-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstatu)
   sname = 'index_of_first_local_nest_u_node_point_in_full_array'
   fname = 'lbhnstatu'
   lname = 'Index of first local nest U-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstatv)
   sname = 'index_of_first_local_nest_v_node_point_in_full_array'
   fname = 'lbhnstatv'
   lname = 'Index of first local nest V-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstatx)
   sname = 'index_of_first_local_nest_x_node_point_in_full_array'
   fname = 'lbhnstatx'
   lname = 'Index of first local nest X-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_lbhnstaty)
   sname = 'index_of_first_local_nest_y_node_point_in_full_array'
   fname = 'lbhnstaty'
   lname = 'Index of first local nest Y-node point in full array'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstatc)
   sname = 'local_number_of_nest_points_at_c_nodes_per_set'
   fname = 'nohnstatc'
   lname = 'Local number of nest points at C-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstatu)
   sname = 'local_number_of_nest_points_at_u_nodes_per_set'
   fname = 'nohnstatu'
   lname = 'Local number of nest points at U-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstatv)
   sname = 'local_number_of_nest_points_at_v_nodes_per_set'
   fname = 'nohnstatv'
   lname = 'Local number of nest points at V-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstatx)
   sname = 'local_number_of_nest_points_at_x_nodes_per_set'
   fname = 'nohnstatx'
   lname = 'Local number of nest points at X-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstaty)
   sname = 'local_number_of_nest_points_at_y_nodes_per_set'
   fname = 'nohnstaty'
   lname = 'Local number of nest points at Y-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstcprocs)
   sname = 'local_number_of_nest_points_at_c_nodes_per_process_and_set'
   fname = 'nohnstcprocs'
   lname = 'Local number of nest points at C-nodes per process and set'
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nonestsets/)
CASE (iarr_nohnstglbc)
   sname = 'number_of_nest_points_at_c_nodes_per_set'
   fname = 'nohnstglbc'
   lname = 'Number of nest points at C-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstglbu)
   sname = 'number_of_nest_points_at_u_nodes_per_set'
   fname = 'nohnstglbu'
   lname = 'Number of nest points at U-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstglbv)
   sname = 'number_of_nest_points_at_v_nodes_per_set'
   fname = 'nohnstglbv'
   lname = 'Number of nest points at V-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstglbx)
   sname = 'number_of_nest_points_at_x_nodes_per_set'
   fname = 'nohnstglbx'
   lname = 'Number of nest points at X-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstglby)
   sname = 'number_of_nest_points_at_y_nodes_per_set'
   fname = 'nohnstglby'
   lname = 'Number of nest points at Y-nodes per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets
CASE (iarr_nohnstuvprocs)
   sname = 'number_of_local_nest_points_at_u_and_v_nodes_per_process_and_set'
   fname = 'nohnstuvprocs'
   lname = 'Number of local nest points at U- and V-nodes per process and set'
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nonestsets/)
CASE (iarr_nohnstxyprocs)
   sname = 'number_of_local_nest_points_at_x_and_y_nodes_per_process_and_set'
   fname = 'nohnstuvprocs'
   lname = 'Number of local nest points at X- and Y-nodes per process and set'
   dtype = int_type
   cnode = ''
   nodim = 2
   ngdims(1:2) = (/nprocs,nonestsets/)
CASE (iarr_novnst)
   sname = 'number_of_vertical_nest_points_per_set'
   fname = 'novnst'
   lname = 'Number of vertical nest points per set'
   dtype = int_type
   cnode = ''
   nodim = 1
   ngdims(1) = nonestsets

!
!3.19 Elliptic parameters
!------------------------
!

CASE (iarr_ellac2d)
   sname = 'anticyclonic_2d_current'
   fname = 'ellac2d'
   lname = 'Anticyclonic 2-D current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellac3d)
   sname = 'anticyclonic_3d_current'
   fname = 'ellac3d'
   lname = 'Anticyclonic 3-D current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellcc2d)
   sname = 'cyclonic_2d_current'
   fname = 'ellcc2d'
   lname = 'Cyclonic 2-D current'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellcc3d)
   sname = 'cyclonic_3d_current'
   fname = 'ellcc3d'
   lname = 'Cyclonic 3-D current'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellinc2d)
   sname = 'inclination_of_2d_tidal_ellipses'
   fname = 'ellinc2d'
   lname = 'Inclination of 2-D tidal ellipses '
   unit = MERGE('degrees','radian ',DegreesOut)
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellinc3d)
   sname = 'inclination_of_3d_tidal_ellipses'
   fname = 'ellinc3d'
   lname = 'Inclination of 3-D tidal ellipses'
   unit = MERGE('degrees','radian ',DegreesOut)
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellip2d)
   sname = 'ellipticity_of_2d_tidal_ellipses'
   fname = 'ellip2d'
   lname = 'Ellipticity of 2-D tidal ellipses'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellip3d)
   sname = 'ellipticity_of_3d_tidal_ellipses'
   fname = 'ellip3d'
   lname = 'Ellipticity of 3-D tidal ellipses'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellmaj2d)
   sname = 'major_axis_of_2d_tidal_ellipses'
   fname = 'ellmaj2d'
   lname = 'Major axis of 2-D tidal ellipses'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellmaj3d)
   sname = 'major_axis_of_3d_tidal_ellipses'
   fname = 'ellmaj3d'
   lname = 'Major axis of 3-D tidal ellipses'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellmin2d)
   sname = 'minor_axis_of_2d_tidal_ellipses'
   fname = 'ellmin2d'
   lname = 'Minor axis of 2-D tidal ellipses'
   unit = 'm s-1'
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellmin3d)
   sname = 'minor_axis_of_3d_tidal_ellipses'
   fname = 'ellmin3d'
   lname = 'Minor axis of 3-D tidal ellipses'
   unit = 'm s-1'
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)
CASE (iarr_ellpha2d)
   sname = 'elliptic_phase_of_2d_tidal_ellipses'
   fname = 'ellpha2d'
   lname = 'Elliptic phase of 2-D tidal ellipses'
   unit = MERGE('degrees','radian ',DegreesOut)
   nodim = 2
   nldims(1:2) = (/ncloc,nrloc/)
   ngdims(1:2) = (/nc,nr/)
CASE (iarr_ellpha3d)
   sname = 'elliptic_phase_of_3d_tidal_ellipses'
   fname = 'ellpha3d'
   lname = 'Elliptic phase of 3-D tidal ellipses'
   unit = MERGE('degrees','radian ',DegreesOut)
   nodim = 3
   nldims(1:3) = (/ncloc,nrloc,nz/)
   ngdims(1:3) = (/nc,nr,nz/)

!
!3.20 Relative coordinates
!-------------------------
!
!---i-coordinate
CASE (iarr_icoordC)
   sname = 'i_interpolation_index_coordinate_from_external_to_c_node_grid'
   fname = 'icoordC'
   lname = 'I interpolation index coordinate from external to C-node grid'
   dtype = int_type
CASE (iarr_icoordCC)
   sname = 'i_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'c_node_grid' 
   fname = 'icoordCC'
   lname = 'I interpolation index coordinate from external C-node to model '//&
         & 'C-node grid'
   dtype = int_type
CASE (iarr_icoordCU)
   sname = 'i_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'u_node_grid'   
   fname = 'icoordCU'
   lname = 'I interpolation index coordinate from external C-node to model '//&
        & 'U-node grid'
   dtype = int_type
CASE (iarr_icoordCV)
   sname = 'i_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'v_node_grid'
   fname = 'icoordCV'
   lname = 'I interpolation index coordinate from external C-node to model '//&
         & 'V-node grid'
   dtype = int_type
CASE (iarr_icoordU)
   sname = 'i_interpolation_index_coordinate_from_external_to_u_node_grid'
   fname = 'icoordU'
   lname = 'I interpolation index coordinate from external to U-node grid'
   dtype = int_type
CASE (iarr_icoordUU)
   sname = 'i_interpolation_index_coordinate_from_external_u_node_to_model_'//&
         & 'u_node_grid'
   fname = 'icoordUU'
   lname = 'I interpolation index coordinate from external U-node to model '//&
         & 'U-node grid'
   dtype = int_type
CASE (iarr_icoordUY)
   sname = 'i_interpolation_index_coordinate_from_external_u_node_to_model_'//&
         & 'y_node_grid'
   fname = 'icoordUY'
   lname = 'I interpolation index coordinate from external U-node to model '//&
         & 'Y-node grid'
   dtype = int_type
CASE (iarr_icoordV)
   sname = 'i_interpolation_index_coordinate_from_external_to_v_node_grid'
   fname = 'icoordV'
   lname = 'I interpolation index coordinate from external to V-node grid'
   dtype = int_type
CASE (iarr_icoordVV)
   sname = 'i_interpolation_index_coordinate_from_external_v_node_to_model_'//&
         & 'v_node_grid'
   fname = 'icoordVV'
   lname = 'I interpolation index coordinate from external V-node to model '//&
         & 'V-node grid'
   dtype = int_type
CASE (iarr_icoordVX)
   sname = 'i_interpolation_index_coordinate_from_external_v_node_to_model_'//&
         & 'x_node_grid'
   fname = 'icoordVX'
   lname = 'I interpolation index coordinate from external V-node to model '//&
         & 'X-node grid'
   dtype = int_type

!---j-coordinate
CASE (iarr_jcoordC)
   sname = 'j_interpolation_index_coordinate_from_external_to_c_node_grid'
   fname = 'jcoordC'
   lname = 'J interpolation index coordinate from external to C-node grid'
   dtype = int_type
CASE (iarr_jcoordCC)
   sname = 'j_interpolation_index_coordinate_from_external_c_node_to_model_'//&
        & 'c_node_grid'
   fname = 'jcoordCC'
   lname = 'J interpolation index coordinate from external C-node to model '//&
         & 'C-node grid'
   dtype = int_type
CASE (iarr_jcoordCU)
   sname = 'j_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'u_node_grid'
   fname = 'jcoordCU'
   lname = 'J interpolation index coordinate from external C-node to model '//&
         & 'U-node grid'
   dtype = int_type
CASE (iarr_jcoordCV)
   sname = 'j_interpolation_index_coordinate_from_external_c_node_to_model_'//&
         & 'v_node_grid'
   fname = 'jcoordCV'
   lname = 'J interpolation index coordinate from external C-node to model '//&
         & 'V-node grid'
   dtype = int_type
CASE (iarr_jcoordU)
   sname = 'j_interpolation_index_coordinate_from_external_to_u_node_grid'
   fname = 'jcoordU'
   lname = 'J interpolation index coordinate from external to U-node grid'
   dtype = int_type
CASE (iarr_jcoordUU)
   sname = 'j_interpolation_index_coordinate_from_external_u_node_to_model_'//&
         & 'u_node_grid'
   fname = 'jcoordUU'
   lname = 'J interpolation index coordinate from external U-node to model '//&
         & 'U-node grid'
   dtype = int_type
CASE (iarr_jcoordUY)
   sname = 'j_interpolation_index_coordinate_from_external_u_node_to_model_'//&
         & 'y_node_grid'
   fname = 'jcoordUY'
   lname = 'J interpolation index coordinate from external U-node to model '//&
         & 'Y-node grid'
   dtype = int_type
CASE (iarr_jcoordV)
   sname = 'j_interpolation_index_coordinate_from_external_to_v_node_grid'
   fname = 'jcoordV'
   lname = 'J interpolation index coordinate from external to V-node grid'
   dtype = int_type
CASE (iarr_jcoordVV)
   sname = 'j_interpolation_index_coordinate_from_external_v_node_to_model_'//&
         & 'v_node_grid'
   fname = 'jcoordVV'
   lname = 'J interpolation index coordinate from external V-node to model '//&
         & 'V-node grid'
   dtype = int_type
CASE (iarr_jcoordVX)
   sname = 'j_interpolation_index_coordinate_from_external_v_node_to_model_'//&
         & 'x_node_grid'
   fname = 'jcoordVX'
   lname = 'J interpolation index coordinate from external V-node to model '//&
         & 'X-node grid'
   dtype = int_type

!---x-coordinate
CASE (iarr_weightsC)
   sname = 'weights_for_interpolation_from_external_to_c_node_grid'
   fname = 'weightsC'
   lname = 'Weights for interpolation from external to C-node grid'
CASE (iarr_weightsCC)
   sname = 'weights_for_interpolation_from_external_c_node_to_model_c_node_grid'
   fname = 'weightsCC'
   lname = 'Weights for interpolation from external C-node to model C-node grid'
CASE (iarr_weightsCU)
   sname = 'weights_for_interpolation_from_external_c_node_to_model_u_node_grid'
   fname = 'weightsCU'
   lname = 'Weights for interpolation from external C-node to model U-node grid'
CASE (iarr_weightsCV)
   sname = 'weights_for_interpolation_from_external_c_node_to_model_v_node_grid'
   fname = 'weightsCV'
   lname = 'Weights for interpolation from external C-node to model V-node grid'
CASE (iarr_weightsU)
   sname = 'weights_for_interpolation_from_external_to_u_node_grid'
   fname = 'weightsU'
   lname = 'Weights for interpolation from external to U-node grid'
CASE (iarr_weightsUU)
   sname = 'weights_for_interpolation_from_external_u_node_to_model_u_node_grid'
   lname = 'Weights for interpolation from external U-node to model U-node grid'
   fname = 'weightsUU'
CASE (iarr_weightsUY)
   sname = 'weights_for_interpolation_from_external_u_node_to_model_y_node_grid'
   fname = 'weightsUY'
   lname = 'Weights for interpolation from external U-node to model Y-node grid'
CASE (iarr_weightsV)
   sname = 'weights_for_interpolation_from_external_to_v_node_grid'
   fname = 'weightsV'
   lname = 'Weights for interpolation from external to V-node grid'
CASE (iarr_weightsVV)
   sname = 'weights_for_interpolation_from_external_v_node_to_model_v_node_grid'
   fname = 'weightsVV'
   lname = 'Weights for interpolation from external V-node to model V-node grid'
CASE (iarr_weightsVX)
   sname = 'weights_for_interpolation_from_external_v_node_to_model_x_node_grid'
   fname = 'weightsVX'
   lname = 'Weights for interpolation from external V-node to model X-node grid'

!
!3.21 Data coordinates
!---------------------
!

CASE (iarr_depmean)
   sname = 'mean_sea_floor_depth_below_sea_level'
   fname = 'depth'
   lname = 'Mean water depth'
   unit = 'm'
   cnode = ''
CASE (iarr_scoord)
   IF (iopt_grid_vtype.EQ.1) THEN
      fname = 'lev'
      sname = 'ocean_sigma_coordinate'
      lname = 'Sigma coordinate'
      nodim = 1
   ELSEIF (iopt_grid_vtype.EQ.2) THEN
      fname = 'lev'
      sname = 'ocean_s_coordinate'
      lname = 's-coordinate'
      nodim = 1
   ELSEIF (iopt_grid_vtype.EQ.3) THEN
      fname = 'glev'
      sname = 'ocean_s_coordinate_g'
      lname = 'General s-coordinate'
      nodim = 3
   ENDIF
CASE (iarr_scoordatc)
   sname = 'ocean_s_coordinate_at_c_nodes'
   fname = 'scoordatc'
   lname = 's-coordinate at C-nodes'
CASE (iarr_scoordatu)
   sname = 'ocean_s_coordinate_at_u_nodes'
   fname = 'scoordatu'
   lname = 's-coordinate at U-nodes'
   unit = 'm'
   cnode = 'U'
CASE (iarr_scoordatv)
   sname = 'ocean_s_coordinate_at_v_nodes'
   fname = 'scoordatv'
   lname = 's-coordinate at V-nodes'
   unit = 'm'
   cnode = 'V'
CASE (iarr_scoordatx)
   sname = 'ocean_s_coordinate_at_x_nodes'
   fname = 'scoordatx'
   lname = 's-coordinate at X-nodes'
   unit = 'm'
   cnode = 'X'
CASE (iarr_scoordaty)
   sname = 'ocean_s_coordinate_at_y_nodes'
   fname = 'scoordaty'
   lname = 's-coordinate at Y-nodes'
   unit = 'm'
   cnode = 'Y'
CASE (iarr_seamask)
   sname = 'sea_point_indices'
   fname = 'seapoint'
   lname = 'Sea point grid indices'
   dtype = int_type
CASE (iarr_statnames)
   sname = 'station_names'
   fname = 'statnames'
   lname = 'Station names'
   dtype = char_type
CASE (iarr_subname)
   sname = 'sub_name'
   fname = 'subname'
   lname = 'Number of sub-variable names'
   dtype = char_type
   cnode = ''
CASE (iarr_time)
   sname = 'time'
   fname = 'time'
   lname = 'Time'
   unit = 'date/time'
   dtype = char_type
   cnode = ''
CASE (iarr_xcoord)
   sname =  MERGE('projection_x_coordinate','longitude              ',&
                & iopt_grid_sph.EQ.0)
   fname =  MERGE('x  ','lon',iopt_grid_sph.EQ.0)
   lname =  MERGE('X-coordinate','longitude   ',iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = ''
CASE (iarr_xcoordatc)
   sname = MERGE('projection_x_coordinate_at_c_nodes',&
               & 'longitude                         ',iopt_grid_sph.EQ.0)
   fname =  'xcoordatc'
   lname =  MERGE('X-coordinate at C-nodes','longitude at C-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
CASE (iarr_xcoordatu)
   sname = MERGE('projection_x_coordinate_at_u_nodes',&
                & 'longitude_at_u_nodes              ',iopt_grid_sph.EQ.0)
   fname =  'xcoordatu'
   lname =  MERGE('X-coordinate at U-nodes','longitude at U-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'U'
CASE (iarr_xcoordatv)
   sname = MERGE('projection_x_coordinate_at_v_nodes',&
               & 'longitude_at_v_nodes              ',iopt_grid_sph.EQ.0)
   fname =  'xcoordatv'
   lname =  MERGE('X-coordinate at V-nodes','longitude at V-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'V'
CASE (iarr_xcoordatx)
   sname = MERGE('projection_x_coordinate_at_x_nodes',&
               & 'longitude_at_x_nodes              ',iopt_grid_sph.EQ.0)
   fname =  'xcoordatx'
   lname =  MERGE('X-coordinate at X-nodes','longitude at X-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'X'
CASE (iarr_xcoordaty)
   sname = MERGE('projection_x_coordinate_at_y_nodes',&
               & 'longitude_at_y_nodes              ',iopt_grid_sph.EQ.0)
   fname =  'xcoordaty'
   lname =  MERGE('X-coordinate at Y-nodes','longitude at Y-nodes   ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m           ','degrees_east',iopt_grid_sph.EQ.0)
   cnode = 'Y'
CASE (iarr_ycoord)
   sname =  MERGE('projection_y_coordinate','latitude               ',&
                & iopt_grid_sph.EQ.0)
   fname =  MERGE('y  ','lat',iopt_grid_sph.EQ.0)
   lname =  MERGE('Y-coordinate','latitude    ',iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = ''
CASE (iarr_ycoordatc)
   sname =  MERGE('projection_y_coordinate_at_c_nodes',&
                & 'latitude                          ',iopt_grid_sph.EQ.0)
   fname =  'ycoordatc'
   lname =  MERGE('Y-coordinate at C-nodes','latitude at C-nodes    ',&
                 & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
CASE (iarr_ycoordatu)
   sname =  MERGE('projection_y_coordinate_at_u_nodes',&
                & 'latitude_at_u_nodes               ',iopt_grid_sph.EQ.0)
   fname =  'ycoordatu'
   lname =  MERGE('Y-coordinate at U-nodes','latitude at U-nodes    ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'U'
CASE (iarr_ycoordatv)
   sname =  MERGE('projection_y_coordinate_at_v_nodes',&
                & 'latitude_at_v_nodes               ',iopt_grid_sph.EQ.0)
   fname =  'ycoordatv'
   lname =  MERGE('Y-coordinate at V-nodes','latitude at V-nodes    ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'V'
CASE (iarr_ycoordatx)
   sname =  MERGE('projection_y_coordinate_at_x_nodes',&
                & 'latitude_at_x_nodes               ',iopt_grid_sph.EQ.0)
   fname =  'ycoordatx'
   lname =  MERGE('Y-coordinate at X-nodes','latitude at X-nodes    ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'X'
CASE (iarr_ycoordaty)
   sname =  MERGE('projection_y_coordinate_at_Y_nodes',&
                & 'latitude_at_y_nodes               ',iopt_grid_sph.EQ.0)
   fname =  'ycoordaty'
   lname =  MERGE('Y-coordinate at Y-nodes','latitude at Y-nodes    ',&
                & iopt_grid_sph.EQ.0)
   unit =  MERGE ('m            ','degrees_north',iopt_grid_sph.EQ.0)
   cnode = 'Y'
CASE (iarr_zcoord)
   sname = 'ocean_z_coordinate'
   fname = 'z'
   lname = 'Z-coordinate'
   unit = 'm'
   cnode = ''
CASE (iarr_zcoordatc)
   sname = 'ocean_z_coordinate_at_c_nodes'
   fname = 'zcoordatc'
   lname = 'Z-coordinate at C-nodes'
   unit = 'm'
CASE (iarr_zcoordatu)
   sname = 'ocean_z_coordinate_at_u_nodes'
   fname = 'zcoordatu'
   lname = 'Z-coordinate at U-nodes'
   unit = 'm'
   cnode = 'U'
CASE (iarr_zcoordatv)
   sname = 'ocean_z_coordinate_at_v_nodes'
   fname = 'zcoordatv'
   lname = 'Z-coordinate at V-nodes'
   unit = 'm'
   cnode = 'V'
CASE (iarr_zcoordatx)
   sname = 'ocean_z_coordinate_at_x_nodes'
   fname = 'zcoordatx'
   lname = 'Z-coordinate at X-nodes'
   unit = 'm'
   cnode = 'X'
CASE (iarr_zcoordaty)
   sname = 'ocean_z_coordinate_at_y_nodes'
   fname = 'zcoordaty'
   lname = 'Z-coordinate at Y-nodes'
   unit = 'm'
   cnode = 'Y'
CASE (iarr_zcoordm)
   sname = 'ocean_z_coordinate_with_respect_to_mean_sea_level'
   fname = 'zmean'
   lname = 'Z-coordinate'
   unit = 'm'
   cnode = ''
CASE (iarr_zetout)
   sname = 'sea_surface_height_above_geoid'
   fname = 'zetout'
   lname = 'Surface elevation'
   unit = 'm'

!
!3.22 Model parameters
!---------------------
!

CASE (iarr_density_ref)
   sname = 'sea_water_reference_density '
   fname = 'density_ref'
   lname = 'Reference density'
   unit = 'kg m-3'
   nodim = 0
CASE (iarr_gacc_mean)
   sname = 'acceleration_due_to_gravity'
   fname = 'gacc_mean'
   lname = 'Mean acceleration of gravity'
   unit = 'm s-2'
   nodim = 0
END SELECT

!
!4. Formats
!----------
!

SELECT CASE (dtype)
   CASE (char_type); formatx = characterfmt
   CASE (int_type); formatx = integerfmt
   CASE (log_type); formatx = logicalfmt
   CASE (real_type); formatx = realfmt
   CASE (rlong_type); formatx = doublefmt
END SELECT
   
!
!5. Store
!--------
!

IF (PRESENT(f90_name)) f90_name = fname
IF (PRESENT(long_name)) long_name = lname
IF (PRESENT(standard_name)) standard_name = sname
IF (PRESENT(vector_name)) vector_name = vname
IF (PRESENT(units)) units = unit
IF (PRESENT(format)) format = formatx
IF (PRESENT(node)) node = cnode
IF (PRESENT(kind_type)) kind_type = dtype
IF (PRESENT(nrank)) nrank = nodim
IF (PRESENT(halo_dims)) halo_dims = nhdims
IF (PRESENT(global_dims)) global_dims = ngdims
IF (PRESENT(local_dims)) local_dims = nldims

IF (PRESENT(varatts)) THEN
   varatts%f90_name = fname
   varatts%long_name = lname
   varatts%standard_name = sname
   varatts%vector_name = vname
   varatts%units = unit
   varatts%format= formatx
   varatts%node = cnode   
   varatts%kind_type = dtype
   varatts%nrank = nodim
   varatts%halo_dims = nhdims
   varatts%global_dims = ngdims
   varatts%local_dims = nldims
   varatts%fill_value = float_fill
ENDIF


RETURN

END SUBROUTINE inquire_var

!========================================================================

FUNCTION inquire_varid(f90_name)
!************************************************************************
!
! *inquire_varid* Returns the key id of a variable given its name
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.12.1
!
! Description -
!
! Module calls - inquire_varid_bio, inquire_varid_sed, inquire_varid_part
!
!************************************************************************
!
USE modids
USE switches
USE biovars_routines, ONLY: inquire_varid_bio
USE partvars_routines, ONLY: inquire_varid_part
USE sedvars_routines, ONLY: inquire_varid_sed

!
!*Arguments
!
CHARACTER (LEN=*), INTENT(IN) :: f90_name
INTEGER :: inquire_varid

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*f90_name* CHAR    FORTRAN name of the variable
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER  :: ivarid


ivarid = 0

SELECT CASE(TRIM(f90_name))

!
!1. General and hydrodynamical variables
!---------------------------------------
!
!1.1 Model grid
!--------------
!

CASE ('alphatc_fld')
   ivarid = iarr_alphatc_fld
CASE ('alphatu_fld')
   ivarid = iarr_alphatu_fld
CASE ('alphatv_fld')
   ivarid = iarr_alphatv_fld
CASE ('coriolatu')
   ivarid = iarr_coriolatu
CASE ('coriolatv')
   ivarid = iarr_coriolatv
CASE ('delxatc')
   ivarid = iarr_delxatc
CASE ('delxatu')
   ivarid = iarr_delxatu
CASE ('delxatuv')
   ivarid = iarr_delxatuv
CASE ('delxatv')
   ivarid = iarr_delxatv
CASE ('delyatc')
   ivarid = iarr_delyatc
CASE ('delyatu')
   ivarid = iarr_delyatu
CASE ('delyatuv')
   ivarid = iarr_delyatuv
CASE ('delyatv')
   ivarid = iarr_delyatv
CASE ('delzatc')
   ivarid = iarr_delzatc
CASE ('delzatu')
   ivarid = iarr_delzatu
CASE ('delzatuv')
   ivarid = iarr_delzatuv
CASE ('delzatuw')
   ivarid = iarr_delzatuw
CASE ('delzatv')
   ivarid = iarr_delzatv
CASE ('delzatvw')
   ivarid = iarr_delzatvw
CASE ('delzatw')
   ivarid = iarr_delzatw
CASE ('dryfac')
   ivarid = iarr_dryfac
CASE ('gaccatc')
   ivarid = iarr_gaccatc
CASE ('gaccatu')
   ivarid = iarr_gaccatu
CASE ('gaccatv')
   ivarid = iarr_gaccatv
CASE ('gangleatc')
   ivarid = iarr_gangleatc
CASE ('garea')
   ivarid = iarr_garea
CASE ('gdelxglb')
   ivarid = iarr_gdelxglb
CASE ('gdelyglb')
   ivarid = iarr_gdelyglb
CASE ('gscoordatc')
   ivarid = iarr_gscoordatc
CASE ('gscoordatu')
   ivarid = iarr_gscoordatu
CASE ('gscoordatuvw')
   ivarid = iarr_gscoordatuvw
CASE ('gscoordatuw')
   ivarid = iarr_gscoordatuw
CASE ('gscoordatv')
   ivarid = iarr_gscoordatv
CASE ('gscoordatvw')
   ivarid = iarr_gscoordatvw
CASE ('gscoordatw')
   ivarid = iarr_gscoordatw
CASE ('gscoordglb')
   ivarid = iarr_gscoordglb
CASE ('gsigcoordatc')
   ivarid = iarr_gsigcoordatc
CASE ('gsigcoordatw')
   ivarid = iarr_gsigcoordatw
CASE ('gxcoord')
   ivarid = iarr_gxcoord
CASE ('gxcoordglb')
   ivarid = iarr_gxcoordglb
CASE ('gxlon')
   ivarid = iarr_gxlon
CASE ('gycoord')
   ivarid = iarr_gycoord
CASE ('gycoordglb')
   ivarid = iarr_gycoordglb
CASE ('gylat')
   ivarid = iarr_gylat
CASE ('indexobu')
   ivarid = iarr_indexobu
CASE ('indexobuprocs')
   ivarid = iarr_indexobuprocs
CASE ('indexobv')
   ivarid = iarr_indexobv
CASE ('indexobvprocs')
   ivarid = iarr_indexobvprocs
CASE ('indexobx')
   ivarid = iarr_indexobx
CASE ('indexobxprocs')
   ivarid = iarr_indexobxprocs
CASE ('indexoby')
   ivarid = iarr_indexoby
CASE ('indexobyprocs')
   ivarid = iarr_indexobyprocs
CASE ('iobu')
   ivarid = iarr_iobu
CASE ('iobuloc')
   ivarid = iarr_iobuloc
CASE ('iobv')
   ivarid = iarr_iobv
CASE ('iobvloc')
   ivarid = iarr_iobvloc
CASE ('iobx')
   ivarid = iarr_iobx
CASE ('iobxloc')
   ivarid = iarr_iobxloc
CASE ('ioby')
   ivarid = iarr_ioby
CASE ('iobyloc')
   ivarid = iarr_iobyloc
CASE ('jobu')
   ivarid = iarr_jobu
CASE ('jobuloc')
   ivarid = iarr_jobuloc
CASE ('jobv')
   ivarid = iarr_jobv
CASE ('jobvloc')
   ivarid = iarr_jobvloc
CASE ('jobx')
   ivarid = iarr_jobx
CASE ('jobxloc')
   ivarid = iarr_jobxloc
CASE ('joby')
   ivarid = iarr_joby
CASE ('jobyloc')
   ivarid = iarr_jobyloc
CASE ('maskatc_int')
   ivarid = iarr_maskatc_int
CASE ('mgvars_nc1procs')
   ivarid = iarr_mgvars_nc1procs
CASE ('mgvars_nc2procs')
   ivarid = iarr_mgvars_nc2procs
CASE ('mgvars_nr1procs')
   ivarid = iarr_mgvars_nr1procs
CASE ('mgvars_nr2procs')
   ivarid = iarr_mgvars_nr2procs
CASE ('ncprocs')
   ivarid = iarr_ncprocs
CASE ('nobuprocs')
   ivarid = iarr_nobuprocs
CASE ('nobvprocs')
   ivarid = iarr_nobvprocs
CASE ('nobxprocs')
   ivarid = iarr_nobxprocs
CASE ('nobyprocs')
   ivarid = iarr_nobyprocs
CASE ('nodeatc')
   ivarid = iarr_nodeatc
CASE ('nodeatu')
   ivarid = iarr_nodeatu
CASE ('nodeatuv')
   ivarid = iarr_nodeatuv
CASE ('nodeatuw')
   ivarid = iarr_nodeatuw
CASE ('nodeatv')
   ivarid = iarr_nodeatv
CASE ('nodeatvw')
   ivarid = iarr_nodeatvw
CASE ('node2du')
   ivarid = iarr_node2du
CASE ('node2duv')
   ivarid = iarr_node2duv
CASE ('node2dv')
   ivarid = iarr_node2dv
CASE ('nosbuprocs')
   ivarid = iarr_nosbuprocs
CASE ('nosbvprocs')
   ivarid = iarr_nosbvprocs
CASE ('nrprocs')
   ivarid = iarr_nrprocs
CASE ('nrvbuprocs')
   ivarid = iarr_nrvbuprocs
CASE ('nrvbvprocs')
   ivarid = iarr_nrvbvprocs
CASE ('rlxobcatu')
   ivarid = iarr_rlxobcatu
CASE ('rlxobcatv')
   ivarid = iarr_rlxobcatv
CASE ('seapoint')
   ivarid = iarr_seapoint
CASE ('seapointglb')
   ivarid = iarr_seapointglb
CASE ('soutobv')
   ivarid = iarr_soutobv
CASE ('soutoby')
   ivarid = iarr_soutoby
CASE ('westobu')
   ivarid = iarr_westobu
CASE ('westobx')
   ivarid = iarr_westobx

!
!1.2 Depths
!----------
!

CASE ('depmeanatc')
   ivarid = iarr_depmeanatc
CASE ('depmeanatu')
   ivarid = iarr_depmeanatu
CASE ('depmeanatuv')
   ivarid = iarr_depmeanatuv
CASE ('depmeanatv')
   ivarid = iarr_depmeanatv
CASE ('depmeanglb')
   ivarid = iarr_depmeanglb
CASE ('deptotatc')
   ivarid = iarr_deptotatc
CASE ('deptotatc_err')
   ivarid = iarr_deptotatc_err
CASE ('deptotatc_old')
   ivarid = iarr_deptotatc_old
CASE ('deptotatu')
   ivarid = iarr_deptotatu
CASE ('deptotatu_old')
   ivarid = iarr_deptotatu_old
CASE ('deptotatuv')
   ivarid = iarr_deptotatuv
CASE ('deptotatv')
   ivarid = iarr_deptotatv
CASE ('deptotatv_old')
   ivarid = iarr_deptotatv_old
CASE ('dzeta')
   ivarid = iarr_dzeta
CASE ('zeta')
   ivarid = iarr_zeta
CASE ('zeta_old')
   ivarid = iarr_zeta_old

!
!1.3 Currents
!------------
!

CASE ('hdvelmag')
   ivarid = iarr_hdvelmag
CASE ('hmvelmag')
   ivarid = iarr_hmvelmag
CASE ('hmvelmag_hadv_cfl')
   ivarid = iarr_hmvelmag_hadv_cfl
CASE ('hvelmag')
   ivarid = iarr_hvelmag
CASE ('hvelmag_hadv_cfl')
   ivarid = iarr_hvelmag_hadv_cfl
CASE ('p2dbcgradatu')
   ivarid = iarr_p2dbcgradatu
CASE ('p2dbcgradatv')
   ivarid = iarr_p2dbcgradatv
CASE ('p3dbcgradatu')
   ivarid = iarr_p3dbcgradatu
CASE ('p3dbcgradatv')
   ivarid = iarr_p3dbcgradatv
CASE ('udevint')
   ivarid = iarr_udevint
CASE ('udfvel')
   ivarid = iarr_udfvel
CASE ('udvel')
   ivarid = iarr_udvel
CASE ('udvel_old')
   ivarid = iarr_udvel_old
CASE ('ufvel')
   ivarid = iarr_ufvel
CASE ('umpred')
   ivarid = iarr_umpred
CASE ('umvel')
   ivarid = iarr_umvel
CASE ('umvel_hadv_cfl')
   ivarid = iarr_umvel_hadv_cfl
CASE ('umvel_old')
   ivarid = iarr_umvel_old
CASE ('uvel')
   ivarid = iarr_uvel
CASE ('uvel_hadv_cfl')
   ivarid = iarr_uvel_hadv_cfl
CASE ('uvel_old')
   ivarid = iarr_uvel_old
CASE ('vdevint')
   ivarid = iarr_vdevint
CASE ('vdfvel')
   ivarid = iarr_vdfvel
CASE ('vdvel')
   ivarid = iarr_vdvel
CASE ('vdvel_old')
   ivarid = iarr_vdvel_old
CASE ('vel2d')
   ivarid = iarr_vel2d
CASE ('vel3d')
   ivarid = iarr_vel3d
CASE ('vfvel')
   ivarid = iarr_vfvel
CASE ('vmpred')
   ivarid = iarr_vmpred
CASE ('vmvel')
   ivarid = iarr_vmvel
CASE ('vmvel_hadv_cfl')
   ivarid = iarr_vmvel_hadv_cfl
CASE ('vmvel_old')
   ivarid = iarr_vmvel_old
CASE ('vvel')
   ivarid = iarr_vvel
CASE ('vvel_hadv_cfl')
   ivarid = iarr_vvel_hadv_cfl
CASE ('vvel_old')   
   ivarid = iarr_vvel_old
CASE ('wphys')
   ivarid = iarr_wphys
CASE ('wvel')
   ivarid = iarr_wvel
CASE ('wvel_vadv_cfl')
   ivarid = iarr_wvel_vadv_cfl
   
!
!1.4 Density
!-----------
!

CASE ('beta_sal')
   ivarid = iarr_beta_sal
CASE ('beta_temp')
   ivarid = iarr_beta_temp
CASE ('dens')
   ivarid = iarr_dens
CASE ('sal')
   ivarid = iarr_sal
CASE ('temp')
   ivarid = iarr_temp

!
!1.5 Diffusion coefficients
!--------------------------
!

CASE ('hdifcoef2datc')
   ivarid = iarr_hdifcoef2datc
CASE ('hdifcoef2d_mom')
   ivarid = iarr_hdifcoef2d_mom
CASE ('hdifcoef2d_scal')
   ivarid = iarr_hdifcoef2d_scal
CASE ('hdifcoef2datuv')
   ivarid = iarr_hdifcoef2datuv
CASE ('hdifcoef3datc')
   ivarid = iarr_hdifcoef3datc
CASE ('hdifcoef3d_mom')
   ivarid = iarr_hdifcoef3d_mom
CASE ('hdifcoef3d_scal')
   ivarid = iarr_hdifcoef3d_mom
CASE ('hdifcoef3datu')
   ivarid = iarr_hdifcoef3datu
CASE ('hdifcoef3datuv')
   ivarid = iarr_hdifcoef3datuv
CASE ('hdifcoef3datv')
   ivarid = iarr_hdifcoef3datv
CASE ('hdifcoef3datw')
   ivarid = iarr_hdifcoef3datv
CASE ('kinvisc')
   ivarid = iarr_kinvisc
CASE ('mom_vdif_cfl')
   ivarid = iarr_mom_vdif_cfl
CASE ('scal_vdif_cfl')
   ivarid = iarr_scal_vdif_cfl
CASE ('vdifcoefmom')
   ivarid = iarr_vdifcoefmom
CASE ('vdifcoefscal')
   ivarid = iarr_vdifcoefscal
CASE ('vdifcoefscal_norot')
   ivarid = iarr_vdifcoefscal_norot
CASE ('vdifcoefscal_rot')
   ivarid = iarr_vdifcoefscal_rot
CASE ('vdifcoeftke')
   ivarid = iarr_vdifcoeftke
CASE ('xslopeatu_geo')
   ivarid = iarr_xslopeatu_geo
CASE ('xslopeatu_ziso')
   ivarid = iarr_xslopeatu_ziso
CASE ('xslopeatw_geo')
   ivarid = iarr_xslopeatw_geo
CASE ('yslopeatv_geo')
   ivarid = iarr_yslopeatv_geo
CASE ('yslopeatv_ziso')
   ivarid = iarr_yslopeatv_ziso
CASE ('yslopeatw_geo')
   ivarid = iarr_yslopeatw_geo
CASE ('2D_hdif_cfl')
   ivarid = iarr_2D_hdif_cfl

!
!1.6 Turbulence
!--------------
!

CASE ('buofreq2')
   ivarid = iarr_buofreq2
CASE ('dissip')
   ivarid = iarr_dissip
CASE ('ricnum')
   ivarid = iarr_ricnum
CASE ('shearfreq2')
   ivarid = iarr_shearfreq2
CASE ('tke')
   ivarid = iarr_tke
CASE ('tke_old')
   ivarid = iarr_tke_old
CASE ('tkezl')
   ivarid = iarr_tkezl
CASE ('zlmix')
   ivarid = iarr_zlmix

!
!1.7 Tides
!---------
!

CASE ('astro_earth')
   ivarid = iarr_astro_earth
CASE ('fnodal_anal')
   ivarid = iarr_fnodal_anal
CASE ('fnodal_astro')
   ivarid = iarr_fnodal_astro
CASE ('fnodal_obc')
   ivarid = iarr_fnodal_obc
CASE ('fxastro')
   ivarid = iarr_fxastro
CASE ('fyastro')
   ivarid = iarr_fyastro
CASE ('index_astro')
   ivarid = iarr_index_astro
CASE ('index_obc')
   ivarid = iarr_index_obc
CASE ('ispec_tides')
   ivarid = iarr_ispec_tides
CASE ('phase_anal')
   ivarid = iarr_phase_anal
CASE ('phase_astro')
   ivarid = iarr_phase_astro
CASE ('phase_obc')
   ivarid = iarr_phase_obc
CASE ('tidal_spectrum')
   ivarid = iarr_tidal_spectrum

!
!1.8 Meteo forcing
!-----------------
!

CASE ('airtemp')
   ivarid = iarr_airtemp
CASE ('atmpres')
   ivarid = iarr_atmpres
CASE ('cloud_cover')
   ivarid = iarr_cloud_cover
CASE ('evapminprec')
   ivarid = iarr_evapminprec
CASE ('evaporation')
   ivarid = iarr_evaporation
CASE ('meteodata')
   ivarid = iarr_meteodata
CASE ('precipitation')
   ivarid = iarr_precipitation
CASE ('relhum')
   ivarid = iarr_relhum
CASE ('sst')
   ivarid = iarr_sst
CASE ('uwindatc')
   ivarid = iarr_uwindatc
CASE ('uwindatc_old')
   ivarid = iarr_uwindatc_old
CASE ('vappres_air')
   ivarid = iarr_vappres_air
CASE ('vwindatc')
   ivarid = iarr_vwindatc
CASE ('vwindatc_old')
   ivarid = iarr_vwindatc_old
CASE ('windatc')
   ivarid = iarr_windatc

!
!1.9 Waves
!---------
!

CASE ('gxcoordglbwav')
   ivarid = iarr_gxcoordglbwav
CASE ('gycoordglbwav')
   ivarid = iarr_gycoordglbwav
CASE ('hbwdissipmag')
   ivarid = iarr_hbwdissipmag
CASE ('hmbwdissipmag')
   ivarid = iarr_hmbwdissipmag
CASE ('hmswdissipmag')
   ivarid = iarr_hmswdissipmag
CASE ('hswdissipmag')
   ivarid = iarr_hswdissipmag
CASE ('maskglbwav')
   ivarid = iarr_maskglbwav
CASE ('ubwdissipatc')
   ivarid = iarr_ubwdissipatc
CASE ('umbwdissipatc')
   ivarid = iarr_umbwdissipatc
CASE ('umswdissipatc')
   ivarid = iarr_umswdissipatc
CASE ('uswdissipatc')
   ivarid = iarr_uswdissipatc
CASE ('vbwdissipatc')
   ivarid = iarr_vbwdissipatc
CASE ('vmbwdissipatc')
   ivarid = iarr_vmbwdissipatc
CASE ('vmswdissipatc')
   ivarid = iarr_vmswdissipatc
CASE ('vswdissipatc')
   ivarid = iarr_vswdissipatc
CASE ('wavedir')
   ivarid = iarr_wavedir
CASE ('waveexcurs')
   ivarid = iarr_waveexcurs
CASE ('wavefreq')
   ivarid = iarr_wavefreq
CASE ('waveheight')
   ivarid = iarr_waveheight
CASE ('wavenum')
   ivarid = iarr_wavenum
CASE ('waveperiod')
   ivarid = iarr_waveperiod
CASE ('wavepres')
   ivarid = iarr_wavepres
CASE ('wavevel')
   ivarid = iarr_wavevel

!
!1.10 Stokes velocities
!----------------------
!

CASE ('hmstokesmag')
   ivarid = iarr_hmstokesmag
CASE ('hmveltotmag')
   ivarid = iarr_hmveltotmag
CASE ('hstokesmag')
   ivarid = iarr_hstokesmag
CASE ('hveltotmag')
   ivarid = iarr_hveltotmag
CASE ('stokessource2du')
   ivarid = iarr_stokessource2du
CASE ('stokessource2dv')
   ivarid = iarr_stokessource2dv
CASE ('udstokesatu')
   ivarid = iarr_udstokesatu
CASE ('umstokesatc')
   ivarid = iarr_umstokesatc
CASE ('umstokesatu')
   ivarid = iarr_umstokesatu
CASE ('umveltot')
   ivarid = iarr_umveltot
CASE ('ustokesatc')
   ivarid = iarr_ustokesatc
CASE ('ustokesatu')
   ivarid = iarr_ustokesatu
CASE ('uveltot')
   ivarid = iarr_uveltot
CASE ('vdstokesatv')
   ivarid = iarr_vdstokesatv
CASE ('vmstokesatc')
   ivarid = iarr_vmstokesatc
CASE ('vmstokesatv')
   ivarid = iarr_vmstokesatv
CASE ('vmveltot')
   ivarid = iarr_vmveltot
CASE ('vstokesatc')
   ivarid = iarr_vstokesatc
CASE ('vstokesatv')
   ivarid = iarr_vstokesatv
CASE ('vveltot')
   ivarid = iarr_vveltot
CASE ('wstokesatw')
   ivarid = iarr_wstokesatw

!
!1.11 Optical arrays
!-------------------
!

CASE ('optattcoef2')
   ivarid = iarr_optattcoef2
CASE ('qrad')
   ivarid = iarr_qrad
CASE ('radiance')
   ivarid = iarr_radiance

!
!1.12 Bottom/surface fluxes
!--------------------------
!

CASE ('bdragcoefatc')
   ivarid = iarr_bdragcoefatc
CASE ('bfricatu')
   ivarid = iarr_bfricatu
CASE ('bfricatv')
   ivarid = iarr_bfricatv
CASE ('bfricvel')
   ivarid = iarr_bfricvel
CASE ('bfricvel_max')
   ivarid = iarr_bfricvel_max
CASE ('bfricvel_wav')
   ivarid = iarr_bfricvel_wav
CASE ('bstresatc')
   ivarid = iarr_bstresatc
CASE ('bstresatc_max')
   ivarid = iarr_bstresatc_max
CASE ('bstresatc_wav')
   ivarid = iarr_bstresatc_wav
CASE ('bstresatu')
   ivarid = iarr_bstresatu
CASE ('bstresatv')
   ivarid = iarr_bstresatv
CASE ('cds')
   ivarid = iarr_cds
CASE ('ces')
   ivarid = iarr_ces
CASE ('chs')
   ivarid = iarr_chs
CASE ('fwave')
   ivarid = iarr_fwave
CASE ('qlatent')
   ivarid = iarr_qlatent
CASE ('qlwave')
   ivarid = iarr_qlwave
CASE ('qnonsol')
   ivarid = iarr_qnonsol
CASE ('qsensible')
   ivarid = iarr_qsensible
CASE ('qtot')
   ivarid = iarr_qtot
CASE ('sfricatc')
   ivarid = iarr_sfricatc
CASE ('ssalflux')
   ivarid = iarr_ssalflux
CASE ('sstresatc')
   ivarid = iarr_sstresatc
CASE ('surdat1d')
   ivarid = iarr_surdat1d
CASE ('ubstresatc')
   ivarid = iarr_ubstresatc
CASE ('ubstresatu')
   ivarid = iarr_ubstresatu
CASE ('usstresatc')
   ivarid = iarr_usstresatc
CASE ('usstresatu')
   ivarid = iarr_usstresatu
CASE ('vbstresatc')
   ivarid = iarr_vbstresatc
CASE ('vbstresatv')
   ivarid = iarr_vbstresatv
CASE ('vsstresatc')
   ivarid = iarr_vsstresatc
CASE ('vsstresatv')
   ivarid = iarr_vsstresatv
CASE ('wavedata')
   ivarid = iarr_wavedata
CASE ('wavethickatc')
   ivarid = iarr_wavethickatc
CASE ('zaroughatc')
   ivarid = iarr_zaroughatc
CASE ('zroughatc')
   ivarid = iarr_zroughatc

!
!1.13 Open boundary forcing
!--------------------------
!

CASE ('floutobu')
   ivarid = iarr_floutobu
CASE ('floutobv')
   ivarid = iarr_floutobv
CASE ('gxslope')
   ivarid = iarr_gxslope
CASE ('gxslope_amp')
   ivarid = iarr_gxslope_amp
CASE ('gxslope_pha')
   ivarid = iarr_gxslope_pha
CASE ('gyslope')
   ivarid = iarr_gyslope
CASE ('gyslope_amp')
   ivarid = iarr_gyslope_amp
CASE ('gyslope_pha')
   ivarid = iarr_gyslope_pha
CASE ('iloczobu')
   ivarid = iarr_iloczobu
CASE ('iloczobv')
   ivarid = iarr_iloczobv
CASE ('indexprof')
   ivarid = iarr_indexprof
CASE ('index2dobuv')
   ivarid = iarr_index2dobuv
CASE ('index2dobxy')
   ivarid = iarr_index2dobxy
CASE ('iobc2dtype')
   ivarid = iarr_iobc2dtype
CASE ('iprofobu')
   ivarid = iarr_iprofobu
CASE ('iprofobv')
   ivarid = iarr_iprofobv
CASE ('iprofobx')
   ivarid = iarr_iprofobx
CASE ('iprofoby')
   ivarid = iarr_iprofoby
CASE ('iqsecobu')
   ivarid = iarr_iqsecobu
CASE ('itypobu')
   ivarid = iarr_itypobu
CASE ('itypobv')
   ivarid = iarr_itypobv
CASE ('itypobx')
   ivarid = iarr_itypobx
CASE ('itypoby')
   ivarid = iarr_itypoby
CASE ('ityp2dobu')
   ivarid = iarr_ityp2dobu
CASE ('ityp2dobv')
   ivarid = iarr_ityp2dobv
CASE ('ityp2dobx')
   ivarid = iarr_ityp2dobx
CASE ('ityp2doby')
   ivarid = iarr_ityp2doby
CASE ('jqsecobv')
   ivarid = iarr_jqsecobv
CASE ('noprofsd')
   ivarid = iarr_noprofsd
CASE ('no2dobuv')
   ivarid = iarr_no2dobuv
CASE ('no2dobxy')
   ivarid = iarr_no2dobxy
CASE ('obcdatuv2d')
   ivarid = iarr_obcdatuv2d
CASE ('obcdatxy2d')
   ivarid = iarr_obcdatxy2d
CASE ('obcsalatu')
   ivarid = iarr_obcsalatu
CASE ('obcsalatv')
   ivarid = iarr_obcsalatv
CASE ('obctmpatu')
   ivarid = iarr_obctmpatu
CASE ('obctmpatv')
   ivarid = iarr_obctmpatv
CASE ('obc2uvatu')
   ivarid = iarr_obc2uvatu
CASE ('obc2uvatu_old')
   ivarid = iarr_obc2uvatu_old
CASE ('obc2uvatv')
   ivarid = iarr_obc2uvatv
CASE ('obc2uvatv_old')
   ivarid = iarr_obc2uvatv_old
CASE ('obc3uvatu')
   ivarid = iarr_obc3uvatu
CASE ('obc3uvatv')
   ivarid = iarr_obc3uvatv
CASE ('prof3dsal')
   ivarid = iarr_prof3dsal
CASE ('prof3dtmp')
   ivarid = iarr_prof3dtmp
CASE ('prof3dveluv')
   ivarid = iarr_prof3dveluv
CASE ('prof3dvelxy')
   ivarid = iarr_prof3dvelxy
CASE ('return_time')
   ivarid = iarr_return_time
CASE ('udatobu')
   ivarid = iarr_udatobu
CASE ('udatobu_amp')
   ivarid = iarr_udatobu_amp
CASE ('udatobu_pha')
   ivarid = iarr_udatobu_pha
CASE ('udatoby')
   ivarid = iarr_udatoby
CASE ('udatoby_amp')
   ivarid = iarr_udatoby_amp
CASE ('udatoby_pha')
   ivarid = iarr_udatoby_pha
CASE ('vdatobv')
   ivarid = iarr_vdatobv
CASE ('vdatobv_amp')
   ivarid = iarr_vdatobv_amp
CASE ('vdatobv_pha')
   ivarid = iarr_vdatobv_pha
CASE ('vdatobx')
   ivarid = iarr_vdatobx
CASE ('vdatobx_amp')
   ivarid = iarr_vdatobx_amp
CASE ('vdatobx_pha')
   ivarid = iarr_vdatobx_pha
CASE ('zdatobu')
   ivarid = iarr_zdatobu
CASE ('zdatobu_amp')
   ivarid = iarr_zdatobu_amp
CASE ('zdatobu_pha')
   ivarid = iarr_zdatobu_pha
CASE ('zdatobv')
   ivarid = iarr_zdatobv
CASE ('zdatobv_amp')
   ivarid = iarr_zdatobv_amp
CASE ('zdatobv_pha')
   ivarid = iarr_zdatobv_pha
CASE ('zeta_amp')
   ivarid = iarr_zeta_amp
CASE ('zeta_pha')
   ivarid = iarr_zeta_pha

!
!1.14 Structure module
!---------------------
!

CASE ('idry')
   ivarid = iarr_idry
CASE ('indexwbaru')
   ivarid = iarr_indexwbaru
CASE ('indexwbaruprocs')
   ivarid = iarr_indexwbaruprocs
CASE ('indexwbarv')
   ivarid = iarr_indexwbarv
CASE ('indexwbarvprocs')
   ivarid = iarr_indexwbarvprocs
CASE ('ithinu')
   ivarid = iarr_ithinu
CASE ('ithinuloc')
   ivarid = iarr_ithinuloc
CASE ('ithinv')
   ivarid = iarr_ithinv
CASE ('ithinvloc')
   ivarid = iarr_ithinvloc
CASE ('iwbaru')
   ivarid = iarr_iwbaru
CASE ('iwbaruloc')
   ivarid = iarr_iwbaruloc
CASE ('iwbarv')
   ivarid = iarr_iwbarv
CASE ('iwbarvloc')
   ivarid = iarr_iwbarvloc
CASE ('jdry')
   ivarid = iarr_jdry
CASE ('jthinu')
   ivarid = iarr_jthinu
CASE ('jthinuloc')
   ivarid = iarr_jthinuloc
CASE ('jthinv')
   ivarid = iarr_jthinv
CASE ('jthinvloc')
   ivarid = iarr_jthinvloc
CASE ('jwbaru')
   ivarid = iarr_jwbaru
CASE ('jwbaruloc')
   ivarid = iarr_jwbaruloc
CASE ('jwbarv')
   ivarid = iarr_jwbarv
CASE ('jwbarvloc')
   ivarid = iarr_jwbarvloc
CASE ('nowbaruprocs')
   ivarid = iarr_nowbaruprocs
CASE ('nowbarvprocs')
   ivarid = iarr_nowbarvprocs
CASE ('oricoefu')
   ivarid = iarr_oricoefu
CASE ('oricoefv')
   ivarid = iarr_oricoefv
CASE ('oriheightu')
   ivarid = iarr_oriheightu
CASE ('oriheightv')
   ivarid = iarr_oriheightv
CASE ('oriwidthu')
   ivarid = iarr_oriwidthu
CASE ('oriwidthv')
   ivarid = iarr_oriwidthv
CASE ('wbarcoefu')
   ivarid = iarr_wbarcoefu
CASE ('wbarcoefv')
   ivarid = iarr_wbarcoefv
CASE ('wbarcrestu')
   ivarid = iarr_wbarcrestu
CASE ('wbarcrestv')
   ivarid = iarr_wbarcrestv
CASE ('wbarmodlu')
   ivarid = iarr_wbarmodlu
CASE ('wbarmodlv')
   ivarid = iarr_wbarmodlv
CASE ('wbarelossu')
   ivarid = iarr_wbarelossu
CASE ('wbarelossv')
   ivarid = iarr_wbarelossv
CASE ('mpvcov')
   ivarid = iarr_mpvcov

!
!1.15 Discharges
!---------------
!

CASE ('disarea')
   ivarid = iarr_disarea
CASE ('disdir')
   ivarid = iarr_disdir
CASE ('disflag')
   ivarid = iarr_disflag
CASE ('dissal')
   ivarid = iarr_disspeed
CASE ('disspeed')
   ivarid = iarr_disspeed
CASE ('distmp')
   ivarid = iarr_distmp
CASE ('disvol')
   ivarid = iarr_disvol
CASE ('idis')
   ivarid = iarr_idis
CASE ('idisloc')
   ivarid = iarr_idisloc
CASE ('indexdisloc')
   ivarid = iarr_indexdisloc
CASE ('indexdisprocs')
   ivarid = iarr_indexdisprocs
CASE ('jdis')
   ivarid = iarr_jdis
CASE ('jdisloc')
   ivarid = iarr_jdisloc
CASE ('kdis')
   ivarid = iarr_kdis
CASE ('kdistype')
   ivarid = iarr_kdistype
CASE ('nodisprocs')
   ivarid = iarr_nodisprocs
CASE ('xdiscoord')
   ivarid = iarr_xdiscoord
CASE ('ydiscoord')
   ivarid = iarr_ydiscoord
CASE ('zdiscoord')
   ivarid = iarr_zdiscoord

!
!1.16 Parameters for parallel mode
!---------------------------------
!

CASE ('idprocs')
   ivarid = iarr_idprocs

!
!1.17 Energy equation, enstrophy, vorticity
!------------------------------------------
!

CASE ('edens0d')
   ivarid = iarr_edens0d
CASE ('edens2d')
   ivarid = iarr_edens2d
CASE ('edens3d')
   ivarid = iarr_edens3d
CASE ('edissip0d')
   ivarid = iarr_edissip0d
CASE ('edissip2d')
   ivarid = iarr_edissip2d
CASE ('edissip3d')
   ivarid = iarr_edissip3d
CASE ('eflux2du')
   ivarid = iarr_eflux2du
CASE ('eflux2dv')
   ivarid = iarr_eflux2dv
CASE ('eflux3du')
   ivarid = iarr_eflux3du
CASE ('eflux3dv')
   ivarid = iarr_eflux3dv
CASE ('eflux3dw')
   ivarid = iarr_eflux3dw
CASE ('ekin0d')
   ivarid = iarr_ekin0d
CASE ('ekin2d')
   ivarid = iarr_ekin2d
CASE ('ekin3d')
   ivarid = iarr_ekin3d
CASE ('enstr0d')
   ivarid = iarr_enstr0d
CASE ('epot0d')
   ivarid = iarr_epot0d
CASE ('epot2d')
   ivarid = iarr_epot2d
CASE ('etot0d')
   ivarid = iarr_etot0d
CASE ('etot2d')
   ivarid = iarr_etot2d
CASE ('etot3d')
   ivarid = iarr_etot3d
CASE ('vortic2d')
   ivarid = iarr_vortic2d
CASE ('vortic3d')
   ivarid = iarr_vortic3d

!
!1.18 Nesting
!------------
!

CASE ('inst2dtype')
   ivarid = iarr_inst2dtype
CASE ('lbhnstatc')
   ivarid = iarr_lbhnstatc
CASE ('lbhnstatu')
   ivarid = iarr_lbhnstatu
CASE ('lbhnstatv')
   ivarid = iarr_lbhnstatv
CASE ('lbhnstatx')
   ivarid = iarr_lbhnstatx
CASE ('lbhnstaty')
   ivarid = iarr_lbhnstaty
CASE ('nohnstatc')
   ivarid = iarr_nohnstatc
CASE ('nohnstatu')
   ivarid = iarr_nohnstatu
CASE ('nohnstatv')
   ivarid = iarr_nohnstatv
CASE ('nohnstatx')
   ivarid = iarr_nohnstatx
CASE ('nohnstaty')
   ivarid = iarr_nohnstaty
CASE ('nohnstcprocs')
   ivarid = iarr_nohnstcprocs
CASE ('nohnstglbc')
   ivarid = iarr_nohnstglbc
CASE ('nohnstglbu')
   ivarid = iarr_nohnstglbu
CASE ('nohnstglbv')
   ivarid = iarr_nohnstglbv
CASE ('nohnstglbx')
   ivarid = iarr_nohnstglbx
CASE ('nohnstglby')
   ivarid = iarr_nohnstglby
CASE ('nohnstuvprocs')
   ivarid = iarr_nohnstuvprocs
CASE ('nohnstxyprocs')
   ivarid = iarr_nohnstuvprocs
CASE ('novnst')
   ivarid = iarr_novnst

!
!1.19 Elliptic parameters
!------------------------
!

CASE ('ellac2d')
   ivarid = iarr_ellac2d
CASE ('ellac3d')
   ivarid = iarr_ellac3d
CASE ('ellcc2d')
   ivarid = iarr_ellcc2d
CASE ('ellcc3d')
   ivarid = iarr_ellcc3d
CASE ('ellinc2d')
   ivarid = iarr_ellinc2d
CASE ('ellinc3d')
   ivarid = iarr_ellinc3d
CASE ('ellip2d')
   ivarid = iarr_ellip2d
CASE ('ellip3d')
   ivarid = iarr_ellip3d
CASE ('ellmaj2d')
   ivarid = iarr_ellmaj2d
CASE ('ellmaj3d')
   ivarid = iarr_ellmaj3d
CASE ('ellmin2d')
   ivarid = iarr_ellmin2d
CASE ('ellmin3d')
   ivarid = iarr_ellmin3d
CASE ('ellpha2d')
   ivarid = iarr_ellpha2d
CASE ('ellpha3d')
   ivarid = iarr_ellpha3d

!
!1.20 Relative coordinates
!-------------------------
!
   
CASE ('icoordC')
   ivarid = iarr_icoordC
CASE ('icoordCC')
   ivarid = iarr_icoordCC
CASE ('icoordCU')
   ivarid = iarr_icoordCU
CASE ('icoordCV')
   ivarid = iarr_icoordCV
CASE ('icoordU')
   ivarid = iarr_icoordU
CASE ('icoordUU')
   ivarid = iarr_icoordUU
CASE ('icoordUY')
   ivarid = iarr_icoordUY
CASE ('icoordV')
   ivarid = iarr_icoordV
CASE ('icoordVV')
   ivarid = iarr_icoordVV
CASE ('icoordVX')
   ivarid = iarr_icoordVX
CASE ('jcoordC')
   ivarid = iarr_jcoordC
CASE ('jcoordCC')
   ivarid = iarr_jcoordCC
CASE ('jcoordCU')
   ivarid = iarr_jcoordCU
CASE ('jcoordCV')
   ivarid = iarr_jcoordCV
CASE ('jcoordU')
   ivarid = iarr_jcoordU
CASE ('jcoordUU')
   ivarid = iarr_jcoordUU
CASE ('jcoordUY')
   ivarid = iarr_jcoordUY
CASE ('jcoordV')
   ivarid = iarr_jcoordV
CASE ('jcoordVV')
   ivarid = iarr_jcoordVV
CASE ('jcoordVX')
   ivarid = iarr_jcoordVX
CASE ('weightsC')
   ivarid = iarr_weightsC
CASE ('weightsCC')
   ivarid = iarr_weightsCC
CASE ('weightsCU')
   ivarid = iarr_weightsCU
CASE ('weightsCV')
   ivarid = iarr_weightsCV
CASE ('weightsU')
   ivarid = iarr_weightsU
CASE ('weightsUU')
   ivarid = iarr_weightsUU
CASE ('weightsUY')
   ivarid = iarr_weightsUY
CASE ('weightsV')
   ivarid = iarr_weightsV
CASE ('weightsVV')
   ivarid = iarr_weightsVV
CASE ('weightsVX')
   ivarid = iarr_weightsVX

!
!1.21 Data coordinates
!---------------------
!

CASE ('depth')
   ivarid = iarr_depmean
CASE ('glev')
   ivarid = iarr_scoord
CASE ('lev')
   ivarid = iarr_scoord
CASE ('scoordatc')
   ivarid = iarr_scoordatc
CASE ('scoordatu')
   ivarid = iarr_scoordatu
CASE ('scoordatv')
   ivarid = iarr_scoordatv
CASE ('scoordatx')
   ivarid = iarr_scoordatx
CASE ('scoordaty')
   ivarid = iarr_scoordaty
CASE ('seamask')
   ivarid = iarr_seapoint
CASE ('statnames')
   ivarid = iarr_statnames
CASE ('subname')
   ivarid = iarr_subname
CASE ('time')
   ivarid = iarr_time
CASE ('x','lon')
   ivarid = iarr_xcoord
CASE ('xcoordatc')
   ivarid = iarr_xcoordatc
CASE ('xcoordatu')
   ivarid = iarr_xcoordatu
CASE ('xcoordatv')
   ivarid = iarr_xcoordatv
CASE ('xcoordatx')
   ivarid = iarr_xcoordatx
CASE ('xcoordaty')
   ivarid = iarr_xcoordaty
CASE ('y','lat')
   ivarid = iarr_ycoord
CASE ('ycoordatc')
   ivarid = iarr_ycoordatc
CASE ('ycoordatu')
   ivarid = iarr_ycoordatu
CASE ('ycoordatv')
   ivarid = iarr_ycoordatv
CASE ('ycoordatx')
   ivarid = iarr_ycoordatx
CASE ('ycoordaty')
   ivarid = iarr_ycoordaty
CASE ('z')
   ivarid = iarr_zcoord
CASE ('zcoordatc')
   ivarid = iarr_zcoordatc
CASE ('zcoordatu')
   ivarid = iarr_zcoordatu
CASE ('zcoordatv')
   ivarid = iarr_zcoordatv
CASE ('zcoordatx')
   ivarid = iarr_zcoordatx
CASE ('zcoordaty')
   ivarid = iarr_zcoordaty
CASE ('zmean')
   ivarid = iarr_zcoordm
CASE ('zetout')
   ivarid = iarr_zetout

!
!1.22 Model parameters
!---------------------
!

CASE ('density_ref')
   ivarid = iarr_density_ref
CASE ('gacc_mean')
   ivarid = iarr_gacc_mean
END SELECT

!
!2. Sediment and morphology
!--------------------------
!

IF (ivarid.EQ.0.AND.iopt_sed.GT.0) ivarid = inquire_varid_sed(f90_name)

!
!3. Biology
!----------
!

IF (ivarid.EQ.0.AND.iopt_biolgy.GT.0) ivarid = inquire_varid_bio(f90_name)

!
!4. Particle model
!-----------------
!

IF (ivarid.EQ.0.AND.iopt_part_model.GT.0) THEN
   ivarid = inquire_varid_part(f90_name)
ENDIF

!
!5. Warning if no match is found
!-------------------------------
!

IF (ivarid.EQ.0.AND.warnflag) THEN
   WRITE(iowarn,'(A)') 'WARNING: no key id exists for variable: '//&
                     & TRIM(f90_name)
ENDIF
   
inquire_varid = ivarid


RETURN

END FUNCTION inquire_varid

!========================================================================

SUBROUTINE set_modfiles_atts(iddesc,ifil,iotype,filepars)
!************************************************************************
!
! *set_modfiles_atts* Obtain the global attributes of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.x
!
! Description -
!
! Module calls - inquire_var, set_biofiles_atts, set_modfiles_name,
!                set_partfiles_atts, set_sedfiles_atts
!
!************************************************************************
!
USE gridpars
USE nestgrids
USE obconds
USE physpars
USE structures
USE switches
USE tide
USE biovars_routines, ONLY: set_biofiles_atts
USE sedvars_routines, ONLY: set_sedfiles_atts
USE partvars_routines, ONLY: set_partfiles_atts

!
!* Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, iotype
TYPE (FileParams), INTENT(INOUT) :: filepars

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*   INTEGER File id
!*ifil*     INTEGER File number
!*iotype*   INTEGER I/O type of file
!             = 1 => input
!             = 2 => output
!*filepars* DERIVED Attributes of forcing file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=12) :: cfil
CHARACTER (LEN=lendesc) :: title
INTEGER :: nhtype, nocoords, nofiles, nosubnames, novars, nvars, subnameid, &
         & timeid


pglev = pglev + 1
procname(pglev) = 'set_modfiles_atts'
IF (loglev1.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ENDIF
ENDIF

!
!1. Initialise parameters
!------------------------
!

nocoords = 0; nosubnames = 0
timeid = 0; subnameid = 0
nofiles = maxdatafiles(iddesc,1)
WRITE (cfil,'(I12)') ifil; cfil = ADJUSTL(cfil)

!
!2. Define attributes
!--------------------
!

SELECT CASE (iddesc)

!
!2.1  Domain decomposition
!------------------------
!

CASE (io_mppmod)
   title = 'Domain decomposition'
   novars = 4

!
!2.2 Initial conditions
!----------------------
!

CASE (io_inicon,io_fincon)

   SELECT CASE (ifil)
      CASE (ics_phys)
         title = 'Physical initial conditions'
         novars = 0
         IF (iopt_mode_2D.GT.0) novars = novars + 3
         IF (iopt_mode_3D.GT.0) novars = novars + 2
         IF (iopt_mode_3D.GT.0.AND.iopt_grid_nodim.EQ.3) novars = novars + 1
         IF (iopt_temp.GT.0) novars = novars + 1
         IF (iopt_sal.GT.0) novars = novars + 1
         IF (iopt_vdif_coef.EQ.3) THEN
            novars = novars + 1
            IF (iopt_turb_ntrans.GT.0.AND.iopt_turb_param.GT.0) THEN
               novars = novars + 1
            ENDIF
         ENDIF
         IF (iopt_bstres_form.EQ.2) THEN
            IF (iopt_bstres_drag.EQ.2.OR.iopt_bstres_drag.EQ.4) THEN
               novars = novars + 1
            ENDIF
         ENDIF
         IF (iopt_astro_pars.EQ.1.AND.nconobc.GT.0) novars = novars + 1
         IF (iopt_weibar.EQ.1) THEN
            IF (numwbaru.GT.0) novars = novars + 1
            IF (numwbarv.GT.0) novars = novars + 1
         ENDIF
         IF  (obcsalatu_init) novars = novars + 1
         IF  (obcsalatv_init) novars = novars + 1
         IF  (obctmpatu_init) novars = novars + 1
         IF  (obctmpatv_init) novars = novars + 1
         IF  (obc2uvatu_init) novars = novars + 1
         IF  (obc2uvatv_init) novars = novars + 1
         IF  (obc3uvatu_init) novars = novars + 1
         IF  (obc3uvatv_init) novars = novars + 1
         nocoords = 1
     CASE (ics_sed,ics_morph)
        CALL set_sedfiles_atts(iddesc,ifil,iotype,filepars)
        GOTO 1000
     CASE (ics_bio)
        CALL set_biofiles_atts(iddesc,ifil,iotype,filepars)
        GOTO 1000
     CASE (ics_part)
        CALL set_partfiles_atts(iddesc,ifil,iotype,filepars)
        GOTO 1000
   END SELECT

!
!2.3 Model grid and data locations
!---------------------------------
!

CASE (io_modgrd)
   title = 'Model grid'
   novars = MERGE(1,3,iopt_grid_htype.EQ.1)
   IF (iopt_grid_vtype.GT.1) novars = novars + 1
   IF (nobu.GT.0) novars = novars + 2
   IF (nobv.GT.0) novars = novars + 2

CASE (io_metabs)
   title = 'Absolute coordinates of meteo grid'
   novars = 2

CASE (io_sstabs)
   title = 'Absolute coordinates of SST grid'
   novars = 2

CASE (io_wavabs)
   title = 'Absolute coordinates of wave grid'
   novars = 2

CASE (io_nstgrd)
   title = 'Nested grid locations: '//TRIM(cfil)
   novars = 0
   IF (nohnstglbc(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbu(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbv(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglbx(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF
   IF (nohnstglby(ifil).GT.0) THEN
      novars = novars + 2
      IF (novnst(ifil).GT.0) novars = novars + 1
   ENDIF

CASE (io_metrel)
   title = 'Relative coordinates of meteo grid'
   novars = 3

CASE (io_sstrel)
   title = 'Relative coordinates of SST grid'
   novars = 3

CASE (io_wavrel)
   title = 'Relative coordinates of wave grid'
   novars = 3

!
!2.4 Open boundaries
!-------------------
!

CASE (io_1uvsur)
   IF (ifil.EQ.1) THEN
      title = 'Type of boundary conditions for 1-D mode'
      SELECT CASE (iopt_sbc_1D)
         CASE (1); novars = 6
         CASE (2); novars = 2
         CASE (3); novars = 4
      END SELECT
   ELSEIF (ifil.EQ.2) THEN
      title = 'Boundary data for 1-D mode'
      novars = 1
      nocoords = 1
   ENDIF

CASE (io_2uvobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for 2-D mode'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nconobc.GT.0) THEN
         IF (nobu.GT.0) novars = novars + 4
         IF (nobv.GT.0) novars = novars + 4
      ENDIF
      IF (nofiles.GT.1) novars = novars + 3
      IF (nqsecobu.GT.0) novars = novars + 1
      IF (nqsecobv.GT.0) novars = novars + 1
   ELSE
      title = 'Open boundary data for 2-D mode: '//TRIM(cfil)
      novars = 1
      nocoords = 1
   ENDIF

CASE (io_2xyobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of tangential open boundary conditions for 2-D mode'
      novars = 0
      IF (nobx.GT.0) THEN
         novars = novars + MERGE(3,1,nconobc.GT.0)
      ENDIF
      IF (noby.GT.0) THEN
         novars = novars + MERGE(3,1,nconobc.GT.0)
      ENDIF
      IF (nofiles.GT.1) novars = novars + 2
   ELSE
      title = 'Tangential open boundary data for 2-D mode: '//TRIM(cfil)
      novars = 1
      nocoords = 1
   ENDIF

CASE (io_3uvobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for 3-D current'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2
   ELSE
      title = 'Open boundary data for 3-D current: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

CASE (io_3xyobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of tangential open boundary conditions for 3-D current'
      novars = 0
      IF (nobx.GT.0) novars = novars + 2
      IF (noby.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2
   ELSE
      title = 'Tangential open boundary data for 3-D current: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

CASE (io_salobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for salinity'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2
   ELSE
      title = 'Open boundary data for salinity: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

CASE (io_tmpobc)
   IF (ifil.EQ.1) THEN
      title = 'Type of open boundary conditions for temperature'
      novars = 0
      IF (nobu.GT.0) novars = novars + 2
      IF (nobv.GT.0) novars = novars + 2
      IF (nofiles.GT.1) novars = novars + 2
   ELSE
      title = 'Open boundary data for temperature: '//TRIM(cfil)
      novars = 1; nocoords = 1
   ENDIF

!
!2.5 Nested output
!-----------------
!

CASE (io_nstspc)
   title = 'Number of nested grid locations'
   novars = 7

CASE (io_2uvnst)
   title = 'Nested boundary data for 2-D mode: '//TRIM(cfil)
   novars = 1
   nocoords = 1

CASE (io_3uvnst)
   title = 'Nested boundary data for 3-D current: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_2xynst)
   title = 'Nested boundary data for 2-D tangential mode: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_3xynst)
   title = 'Nested boundary data for 3-D tangential current: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_salnst)
   title = 'Nested boundary data for salinity: '//TRIM(cfil)
   novars = 1; nocoords = 1

CASE (io_tmpnst)
   title = 'Nested boundary data for temperature: '//TRIM(cfil)
   novars = 1; nocoords = 1

!
!2.6 Surface data
!----------------
!

CASE (io_metsur)
   title = 'Meteo input data'
   nhtype = surfacegrids(igrd_meteo,ifil)%nhtype
   nvars = 0
   IF (iopt_meteo_stres.EQ.1) nvars = nvars + 2
   IF (iopt_meteo_pres.EQ.1) nvars = nvars + 1
   IF (iopt_meteo_heat.EQ.1) THEN
      SELECT CASE (iopt_meteo_data)
         CASE (1); nvars = nvars + 3
         CASE (2); nvars = nvars + 2
      END SELECT
   ENDIF
   IF (iopt_meteo_precip.EQ.1) nvars = nvars + 1
   IF (nhtype.EQ.0) THEN
      novars = 1; nosubnames = nvars; subnameid = 1
   ELSE
      novars = nvars
   ENDIF
   nocoords = 1

CASE (io_sstsur)
   title = 'SST input data'
   novars = 1; nocoords = 1

CASE (io_wavsur)
   title = 'Surface wave input data'
   nhtype = MERGE(surfacegrids(igrd_waves,1)%nhtype,4,iopt_waves_couple.EQ.0)
   nvars = 3
   IF (iopt_waves_form.EQ.2) THEN
      nvars = nvars + 2
      IF (iopt_waves_curr.EQ.1) THEN
         nvars = nvars + 2
         IF (iopt_waves_pres.EQ.1) nvars = nvars + 1
      ENDIF
   ENDIF
   IF (iopt_waves_dissip.EQ.1) nvars = nvars + 4
   IF (nhtype.EQ.0) THEN
      novars = 1; nosubnames = nvars; subnameid = 1
   ELSE
      novars = nvars
   ENDIF
   nocoords = 1

!
!2.7 Structures
!--------------
!

CASE (io_drycel)
   title = 'Dry cells'
   novars = 2

CASE (io_thndam)
   title = 'Thin dams'
   novars = 0
   IF (numthinu.GT.0) novars = novars + 2
   IF (numthinv.GT.0) novars = novars + 2

CASE (io_weibar)
   title = 'Weirs/barriers'
   novars = 0
   IF (numwbaru.GT.0) novars = novars + 8
   IF (numwbarv.GT.0) novars = novars + 8

CASE (io_mpvcov)
   title = 'MPV'
   novars   = 1

!
!2.8 Discharges
!--------------
!

CASE (io_disspc)
   title = 'Discharge location type and flagging'
   novars = 1

CASE (io_disloc)
   title = 'Discharge locations'
   novars = 3; nocoords = 1

CASE (io_disvol)
   title = 'Volume discharge rates'
   novars = 1; nocoords = 1

CASE (io_discur)
   title = 'Discharge area and direction'
   novars = 2; nocoords = 1

CASE (io_dissal)
   title = 'Salinity discharge'
   novars = 1; nocoords = 1

CASE (io_distmp)
   title = 'Temperature discharge'
   novars = 1; nocoords = 1

!
!2.9 Sediment model files
!------------------------
!
   
CASE (io_sedspc,io_darspc,io_sedobc,io_sednst)
   IF (ifil.EQ.1) THEN
      CALL set_sedfiles_atts(iddesc,ifil,iotype,filepars)
   ELSE
      CALL set_sedfiles_atts(iddesc,ifil,iotype,filepars)
   ENDIF
   GOTO 1000

!
!2.10 Biological model files
!---------------------------
!

CASE (io_parabs,io_parrel,io_bioobc,io_bionst,io_parsur)
   CALL set_biofiles_atts(iddesc,ifil,iotype,filepars)
   GOTO 1000

!
!2.11 Particle model files
!-------------------------
!   

CASE (io_parspc,io_parcld,io_pargrd,io_parphs)
   CALL set_partfiles_atts(iddesc,ifil,iotype,filepars)
   GOTO 1000
   
END SELECT

!
!3. Store attributes
!-------------------
!

filepars%title = title
filepars%nocoords = nocoords
filepars%novars = novars
filepars%nosubnames = nosubnames
filepars%subnameid = subnameid
filepars%timeid = MERGE(0,1+subnameid,nocoords.EQ.0)
filepars%maxrecs = 0 

1000 CONTINUE

IF (loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF
pglev = pglev - 1


RETURN

END SUBROUTINE set_modfiles_atts

!========================================================================

SUBROUTINE set_modfiles_name(iddesc,ifil,iotype)
!************************************************************************
!
! *set_modfiles_name* Obtain the name of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.x
!
! Description -
!
! Module calls - file_suffix
!
!************************************************************************
!
USE utility_routines, ONLY: file_suffix

!* Arguments
!
INTEGER, INTENT(IN) :: iddesc,ifil, iotype

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*   INTEGER File id
!*ifil*     INTEGER File number
!*iotype*   INTEGER I/O type of file
!             = 1 => input
!             = 2 => output
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=3) :: suffix
CHARACTER (LEN=12) :: cfil
CHARACTER (LEN=leniofile) :: namedesc


pglev = pglev + 1
procname(pglev) = 'set_modfiles_name'
IF (loglev1.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ENDIF
ENDIF

!
!1. Initialise parameters
!------------------------
!

WRITE (cfil,'(I12)') ifil; cfil = ADJUSTL(cfil)

!
!2. Define attributes
!--------------------
!

SELECT CASE (iddesc)

!
!2.1 Domain decomposition
!------------------------
!

CASE (io_mppmod)
   namedesc = 'mppmod'

!
!2.2 Initial/final conditions
!----------------------------
!

CASE (io_inicon)
   SELECT CASE (ifil)
      CASE (ics_phys)
         namedesc = 'phsics'
      CASE (ics_sed)
         namedesc = 'sedics'
      CASE (ics_morph)
         namedesc = 'morics'
      CASE (ics_part)
         namedesc = 'parics'
      CASE (ics_bio)
         namedesc = 'bioics'
   END SELECT

CASE (io_fincon)
   SELECT CASE (ifil)
      CASE (ics_phys)
         namedesc = 'phsfin'
      CASE (ics_sed)
         namedesc = 'sedfin'
      CASE (ics_morph)
         namedesc = 'morfin'
      CASE (ics_bio)
         namedesc = 'biofin'
      CASE (ics_part)
         namedesc = 'parfin'
   END SELECT

!
!2.3 Model grid and data locations
!---------------------------------
!

CASE (io_modgrd)
   namedesc = 'modgrd'
CASE (io_metabs)
   namedesc = 'metabs'
CASE (io_sstabs)
   namedesc = 'sstabs'
CASE (io_wavabs)
   namedesc = 'wavabs'
CASE (io_parabs)
   namedesc = 'parabs'
CASE (io_nstgrd)
   namedesc = 'nstgrd'//TRIM(cfil)
CASE (io_metrel)
   namedesc = 'metrel'
CASE (io_sstrel)
   namedesc = 'sstrel'
CASE (io_wavrel)
   namedesc = 'wavrel'
CASE (io_parrel)
   namedesc = 'parrel'

!
!2.4 Open boundaries
!-------------------
!

CASE (io_1uvsur)
   namedesc = '1uvsur'//TRIM(cfil)
CASE (io_2uvobc)
   namedesc = '2uvobc'//TRIM(cfil)
CASE (io_2xyobc)
   namedesc = '2xyobc'//TRIM(cfil)
CASE (io_3uvobc)
   namedesc = '3uvobc'//TRIM(cfil)
CASE (io_3xyobc)
   namedesc = '3xyobc'//TRIM(cfil)
CASE (io_salobc)
   namedesc = 'salobc'//TRIM(cfil)
CASE (io_tmpobc)
   namedesc = 'tmpobc'//TRIM(cfil)
CASE (io_sedobc)
   namedesc = 'sedobc'//TRIM(cfil)
CASE (io_bioobc)
   namedesc = 'bioobc'//TRIM(cfil)

!
!2.5 Nested output
!-----------------
!

CASE (io_nstspc)
   namedesc = 'nstspc'
CASE (io_2uvnst)
   namedesc = '2uvnst'//TRIM(cfil)
CASE (io_3uvnst)
   namedesc = '3uvnst'//TRIM(cfil)
CASE (io_2xynst)
   namedesc = '2xynst'//TRIM(cfil)
CASE (io_3xynst)
   namedesc = '3xynst'//TRIM(cfil)
CASE (io_salnst)
   namedesc = 'salnst'//TRIM(cfil)
CASE (io_tmpnst)
   namedesc = 'tmpnst'//TRIM(cfil)
CASE (io_sednst)
   namedesc = 'sednst'//TRIM(cfil)
CASE (io_bionst)
   namedesc = 'bionst'//TRIM(cfil)

!
!2.6 Surface data
!----------------
!

CASE (io_metsur)
   namedesc = 'metsur'
CASE (io_sstsur)
   namedesc = 'sstsur'
CASE (io_wavsur)
   namedesc = 'wavsur'
CASE (io_parsur)
   namedesc = 'parsur'
   
!
!2.7 Structures
!--------------
!
CASE (io_drycel)
   namedesc = 'drycel'
CASE (io_thndam)
   namedesc = 'thndam'
CASE (io_weibar)
   namedesc = 'weibar'

!
!2.8 Discharges
!--------------
!

CASE (io_disspc)
   namedesc = 'disspc'
CASE (io_disloc)
   namedesc = 'disloc'//TRIM(cfil)
CASE (io_disvol)
   namedesc = 'disvol'//TRIM(cfil)
CASE (io_discur)
   namedesc = 'discur'//TRIM(cfil)
CASE (io_dissal)
   namedesc = 'dissal'//TRIM(cfil)
CASE (io_distmp)
   namedesc = 'distmp'//TRIM(cfil)
CASE (io_dissed)
   namedesc = 'distmp'//TRIM(cfil)

!
!2.9 Specifiers
!--------------
!

CASE (io_sedspc)
   namedesc = 'sedspc'
CASE (io_darspc)
   namedesc = 'darspc'
CASE (io_parspc)
   namedesc = 'parspc'

!
!2.10 Particle model
!-------------------
!

CASE (io_parcld)
   namedesc = 'parcld'//TRIM(cfil)
CASE (io_pargrd)
   namedesc = 'pargrd'
CASE (io_parphs)
   namedesc = 'parphs'
END SELECT

!
!3. Store attributes
!-------------------
!

IF (TRIM(modfiles(iddesc,ifil,iotype)%filename).EQ.'') THEN
   suffix = file_suffix(modfiles(iddesc,ifil,iotype)%form)
   SELECT CASE (modfiles(iddesc,ifil,iotype)%status)
      CASE ('N')
         modfiles(iddesc,ifil,iotype)%filename = ''
      CASE ('R','W')
         modfiles(iddesc,ifil,iotype)%filename = TRIM(intitle)//'.'//&
                        & TRIM(namedesc)//'.'//TRIM(suffix)
   END SELECT
ENDIF

IF (loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF
pglev = pglev - 1


RETURN

END SUBROUTINE set_modfiles_name

!========================================================================

SUBROUTINE set_modvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts,noprofs)
!************************************************************************
!
! *set_modvars_atts* Obtain the variable attributes of a model file
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.x
!
! Description -
!
! Module calls - inquire_var, set_biovars_atts, set_partvars_atts,
!                set_sedvars_atts
!
!************************************************************************
!
USE datatypes
USE gridpars
USE modids
USE nestgrids
USE obconds
USE paralpars
USE physpars
USE structures
USE switches
USE tide
USE biovars_routines, ONLY: set_biovars_atts
USE partvars_routines, ONLY: set_partvars_atts
USE sedvars_routines, ONLY: set_sedvars_atts

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, iotype, numvars
INTEGER, INTENT(IN), OPTIONAL :: noprofs
TYPE (FileParams), INTENT(INOUT) :: filepars
TYPE (VariableAtts), INTENT(INOUT), DIMENSION(numvars) :: varatts

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER File id
!*ifil*      INTEGER File number
!*iotype*    INTEGER I/O type of file
!              = 1 => input
!              = 2 => output
!*filepars*  DERIVED Attributes of the forcing file
!*numvars*   INTEGER Number of variables in the data file
!*varatts*   DERIVED Attributes of the variables in the forcing file
!*noprofs*   INTEGER Number of open boundary profiles in the forcing file
!
!------------------------------------------------------------------------------
!
!*Local variables
!
CHARACTER (LEN=lenname), DIMENSION(MaxIODims) :: dimnames
CHARACTER (LEN=lenname), DIMENSION(MaxSubVars) :: subname
INTEGER :: idim, idtime, igrd, isub, ivar, n, nhtype, nocoords, nodims, &
         & nofiles, nosubnames, novars, nvars, n1dat, n2dat
INTEGER, DIMENSION(numvars) :: ivarid, kind_type, nrank
INTEGER, DIMENSION(MaxIODims) :: dimvals
INTEGER, DIMENSION(numvars,MaxIODims) :: iddims


pglev = pglev + 1
procname(pglev) = 'set_modvars_atts'

IF (loglev1.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, TRIM(procname(pglev))
   ENDIF
ENDIF

!
!1. Initialise
!-------------
!

kind_type = MERGE(real_type,rlong_type,filepars%rtype.EQ.real_type)
idtime = 0; ivarid = 0; nrank = -1; iddims = 0; nodims = 0
nocoords = filepars%nocoords
dimvals = 0; dimnames = ''
nofiles = maxdatafiles(iddesc,1)
novars = filepars%novars
nosubnames = filepars%nosubnames

!
!2. Define variable attributes
!-----------------------------
!

SELECT CASE (iddesc)
   
!
!2.1 Domain decomposition
!------------------------
!

CASE (io_mppmod)
   nodims = 2
   dimvals(1:2) = (/nprocs,nomglevels/)
   dimnames(1:2) = (/'nprocs    ','nomglevels'/)
   ivarid(1:4) = (/iarr_mgvars_nc1procs,iarr_mgvars_nc2procs,&
                 & iarr_mgvars_nr1procs,iarr_mgvars_nr2procs/)
   kind_type(1:4) = int_type
   nrank(1:4) = 2
   iddims(1:4,1) = 1
   iddims(1:4,2) = 2

!
!2.2 Initial conditions
!----------------------
!

CASE (io_inicon,io_fincon)

   SELECT CASE (ifil)

      CASE (ics_phys)
         nodims = 14
         dimvals(1:14) = (/nc,nr,nc-1,nr-1,nz,nz+1,nconobc,numwbaru,numwbarv,&
                         & nobu,nobv,1,2,3/)
         dimnames(1:14) = (/'nc      ','nr      ','nc-1    ','nr-1    ',&
                          & 'nz      ','nz+1    ','nconobc ','numwbaru',&
                          & 'numwbarv','nobu    ','nobv    ','one     ',&
                          & 'two     ','three   '/)
            
!        ---2-D mode
         ivar = 2
         IF (iopt_mode_2D.GT.0) THEN
            ivarid(ivar:ivar+2) = (/iarr_udvel,iarr_vdvel,iarr_zeta/)
            nrank(ivar:ivar+2) = 2
            iddims(ivar,1:2) = (/1,4/)
            iddims(ivar+1,1:2) = (/3,2/) 
            iddims(ivar+2,1:2) = (/3,4/)
            ivar = ivar + 3
         ENDIF

!        ---3-D currents
         IF (iopt_mode_3D.GT.0) THEN
            ivarid(ivar:ivar+1) = (/iarr_uvel,iarr_vvel/)
            nrank(ivar:ivar+1) = 3
            iddims(ivar,1:3) = (/1,4,5/)
            iddims(ivar+1,1:3) = (/3,2,5/)
            ivar = ivar + 2
            IF (iopt_grid_nodim.EQ.3) THEN
               ivarid(ivar) = iarr_wvel
               nrank(ivar) = 3
               iddims(ivar,1:3) = (/3,4,6/)
               ivar = ivar + 1
            ENDIF
         ENDIF

!        ---density arrays
         IF (iopt_temp.GT.0) THEN
            ivarid(ivar) = iarr_temp
            nrank(ivar) = 3
            iddims(ivar,1:3) = (/3,4,5/)
            ivar = ivar + 1
         ENDIF
         IF (iopt_sal.GT.0) THEN         
            ivarid(ivar) = iarr_sal
            nrank(ivar) = 3
            iddims(ivar,1:3) = (/3,4,5/)
            ivar = ivar + 1
         ENDIF
         
!        ---turbulence arrays
         IF (iopt_vdif_coef.EQ.3) THEN
            ivarid(ivar) = iarr_tke
            nrank(ivar) = 3
            iddims(ivar,1:3) = (/3,4,6/)
            ivar = ivar + 1
            IF (iopt_turb_ntrans.GT.0) THEN
               IF (iopt_turb_param.GT.0) THEN
                  IF (iopt_turb_param.EQ.1) THEN
                     ivarid(ivar) = iarr_zlmix
                  ELSEIF (iopt_turb_param.EQ.2) THEN
                     ivarid(ivar) = iarr_dissip
                  ENDIF
                  nrank(ivar) = 3
                  iddims(ivar,1:3) = (/3,4,6/)
                  ivar = ivar + 1
               ENDIF
            ENDIF
         ENDIF

!        ---bottom stress arrays
         IF (iopt_bstres_form.EQ.2) THEN
            IF (iopt_bstres_drag.EQ.2) THEN
               ivarid(ivar) = iarr_bdragcoefatc
               nrank(ivar) = 2
               iddims(ivar,1:2) = (/3,4/)
               ivar = ivar + 1
            ELSEIF (iopt_bstres_drag.EQ.4) THEN
               ivarid(ivar) = iarr_zroughatc
               nrank(ivar) = 2
               iddims(ivar,1:2) = (/3,4/)
               ivar = ivar + 1
            ENDIF

         ENDIF
         
!        ---tidal arrays
         IF (iopt_astro_pars.EQ.1.AND.nconobc.GT.0) THEN
            ivarid(ivar) = iarr_phase_obc
            nrank(ivar) = 1
            iddims(ivar,1) = 7
            ivar = ivar + 1
         ENDIF

!        ---energy losses at weirs
         IF (iopt_weibar.EQ.1) THEN
            IF (numwbaru.GT.0) THEN
               ivarid(ivar) = iarr_wbarelossu
               nrank(ivar) = 1
               iddims(ivar,1) = 8
               ivar = ivar + 1
            ENDIF
            IF (numwbarv.GT.0) THEN
               ivarid(ivar) = iarr_wbarelossv
               nrank(ivar) = 1
               iddims(ivar,1) = 9
               ivar = ivar + 1
            ENDIF
         ENDIF
         
!        ---open boundary arrays
         IF (obcsalatu_init) THEN
            ivarid(ivar) = iarr_obcsalatu
            nrank(ivar) = 4
            iddims(ivar,1:4) = (/10,5,14,12/)
            ivar = ivar + 1
         ENDIF
         IF (obcsalatv_init) THEN
            ivarid(ivar) = iarr_obcsalatv
            nrank(ivar) = 4
            iddims(ivar,1:4) = (/11,5,14,12/)
            ivar = ivar + 1
         ENDIF
         IF (obctmpatu_init) THEN
            ivarid(ivar) = iarr_obctmpatu
            nrank(ivar) = 4
            iddims(ivar,1:4) = (/10,5,14,12/)
            ivar = ivar + 1
         ENDIF
         IF (obctmpatv_init) THEN
            ivarid(ivar) = iarr_obctmpatv
            nrank(ivar) = 4
            iddims(ivar,1:4) = (/11,5,14,12/)
            ivar = ivar + 1
         ENDIF
         IF (obc3uvatu_init) THEN
            ivarid(ivar) = iarr_obc3uvatu
            nrank(ivar) = 3
            iddims(ivar,1:3) = (/10,5,13/)
            ivar = ivar + 1
         ENDIF
         IF (obc3uvatv_init) THEN
            ivarid(ivar) = iarr_obc3uvatv
            nrank(ivar) = 3
            iddims(ivar,1:3) = (/11,5,13/)
            ivar = ivar + 1
         ENDIF
         IF (obc2uvatu_init) THEN
            ivarid(ivar) = iarr_obc2uvatu
            nrank(ivar) = 2
            iddims(ivar,1:2) = (/10,13/)
            ivar = ivar + 1
         ENDIF
         IF (obc2uvatv_init) THEN
            ivarid(ivar) = iarr_obc2uvatv
            nrank(ivar) = 2
            iddims(ivar,1:2) = (/11,13/)
            ivar = ivar + 1
         ENDIF

      CASE (ics_sed,ics_morph)
         CALL set_sedvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts)
         GOTO 1000
      CASE (ics_bio)
         CALL set_biovars_atts(iddesc,ifil,iotype,filepars,numvars,varatts)
         GOTO 1000
      CASE (ics_part)
         CALL set_partvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts)
         GOTO 1000
   END SELECT

!
!2.3 Model grid and data locations
!---------------------------------
!
!2.3.1 Model
!-----------
!

CASE (io_modgrd)
   nhtype = iopt_grid_htype
   ivar = 1
   nodims = 7
   dimvals(1:7) = (/nc,nr,nc-1,nr-1,nz+1,nobu,nobv/)
   dimnames(1:7) = (/'nc  ','nr  ','nc-1','nr-1','nz+1','nobu','nobv'/)
   IF (nhtype.EQ.2) THEN
      ivarid(1:2) = (/iarr_gdelxglb,iarr_gdelyglb/)
      nrank(1:2) = 1
      iddims(1,1) = 1
      iddims(2,1) = 2
      ivar = ivar + 2
   ELSEIF (nhtype.EQ.3) THEN
      ivarid(1:2) = (/iarr_gxcoordglb,iarr_gycoordglb/)
      nrank(1:2) = 2
      iddims(1,1:2) = (/1,2/)
      iddims(2,1:2) = (/1,2/)
      ivar = ivar + 2
   ENDIF
   IF (iopt_grid_vtype.EQ.2) THEN
      ivarid(ivar) = iarr_gsigcoordatw
      nrank(ivar) = 1
      iddims(ivar,1) = 5
      ivar = ivar + 1
   ELSEIF (iopt_grid_vtype.EQ.3) THEN
      ivarid(ivar) = iarr_gscoordglb
      nrank(ivar) = 3
      iddims(ivar,1:3) = (/3,4,5/)
      ivar = ivar + 1
   ENDIF
   ivarid(ivar) = iarr_depmeanglb
   nrank(ivar) = 2
   iddims(ivar,1:2) = (/3,4/)
   ivar = ivar + 1
   IF (nobu.GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_iobu,iarr_jobu/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 6
      kind_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
   ENDIF
   IF (nobv.GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_iobv,iarr_jobv/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 7
      kind_type(ivar:ivar+1) = int_type
   ENDIF

!
!2.3.2 Surface grids
!-------------------
!
!---using absolute coordinates
CASE (io_metabs,io_sstabs,io_wavabs)
   IF (iddesc.EQ.io_metabs) igrd = igrd_meteo 
   IF (iddesc.EQ.io_sstabs) igrd = igrd_sst  
   IF (iddesc.EQ.io_wavabs) igrd = igrd_waves  
   n1dat = surfacegrids(igrd,1)%n1dat
   n2dat = surfacegrids(igrd,1)%n2dat
   nodims = 2
   dimvals(1:2) = (/n1dat,n2dat/)
   dimnames(1:2) = (/'n1dat','n2dat'/) 
   ivarid(1:2) = (/iarr_xcoord,iarr_ycoord/)
   nrank(1:2) = 2
   iddims(1:2,1) = 1
   iddims(1:2,2) = 2

!---using relative coordinates   
CASE (io_metrel,io_sstrel,io_wavrel)
   nodims = 3
   dimvals(1:3) = (/nc,nr,2/)
   dimnames(1:3) = (/'nc ','nr ','two'/)
   ivarid(1:3) = (/iarr_icoordC,iarr_jcoordC,iarr_weightsC/)
   nrank(1:2) = 2
   iddims(1:2,1) = 1
   iddims(1:2,2) = 2
   kind_type(1:2) = int_type
   nrank(3) = 4
   iddims(3,1:4) = (/3,3,1,2/)

!
!2.3.3 Nested grids
!------------------
!

CASE (io_nstgrd)

   nodims = 6
   dimvals(1:6) = (/nohnstglbc(ifil),nohnstglbu(ifil),nohnstglbv(ifil),&
                  & nohnstglbx(ifil),nohnstglby(ifil),novnst(ifil)/)
   dimnames(1:6) = (/'nohnstglbc','nohnstglbu','nohnstglbv','nohnstglbx',&
                   & 'nohnstglby','novnst    '/)
   ivar = 1
   
!  ---C-nodes
   IF (nohnstglbc(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordatc,iarr_ycoordatc/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 1
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_scoordatc
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/1,6/)
         ivar = ivar + 1
      ENDIF
   ENDIF
   
!  ---U-nodes
   IF (nohnstglbu(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordatu,iarr_ycoordatu/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 2
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_scoordatu
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/2,6/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---V-nodes
   IF (nohnstglbv(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordatv,iarr_ycoordatv/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 3
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_scoordatv
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/3,6/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---X-nodes
   IF (nohnstglbx(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordatx,iarr_ycoordatx/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 4
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_scoordatx
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/4,6/)
         ivar = ivar + 1
      ENDIF
   ENDIF

!  ---Y-nodes
   IF (nohnstglby(ifil).GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_xcoordaty,iarr_ycoordaty/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 5
      ivar = ivar + 2
      IF (novnst(ifil).GT.0) THEN
         ivarid(ivar) = iarr_scoordaty
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/5,6/)
         ivar = ivar + 1
      ENDIF
   ENDIF
   
!
!2.4 Open (surface) boundary conditions
!--------------------------------------
!
!2.4.1 1-D case
!--------------
!

CASE (io_1uvsur)
   
   IF (ifil.EQ.1) THEN
      nodims = 1
      dimvals(1) = nconobc
      dimnames(1) = 'nconobc'
      ivar = 1
      IF (iopt_sbc_1D.EQ.1.OR.iopt_sbc_1D.EQ.3) THEN
         ivarid(ivar:ivar+3) = (/iarr_gxslope_amp,iarr_gxslope_pha,&
                               & iarr_gyslope_amp,iarr_gyslope_pha/)
         nrank(ivar:ivar+3) = 1
         iddims(ivar:ivar+3,1) = 1
         ivar = ivar + 4
      ENDIF
      IF (iopt_sbc_1D.EQ.1.OR.iopt_sbc_1D.EQ.2) THEN
         ivarid(ivar:ivar+1) = (/iarr_zeta_amp,iarr_zeta_pha/)
         nrank(ivar:ivar+1) = 1
         iddims(ivar:ivar+1,1) = 1
      ENDIF
      
   ELSEIF (ifil.EQ.2) THEN
      SELECT CASE (iopt_sbc_1D)
         CASE (1); nvars = 3
         CASE (2); nvars = 1
         CASE (3); nvars = 2
      END SELECT
      nodims = 1
      dimvals(1) = nvars
      dimnames(1) = 'novars'
      ivarid(2) = iarr_surdat1d
      nrank(2) = 1
      iddims(2,1) = 1
      
   ENDIF

!
!2.4.2 2-D case (normal)
!-----------------------
!

CASE (io_2uvobc)
   
   IF (ifil.EQ.1) THEN
      ivar = 1
      nodims = 8
      dimvals(1:8) = (/nobu,nobv,nofiles-1,nconobc,nobu+nobv,nqsecobu,&
                     & nqsecobv,2/)
      dimnames(1:2) = (/'nobu','nobv'/)
      dimnames(3) = 'nofiles-1'
      dimnames(4:7) = (/'nconobc ','nodat   ','nqsecobu','nqsecobv'/)
      dimnames(8) = 'two'

!     ---type of open boundary conditions      
      IF (nobu.GT.0) THEN
         ivarid(ivar) = iarr_ityp2dobu 
         nrank(ivar) = 1
         iddims(ivar,1) = 1
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (nobv.GT.0) THEN
         ivarid(ivar) = iarr_ityp2dobv 
         nrank(ivar) = 1
         iddims(ivar,1) = 2
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF

!     ---position of elevation points      
      IF (nobu.GT.0) THEN
         ivarid(ivar) = iarr_iloczobu 
         nrank(ivar) = 1
         iddims(ivar,1) = 1
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (nobv.GT.0) THEN
         ivarid(ivar) = iarr_iloczobv 
         nrank(ivar) = 1
         iddims(ivar,1) = 2
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      
!     ---number of data and index maps per input file
      IF (nofiles.GT.1) THEN
         ivarid(ivar:ivar+2) = (/iarr_no2dobuv,iarr_iobc2dtype,&
                               & iarr_index2dobuv/)
         nrank(ivar:ivar+2) = (/1,1,2/)
         iddims(ivar:ivar+1,1) = 3
         iddims(ivar+2,1:2) = (/5,3/)
         kind_type(ivar:ivar+2) = int_type
         ivar = ivar + 3
      ENDIF

!     ---amplitudes and phases
      IF (nobu.GT.0.AND.nconobc.GT.0) THEN
         ivarid(ivar:ivar+3) = (/iarr_udatobu_amp,iarr_zdatobu_amp,&
                               & iarr_udatobu_pha,iarr_zdatobu_pha/)
         nrank(ivar:ivar+3) = 2
         iddims(ivar:ivar+3,1) = 1
         iddims(ivar:ivar+3,2) = 4
         ivar = ivar + 4
      ENDIF
      IF (nobv.GT.0.AND.nconobc.GT.0) THEN
         ivarid(ivar:ivar+3) = (/iarr_vdatobv_amp,iarr_zdatobv_amp,&
                               & iarr_vdatobv_pha,iarr_zdatobv_pha/)
         nrank(ivar:ivar+3) = 2
         iddims(ivar:ivar+3,1) = 2
         iddims(ivar:ivar+3,2) = 4
         ivar = ivar + 4
      ENDIF

!     ---discharges along open boundary sections
      IF (nqsecobu.GT.0) THEN
         ivarid(ivar) = iarr_iqsecobu
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/6,8/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (nqsecobv.GT.0) THEN
         ivarid(ivar) = iarr_jqsecobv
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/7,8/)
         kind_type(ivar) = int_type
      ENDIF

   ELSE
!     ---data file
      nodims = 2
      nvars = MERGE(2,1,iobc2dtype(ifil).EQ.1)
      dimvals(1:2) = (/no2dobuv(ifil),nvars/)
      dimnames(1) = 'nodat'
      dimnames(2) = 'novars'
      ivarid(2) = iarr_obcdatuv2d
      nrank(2) = 2
      iddims(2,1:2) = (/1,2/)
   ENDIF

!
!2.4.3 2-D case (tangential)
!---------------------------
!

CASE (io_2xyobc)
   IF (ifil.EQ.1) THEN
      nodims = 5
      dimvals(1:5) = (/nobx,noby,nofiles-1,nconobc,nobx+noby/)
      dimnames(1:2) = (/'nobx','noby'/)
      dimnames(3) = 'nofiles-1'
      dimnames(4:5) = (/'nconobc','nodat  '/)
      ivar = 1
      IF (nobx.GT.0) THEN
         ivarid(ivar) = iarr_ityp2dobx
         nrank(ivar) = 1
         iddims(ivar,1) = 1 
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (noby.GT.0) THEN
         ivarid(ivar) = iarr_ityp2doby
         nrank(ivar) = 1
         iddims(ivar,1) = 2
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      
!     ---number of data and index maps per input file
      IF (nofiles.GT.1) THEN
         ivarid(ivar:ivar+1) = (/iarr_no2dobxy,iarr_index2dobxy/)
         nrank(ivar:ivar+1) = (/1,2/)
         iddims(ivar,1) = 3
         iddims(ivar+1,1:2) = (/5,3/)
         kind_type(ivar:ivar+1) = int_type
         ivar = ivar + 2
      ENDIF

!     ---amplitudes and phases
      IF (nobx.GT.0.AND.nconobc.GT.0) THEN
         ivarid(ivar:ivar+1) = (/iarr_vdatobx_amp,iarr_vdatobx_pha/)
         nrank(ivar:ivar+1) = 2
         iddims(ivar:ivar+1,1) = 1
         iddims(ivar:ivar+1,2) = 4
         ivar = ivar + 2
      ENDIF
      IF (noby.GT.0.AND.nconobc.GT.0) THEN
         ivarid(ivar:ivar+1) = (/iarr_udatoby_amp,iarr_udatoby_pha/)
         nrank(ivar:ivar+1) = 2
         iddims(ivar:ivar+1,1) = 2
         iddims(ivar:ivar+1,2) = 4
      ENDIF

   ELSE
!     ---data file
      nodims = 2
      dimvals(1:2) = (/no2dobxy(ifil),1/)
      dimnames(1) = 'nodat'
      dimnames(2) = 'novars'
      ivarid(2) = iarr_obcdatxy2d
      nrank(2) = 2
      iddims(2,1:2) = (/1,2/)

   ENDIF

!
!2.4.4 3-D case (normal)
!-----------------------
!

CASE (io_3uvobc,io_salobc,io_tmpobc)   
   IF (ifil.EQ.1) THEN
      nodims = 5
      dimvals(1:5) = (/nobu,nobv,nofiles-1,nobu+nobv,1/)
      dimnames(1:2) = (/'nobu','nobv'/)
      dimnames(3) = 'nofiles-1'
      dimnames(4) = 'nodat'
      dimnames(5) = 'novars' 
      ivar = 1

!     ---type of boundary condition      
      IF (nobu.GT.0) THEN
         ivarid(ivar) = iarr_itypobu
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/1,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (nobv.GT.0) THEN
         ivarid(ivar) = iarr_itypobv
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/2,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF

!     ---profile numbers
      IF (nobu.GT.0) THEN
         ivarid(ivar) = iarr_iprofobu
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/1,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (nobv.GT.0) THEN
         ivarid(ivar) = iarr_iprofobv
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/2,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF

!     ---number of profiles per file and index mapping
      IF (nofiles.GT.1) THEN
         ivarid(ivar) = iarr_noprofsd
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/3,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
         ivarid(ivar) = iarr_indexprof
         nrank(ivar)  = 3
         iddims(ivar,1:3) = (/4,3,5/) 
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      
   ELSE
!     ---data file
      nodims = 3
      dimvals(1:3) = (/noprofs,nz,1/)
      dimnames(1:3) = (/'noprofs','nz     ','novars '/)
      IF (iddesc.EQ.io_3uvobc) THEN
         ivarid(2) = iarr_prof3dveluv
      ELSEIF (iddesc.EQ.io_salobc) THEN
         ivarid(2) = iarr_prof3dsal
      ELSEIF (iddesc.EQ.io_tmpobc) THEN
         ivarid(2) = iarr_prof3dtmp
      ENDIF
      nrank(2) = 3
      iddims(2,1:3) = (/1,2,3/)
      
   ENDIF

!
!2.4.5 3-D case (tangential)
!---------------------------
!

CASE (io_3xyobc)
   IF (ifil.EQ.1) THEN
      nodims = 5
      dimvals(1:5) = (/nobx,noby,nofiles-1,nobx+noby,1/)
      dimnames(1:2) = (/'nobx  ','noby  '/)
      dimnames(3) = 'nofiles-1'
      dimnames(4) = 'nodat'
      dimnames(5) = 'novars' 
      ivar = 1

!     ---type of boundary condition      
      IF (nobx.GT.0) THEN
         ivarid(ivar) = iarr_itypobx
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/1,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (noby.GT.0) THEN
         ivarid(ivar) = iarr_itypoby
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/2,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF

!     ---profile numbers      
      IF (nobx.GT.0) THEN
         ivarid(ivar) = iarr_iprofobx
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/1,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF
      IF (noby.GT.0) THEN
         ivarid(ivar) = iarr_iprofoby
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/2,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
      ENDIF

!     ---number of profiles per file and index mapping      
      IF (nofiles.GT.1) THEN
         ivarid(ivar) = iarr_noprofsd
         nrank(ivar) = 2
         iddims(ivar,1:2) = (/3,5/)
         kind_type(ivar) = int_type
         ivar = ivar + 1
         ivarid(ivar) = iarr_indexprof 
         nrank(ivar)  = 3
         iddims(ivar,1:3) = (/4,3,5/)
         kind_type(ivar) = int_type
      ENDIF
      
   ELSE
!     ---data file
      nodims = 3
      dimvals(1:3) = (/noprofs,nz,1/)
      dimnames(1:3) = (/'noprofs','nz     ','novars '/)
      ivarid(2) = iarr_prof3dvelxy
      nrank(2) = 3
      iddims(2,1:3) = (/1,2,3/)
   ENDIF

!
!2.5 Nested output
!-----------------
!

CASE (io_nstspc)
   nodims = 1
   dimvals(1) = nonestsets
   dimnames(1) = 'nonestsets'
   ivarid(1:7) = (/iarr_nohnstglbc,iarr_nohnstglbu,iarr_nohnstglbv,&
                 & iarr_nohnstglbx,iarr_nohnstglby,iarr_novnst,iarr_inst2dtype/)
   nrank(1:7) = 1
   iddims(1:7,1) = 1
   kind_type(1:7) = int_type

CASE (io_2uvnst)
   nodims = 2
   nvars = MERGE(2,1,inst2dtype(ifil).EQ.1)
   dimvals(1:2) = (/noprofs,nvars/)
   dimnames(1) = 'noprofs'
   dimnames(2) = 'novars'
   ivarid(2) = iarr_obcdatuv2d
   nrank(2) = 2
   iddims(2,1:2) = (/1,2/)
   
CASE (io_2xynst)
   nodims = 2
   dimvals(1:2) = (/noprofs,1/)
   dimnames(1) = 'noprofs'
   dimnames(2) = 'novars'
   ivarid(2) = iarr_obcdatxy2d
   nrank(2) = 2
   iddims(2,1:2) = (/1,2/)

CASE (io_3uvnst)
   nodims = 3
   dimvals(1:3) = (/noprofs,novnst(ifil),1/)
   dimnames(1:3) = (/'noprofs','nz     ','novars '/)
   ivarid(2) = iarr_prof3dveluv
   nrank(2) = 3
   iddims(2,1:3) = (/1,2,3/)

CASE (io_3xynst)
   nodims = 3
   dimvals(1:3) = (/noprofs,novnst(ifil),1/)
   dimnames(1:3) = (/'noprofs','nz     ','novars '/)
   ivarid(2) = iarr_prof3dvelxy
   nrank(2) = 3
   iddims(2,1:3) = (/1,2,3/)

CASE (io_salnst)
   nodims = 3
   dimvals(1:3) = (/noprofs,novnst(ifil),1/)
   dimnames(1:3) = (/'noprofs','nz     ','novars '/)
   ivarid(2) = iarr_prof3dsal
   nrank(2) = 3
   iddims(2,1:3) = (/1,2,3/)

CASE (io_tmpnst)
   nodims = 3
   dimvals(1:3) = (/noprofs,novnst(ifil),1/)
   dimnames(1:3) = (/'noprofs','nz     ','novars '/)
   ivarid(2) = iarr_prof3dtmp
   nrank(2) = 3
   iddims(2,1:3) = (/1,2,3/)

!
!2.6 Surface data
!----------------
!
!---meteo
CASE (io_metsur)

   nhtype = surfacegrids(igrd_meteo,1)%nhtype
   n1dat = MERGE(nc,surfacegrids(igrd_meteo,1)%n1dat,nhtype.EQ.4)
   n2dat = MERGE(nr,surfacegrids(igrd_meteo,1)%n2dat,nhtype.EQ.4)

   IF (nhtype.EQ.0) THEN

      nodims = 1
      dimnames(1) = 'nosubnames' 
      isub = 1
      IF (iopt_meteo_stres.EQ.1) THEN
         IF (iopt_meteo_data.EQ.1) THEN
            subname(isub:isub+1) = (/'uwindatc','vwindatc'/)
            isub = isub + 2
         ELSEIF (iopt_meteo_data.EQ.2) THEN
            subname (isub:isub+1) = (/'usstresatc','vsstresatc'/)
            isub = isub + 2
         ENDIF
      ENDIF
      IF (iopt_meteo_pres.EQ.1) THEN
         subname (isub) = 'atmpres'
         isub = isub + 1
      ENDIF
      IF (iopt_meteo_heat.EQ.1) THEN
         SELECT CASE (iopt_meteo_data)
            CASE (1)
               subname(isub:isub+2) = (/'airtemp    ','relhum     ',&
                                      & 'cloud_cover'/)
               isub = isub + 3
            CASE (2)
               subname(isub:isub+1) = (/'qnonsol','qrad   '/)
               isub = isub + 2
         END SELECT
      ENDIF
      SELECT CASE (iopt_meteo_precip)
         CASE (1)
            subname(isub) = 'evapminprec'
            isub = isub + 1
         CASE (2)
            subname(isub) = 'precipitation'
            isub = isub + 1
      END SELECT
      nosubnames = isub - 1
      nodims = 1
      dimvals(1) = nosubnames
      dimnames(1) = 'nosubnames'
      ivarid(3) = iarr_meteodata
      nrank(3) = 1
      iddims(3,1) = 1
      
   ELSE
      
      nodims = 2
      dimvals(1:2) = (/n1dat,n2dat/)
      dimnames(1:2) = (/'n1dat','n2dat'/) 
      nrank(2:numvars) = 2
      iddims(2:numvars,1) = 1
      iddims(2:numvars,2) = 2
      ivar = 2
      IF (iopt_meteo_stres.EQ.1) THEN
         IF (iopt_meteo_data.EQ.1) THEN
            ivarid(ivar:ivar+1) = (/iarr_uwindatc,iarr_vwindatc/)
            ivar = ivar + 2
         ELSEIF (iopt_meteo_data.EQ.2) THEN
            ivarid(ivar:ivar+1) = (/iarr_usstresatc,iarr_vsstresatc/)
            ivar = ivar + 2
         ENDIF
      ENDIF
      IF (iopt_meteo_pres.EQ.1) THEN
         ivarid(ivar) = iarr_atmpres
         ivar = ivar + 1
      ENDIF
      IF (iopt_meteo_heat.EQ.1) THEN
         SELECT CASE (iopt_meteo_data)
            CASE (1)
               ivarid(ivar:ivar+2) = (/iarr_airtemp,iarr_relhum,&
                                     & iarr_cloud_cover/)
               ivar = ivar + 3
            CASE (2)
               ivarid(ivar:ivar+1) = (/iarr_qnonsol,iarr_qrad/)
               ivar = ivar + 2
         END SELECT
      ENDIF
      SELECT CASE (iopt_meteo_precip)
         CASE (1)
            ivarid(ivar) = iarr_evapminprec
         CASE (2)
            ivarid(ivar) = iarr_precipitation
      END SELECT

   ENDIF
      
!---sst
CASE (io_sstsur)
   
   nhtype = surfacegrids(igrd_sst,1)%nhtype
   n1dat = MERGE(nc,surfacegrids(igrd_sst,1)%n1dat,nhtype.EQ.4)
   n2dat = MERGE(nr,surfacegrids(igrd_sst,1)%n2dat,nhtype.EQ.4)
   ivarid(2) = iarr_sst
   IF (nhtype.EQ.0) THEN
      nodims = 0
      nrank(2) = 0
   ELSE
      dimvals(1:2) = (/n1dat,n2dat/)
      dimnames(1:2) = (/'n1dat','n2dat'/) 
      nrank(2) = 2
      iddims(2,1:2) = (/1,2/)
   ENDIF
   
!---waves
CASE (io_wavsur)

   nhtype = MERGE(surfacegrids(igrd_waves,1)%nhtype,4,iopt_waves_couple.EQ.0)
   n1dat = MERGE(nc,surfacegrids(igrd_waves,1)%n1dat,nhtype.EQ.4)
   n2dat = MERGE(nr,surfacegrids(igrd_waves,1)%n2dat,nhtype.EQ.4)


   IF (nhtype.EQ.0) THEN

      nodims = 1
      dimnames(1) = 'nosubnames'

      isub = 1
      subname(isub:isub+2) = (/'waveheight','waveperiod','wavedir   '/)
      isub = isub + 3
      IF (iopt_waves_form.EQ.2) THEN
         subname(isub:isub+1) = (/'wavevel   ','waveexcurs'/)
         isub = isub + 2
         IF (iopt_waves_curr.EQ.1) THEN
            subname(isub:isub+1) = (/'umstokesatc','vmstokesatc'/)
            isub = isub + 2
            IF (iopt_waves_pres.EQ.1) THEN
               subname(isub) = 'wavepres'
               isub = isub + 1
            ENDIF
         ENDIF
      ENDIF
      IF (iopt_waves_dissip.EQ.1) THEN
         subname(isub:isub+3) = (/'umswdissipatc','vmswdissipatc',&
                              &   'umbwdissipatc','vmbwdissipatc'/)
         isub = isub + 4
      ENDIF
      nosubnames = isub - 1
      nodims = 1
      dimvals(1) = nosubnames
      dimnames(1) = 'nosubnames'
      ivarid(3) = iarr_wavedata
      nrank(3) = 1
      iddims(3,1) = 1
      
   ELSE

      nodims = 2
      dimvals(1:2) = (/n1dat,n2dat/)
      dimnames(1:2) = (/'n1dat','n2dat'/) 
      nrank(2:numvars) = 2
      iddims(2:numvars,1) = 1
      iddims(2:numvars,2) = 2
      ivar = 2
      ivarid(ivar:ivar+2) = (/iarr_waveheight,iarr_waveperiod,iarr_wavedir/)
      ivar = ivar + 3
      IF (iopt_waves_form.EQ.2) THEN
         ivarid(ivar:ivar+1) = (/iarr_wavevel,iarr_waveexcurs/)
         ivar = ivar + 2
         IF (iopt_waves_curr.EQ.1) THEN
            ivarid(ivar:ivar+1) = (/iarr_umstokesatc,iarr_vmstokesatc/)
            ivar = ivar + 2
            IF (iopt_waves_pres.EQ.1) THEN
               ivarid(ivar) = iarr_wavepres
               ivar = ivar + 1
            ENDIF
         ENDIF
      ENDIF
      IF (iopt_waves_dissip.EQ.1) THEN
         ivarid(ivar:ivar+3) = (/iarr_umswdissipatc,iarr_vmswdissipatc,&
                               & iarr_umbwdissipatc,iarr_vmbwdissipatc/)
      ENDIF

   ENDIF
   
!
!2.7 Structures
!--------------
!

CASE (io_drycel)
   nodims = 1
   dimvals(1) = numdry
   dimnames(1) = 'numdry'
   ivarid(1:2) = (/iarr_idry,iarr_jdry/)
   nrank(1:2) = 1
   iddims(1:2,1) = 1
   kind_type(1:2) = int_type

CASE (io_thndam)
   nodims = 2
   dimvals(1:2) = (/numthinu,numthinv/)
   dimnames(1:2) = (/'numthinu','numthinv'/)
   ivar = 1
   IF (numthinu.GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_ithinu,iarr_jthinu/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 1
      kind_type(ivar:ivar+1) = int_type
      ivar = ivar + 2
   ENDIF
   IF (numthinv.GT.0) THEN
      ivarid(ivar:ivar+1) = (/iarr_ithinv,iarr_jthinv/)
      nrank(ivar:ivar+1) = 1
      iddims(ivar:ivar+1,1) = 2
      kind_type(ivar:ivar+1) = int_type
   ENDIF

CASE (io_weibar)
   nodims = 2
   dimvals(1:2) = (/numwbaru,numwbarv/)
   dimnames(1:2) = (/'numwbaru','numwbarv'/)
   ivar = 1
   IF (numwbaru.GT.0) THEN
      ivarid(ivar:ivar+7) = (/iarr_iwbaru,iarr_jwbaru,iarr_oricoefu,&
                            & iarr_oriheightu,iarr_oriwidthu,iarr_wbarcoefu,&
                            & iarr_wbarcrestu,iarr_wbarmodlu/)
      nrank(ivar:ivar+7) = 1
      iddims(ivar:ivar+7,1) = 1
      kind_type(ivar:ivar+1) = int_type
      ivar = ivar + 8
   ENDIF
   IF (numwbarv.GT.0) THEN
      ivarid(ivar:ivar+7)=(/iarr_iwbarv,iarr_jwbarv,iarr_oricoefv,&
                          & iarr_oriheightv,iarr_oriwidthv,iarr_wbarcoefv,&
                          & iarr_wbarcrestv,iarr_wbarmodlv/)
      nrank(ivar:ivar+7) = 1
      iddims(ivar:ivar+7,1) = 2
      kind_type(ivar:ivar+1) = int_type
   ENDIF

CASE (io_mpvcov)
   nodims        = 2
   dimvals(1:2)  = (/nc-1,nr-1/)
   dimnames(1:2) = (/'nc-1','nr-1'/)
   ivar          = 1
   ivarid(1)     = iarr_mpvcov
   nrank(1)      = 2
   iddims(1,1:2) = (/1,2/)
!
!2.8 Discharges
!--------------
!

CASE (io_disspc)
   nodims = 1
   dimvals(1) = numdis
   dimnames(1) = 'numdis'
   ivarid(1) = iarr_kdistype
   nrank(1) = 1
   iddims(1,1) = 1
   kind_type(1) = int_type

CASE (io_disloc)
   nodims = 1
   dimvals(1) = numdis
   dimnames(1) = 'numdis'
   ivarid(2:4) = (/iarr_xdiscoord,iarr_ydiscoord,iarr_zdiscoord/)
   nrank(2:4) = 1
   iddims(2:4,1) = 1

CASE (io_disvol)
   nodims = 1
   dimvals(1) = numdis
   dimnames(1) ='numdis'
   ivarid(2) = iarr_disvol
   nrank(2) = 1
   iddims(2,1) = 1

CASE (io_discur)
   nodims = 1
   dimvals(1) = numdis
   dimnames(1) = 'numdis'
   ivarid(2:3) = (/iarr_disarea,iarr_disdir/)
   nrank(2:3) = 1
   iddims(2:3,1) = 1

CASE (io_dissal)
   nodims = 1
   dimvals(1) = numdis
   dimnames(1) = 'numdis'
   ivarid(2) = iarr_dissal
   nrank(2) = 1
   iddims(2,1) = 1

CASE (io_distmp)
   nodims = 1
   dimvals(1) = numdis
   dimnames(1) = 'numdis'
   ivarid(2) = iarr_distmp
   nrank(2) = 1
   iddims(2,1) = 1

!
!2.9 Sediment model files
!------------------------
!

CASE (io_sedobc)
   IF (ifil.EQ.1) THEN
      CALL set_sedvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts)
   ELSE
      CALL set_sedvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts,&
                          & noprofs=noprofs)
   ENDIF
   GOTO 1000
CASE (io_sednst)
   CALL set_sedvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts,&
                       & noprofs=noprofs)
   GOTO 1000
CASE (io_sedspc,io_darspc,io_dissed)
   CALL set_sedvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts)
   GOTO 1000
   
!
!2.10 Biological model files
!---------------------------
!

CASE (io_bioobc)
   IF (ifil.EQ.1) THEN
      CALL set_biovars_atts(iddesc,ifil,iotype,filepars,numvars,varatts)
   ELSE
      CALL set_biovars_atts(iddesc,ifil,iotype,filepars,numvars,varatts,&
                          & noprofs=noprofs)
   ENDIF
   GOTO 1000
CASE (io_bionst)
   CALL set_sedvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts,&
                       & noprofs=noprofs)
   GOTO 1000
CASE (io_parabs,io_parrel,io_parsur)
   CALL set_biovars_atts(iddesc,ifil,iotype,filepars,numvars,varatts)
   GOTO 1000
   
!
!2.11 Particle model files
!------------------------
!

CASE (io_parspc,io_pargrd,io_parphs,io_parcld)
   CALL set_partvars_atts(iddesc,ifil,iotype,filepars,numvars,varatts)
   GOTO 1000
   
END SELECT

!
!3. Store file attributes
!-------------------------
!

filepars%nodims = nodims
filepars%dimvals = dimvals
filepars%dimnames = dimnames

!
!4. Store variable attributes (first part)
!-----------------------------------------
!
!4.1 Sub-variable names
!----------------------
!

IF (filepars%subnameid.EQ.1) THEN
   varatts(1)%ivarid = iarr_subname
   varatts(1)%nrank = 2
   varatts(1)%kind_type = MERGE(char_type,char_NF90,filepars%form.NE.'N')
   varatts(1)%global_dims(1:2) = (/lenname,nosubnames/)
   varatts(1)%nosubnames = nosubnames
   varatts(1)%subname(1:nosubnames) = subname(1:nosubnames)
ENDIF

!
!4.2 Time coordinate
!-------------------
!

idtime = filepars%timeid
IF (nocoords.EQ.1) THEN
   varatts(idtime)%ivarid = iarr_time
   varatts(idtime)%nrank = 1
   varatts(idtime)%kind_type = MERGE(char_type,char_NF90,filepars%form.NE.'N')
   varatts(idtime)%global_dims(1) = lentime
ENDIF

!
!4.3 Other variables
!-------------------
!

ivar_430: DO ivar=idtime+1,numvars
   varatts(ivar)%ivarid = ivarid(ivar)
   varatts(ivar)%nrank = nrank(ivar)
   varatts(ivar)%iddims = iddims(ivar,:)
   varatts(ivar)%kind_type = kind_type(ivar)
   n_431: DO n=1,nrank(ivar)
      idim = iddims(ivar,n)
      varatts(ivar)%global_dims(n) = dimvals(idim)
   ENDDO n_431
ENDDO ivar_430

1000 CONTINUE

!
!5. Store variable attributes (second part)
!------------------------------------------
!
!5.1 Sub-variable names
!---------------------
!

IF (filepars%subnameid.EQ.1) THEN
   CALL inquire_var(iarr_subname,f90_name=varatts(1)%f90_name,&
                  & standard_name=varatts(1)%standard_name,&
                  & long_name=varatts(1)%long_name,units=varatts(1)%units,&
                  & format = varatts(1)%format)
ENDIF

!
!5.2 Time coordinate
!-------------------
!

idtime = filepars%timeid
IF (filepars%nocoords.EQ.1) THEN
   CALL inquire_var(iarr_time,f90_name=varatts(idtime)%f90_name,&
                  & standard_name=varatts(idtime)%standard_name,&
                  & long_name=varatts(idtime)%long_name,&
                  & units=varatts(idtime)%units,&
                  & format = varatts(idtime)%format)
ENDIF

!
!5.3 Other variables
!-------------------
!

ivar_530: DO ivar=idtime+1,numvars
   CALL inquire_var(varatts(ivar)%ivarid,f90_name=varatts(ivar)%f90_name,&
                  & standard_name=varatts(ivar)%standard_name,&
                  & long_name=varatts(ivar)%long_name,&
                  & units=varatts(ivar)%units,format=varatts(ivar)%format)
ENDDO ivar_530

IF (loglev2.GE.pglev) THEN
   IF (pglev.LT.10) THEN
      WRITE (iolog,logfmt1) REPEAT(' ',pglev-1), pglev, logexit
   ELSE
      WRITE (iolog,logfmt2) REPEAT(' ',pglev-1), pglev, logexit
   ENDIF
ENDIF
pglev = pglev - 1


RETURN

END SUBROUTINE set_modvars_atts

!========================================================================

SUBROUTINE set_modvars_init_flags_2d
!************************************************************************
!
! *set_modvars_init_flags_2d* Set flags for open boundary work space arrays in
!                             the final condition file (2-D case)
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.12.1
!
! Description -
!
! Module calls - any_value_in_list
!
!************************************************************************
!
USE gridpars
USE obconds
USE physpars
USE utility_routines, ONLY: any_value_in_list


IF (nobu.GT.0) THEN
   obc2uvatu_init = any_value_in_list(ityp2dobu,7,list=(/6,7,9,10,11,12,13/))
ELSE
   obc2uvatu_init = .FALSE.
ENDIF

IF (nobv.GT.0) THEN
   obc2uvatv_init = any_value_in_list(ityp2dobv,7,list=(/6,7,9,10,11,12,13/))
ELSE
   obc2uvatv_init = .FALSE.
ENDIF


RETURN

END SUBROUTINE set_modvars_init_flags_2d
  
!========================================================================

SUBROUTINE set_modvars_init_flags_3d(iddesc,itypobu,itypobv)
!************************************************************************
!
! *set_modvars_init_flags_3d* Set flags for open boundary work space arrays in
!                             the final condition file (3-D) case
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)modvars_routines.f90  V2.12.1
!
! Description -
!
! Module calls - any_value_in_list
!
!************************************************************************
!
USE iopars
USE gridpars  
USE physpars
USE switches
USE utility_routines, ONLY: any_value_in_list

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc
INTEGER, INTENT(IN), DIMENSION(nobu,1) :: itypobu
INTEGER, INTENT(IN), DIMENSION(nobv,1) :: itypobv

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*   INTEGER File id
!*itypobu*  INTEGER Type of U-open boundary condition
!*itypobv*  INTEGER Type of V-open boundary condition
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: uflag, vflag
INTEGER :: l1, l2

uflag = .FALSE.; vflag = .FALSE.

IF (iddesc.EQ.io_3uvobc) THEN
   l1 = 2; l2 = 4
ELSE
   l1 = 2; l2 = 3
ENDIF

IF (nobu.GT.0) THEN
   uflag = any_value_in_list(itypobu(:,1),2,list=(/l1,l2/))
ENDIF

IF (nobv.GT.0) THEN
   vflag = any_value_in_list(itypobv(:,1),2,list=(/l1,l2/))
ENDIF

SELECT CASE (iddesc)
   CASE (io_3uvobc)
      obc3uvatu_init = uflag.AND.iopt_obc_3D.EQ.1
      obc3uvatv_init = vflag.AND.iopt_obc_3D.EQ.1
   CASE (io_salobc)
      obcsalatu_init = uflag.AND.iopt_obc_sal.EQ.1.AND.iopt_sal.EQ.2
      obcsalatv_init = vflag.AND.iopt_obc_sal.EQ.1.AND.iopt_sal.EQ.2
   CASE (io_tmpobc)
      obctmpatu_init = uflag.AND.iopt_obc_temp.EQ.1.AND.iopt_temp.EQ.2
      obctmpatv_init = vflag.AND.iopt_obc_temp.EQ.1.AND.iopt_temp.EQ.2
END SELECT

      
RETURN

END SUBROUTINE set_modvars_init_flags_3d

END MODULE modvars_routines
