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
! limitations under the Licence.

!************************************************************************
!
! *Usrdef_Model* User-defined model setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.12
!
! $Date: 2022-01-03 14:13:59 +0100 (Mo, 03 Jan 2022) $
!
! $Revision: 1447 $
!
! Description - test case river
!
! Reference -
!
! Subroutines - usrdef_init_params, usrdef_mod_params, usrdef_grid,
!               usrdef_partition, usrdef_phsics, usrdef_1dsur_spec,
!               usrdef_2dobc_spec, usrdef_profobc_spec, usrdef_1dsur_data,
!               usrdef_2dobc_data, usrdef_profobc_data
!
!************************************************************************
!

!========================================================================

SUBROUTINE usrdef_init_params
!************************************************************************
!
! *usrdef_init_params* Define parameters for monitoring
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test case river
!
! Reference -
!
! Calling program - simulation_start
!
!************************************************************************
!
USE iopars
USE paralpars

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: iproc


!
!1. Cold/warm start
!------------------
!

IF (ciffile%status.EQ.'W') then
   write(*,*) "It s a cold start"
   cold_start = .TRUE.
else
   write(*,*) "It s a warm start"
endif   

!
!3. Log files
!------------
!
!---program leveling in log files
iproc_310: DO iproc=1,npworld
   levprocs_ini(iproc) = 3
   levprocs_run(iproc) = 3
ENDDO iproc_310

!
!6. Timing
!---------
!

levtimer = 3

!
!7. Parallel setup
!-----------------
!

IF (npworld.GT.1) nprocscoh = 9


RETURN

END SUBROUTINE usrdef_init_params

!========================================================================

SUBROUTINE usrdef_mod_params
!************************************************************************
!
! *usrdef_mod_params* Define parameters for physical model
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.12
!
! Description - test case river
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE gridpars
USE iopars
USE paralpars
USE physpars
USE switches
USE syspars
USE tide
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: sflag

!
!*Local variables
!
CHARACTER (LEN=1) :: modform 
REAL              :: runid

procname(pglev+1) = 'usrdef_mod_params'
CALL log_timer_in()

modform = 'N'

!
!1. Setup flag
!-------------
!

!AC 280922
!sflag = MERGE(.TRUE.,.FALSE.,runtitle(6:6).EQ.'0')

!
!2. Process numbers
!------------------
!

IF (npworld.GT.1) THEN
   nprocsx = 3
!  ---number of processes in Y-direction
   nprocsy = 3
ENDIF

!
!3. Switches
!-----------
!
!---Cartesian or spherical (0/1)
iopt_grid_sph = 1

!AC 280922
!IF (.NOT.sflag) THEN
!  ---equation of state (0/1/2)
iopt_dens = 1
!  ---formulation for baroclinic pressure gradient (0/1/2)
iopt_dens_grad = 1
!ENDIF

!---salinity equation (0/1)
iopt_sal  = 2

!---temperature equation
iopt_temp = 2   
iopt_temp_optic = 1

!---advection scheme for scalars (0/1/2/3/4)
iopt_adv_scal = 3
iopt_adv_turb = 0 ! 0:disabled -  1:Upwind - 2: Lax-Wendroff - 3:TVD
iopt_adv_tvd  = 2 ! 2:Monotone limiter - 1:SuperBee

!---advection scheme for 2-D/3-D currents (0,1,2,3,4)
!IF (runtitle(7:7).EQ.'D') THEN
iopt_adv_2D = 3
iopt_adv_3D = 3
!ELSE
!   iopt_adv_2D = 1
!   iopt_adv_3D = 1
!ENDIF


!---horizontal diffusion scheme
iopt_hdif_coef = 2
iopt_hdif_2D = 1
iopt_hdif_3D = 1
iopt_hdif_turb = 0

!---type of vertical diffusion scheme (0/1/2/3)
iopt_vdif_coef = 3

!---type of algebraic turbulence closure (1/2/3/4/5)
iopt_turb_alg = 1

!---type of background mixing scheme (0/1/2)
iopt_turb_iwlim = 0

!---2-D mode o.b.c (0/1)
iopt_obc_2D = 1 !1
iopt_obc_3D = 0

!---nodal corrections for tides
iopt_astro_pars = 1
!
!4. Date/time parameters
!-----------------------
!
!---Start/End date (YYYY/MM/DD HH:MM:SS,mmm)
CStartDateTime(1:19) = '2007/01/01;00:00:00'
CEndDateTime(1:19)   = '2008/01/01;00:00:00'

!---time step
read( runtitle(6:6),*) runid

timestep = 4.0

!IF (runtitle(6:6).EQ.'1') then
!   timestep = 30.0 !15.0
!else if (runtitle(6:6).EQ.'2') then
!   timestep = 15.0 !15.0   
!else if (runtitle(6:6).EQ.'3') then
!   timestep = 5.0 !15.0   
!else if (runtitle(6:6).EQ.'4') then
!   timestep = 2.0 !15.0
!endif


!---counter for 3-D mode
ic3d = 10

!
!5. Physical model constants
!---------------------------
!
!---grid dimensions
nc = 321; nr = 321; nz = 15

!---number of river/open sea boundaries
nosbu = 2*(nr-1)
nosbv = 2*(nc-1)

dlat_ref= 50.0
dlon_ref= 4.0

!---number of tidal constituents
nconobc = 1

!---acceleration of gravity [m/s^2]
gacc_ref = 9.81

!---uniform mean water depth
depmean_cst = 20.0

!---uniform bottom roughness length
zrough_cst = 0.006

!---tidal indices
index_obc(1) = icon_M2 ! icon_M2

!
!6. Model I/O file properties
!----------------------------
!
!6.1 Input
!---------
!
!---model grid
modfiles(io_modgrd,1,1)%status = 'N'

!---initial conditions
modfiles(io_inicon,ics_phys,1)%status = 'N'
modfiles(io_inicon,1,1)%form = 'N'   

!---open boundary conditions (2-D)
modfiles(io_2uvobc,1,1)%status = 'N'


! METEO
iopt_meteo        = 1
iopt_meteo_data   = 1
iopt_meteo_stres  = 1
iopt_meteo_heat   = 1
iopt_meteo_precip = 0

iopt_sflux_qshort = 1

! iopt_obc_invbar = 1
 
!! Meteo files (ERA5) !!
IF (iopt_meteo.EQ.1)THEN
!   WRITE (cyear,'(I4.4)') iyear
   modfiles(io_metsur,1,1)%status = 'N'
   modfiles(io_metsur,1,1)%form   = 'N'
   modfiles(io_metsur,1,1)%filename = '/home/ulg/mast/acapet/Coherens_Forcings/Shading_2007/BCZ_2007.nc' !/home/acapet/Shading_test/BCZ_2007.nc'
   modfiles(io_metsur,1,1)%tlims = (/0,int_fill,450/)

!---meteo grid (ecmwf)                                                                                                                                                                                            
   surfacegrids(igrd_meteo,1)%nhtype = 1
   surfacegrids(igrd_meteo,1)%n1dat = 3
   surfacegrids(igrd_meteo,1)%n2dat = 3
   surfacegrids(igrd_meteo,1)%x0dat = 2.0
   surfacegrids(igrd_meteo,1)%y0dat = 51.1
   surfacegrids(igrd_meteo,1)%delxdat = 0.25
   surfacegrids(igrd_meteo,1)%delydat = 0.25
ENDIF



!6.2 Output
!----------
!
!---initial conditions
modfiles(io_fincon,ics_phys,2)%status = 'W'
modfiles(io_fincon,ics_phys,2)%form = modform
!---restart times
!AC2809022
!IF (sflag) THEN
norestarts = 2
ntrestart(1:2) = (/0,int_fill/)
!ELSE
!   norestarts = 0
!ENDIF

!
!7. Surface grid parameters
!--------------------------

! For Orthogonal coordinates
!surfacegrids(igrd_model,1)%delxdat = 125.0 
!surfacegrids(igrd_model,1)%delydat = 125.0
!surfacegrids(igrd_model,1)%x0dat = 0.0
!surfacegrids(igrd_model,1)%y0dat = 0.0
! For Spherical coordinates
surfacegrids(igrd_model,1)%delxdat = 0.0017884 
surfacegrids(igrd_model,1)%delydat = 0.0011236
surfacegrids(igrd_model,1)%x0dat = 2.05
surfacegrids(igrd_model,1)%y0dat = 51.15

!
!8. User output
!--------------


!---user output (0/1)
!AC 280922
iopt_out_tsers = 1
iopt_out_avrgd = 0

!---title for forcing files
intitle = runtitle(1:6)//runtitle(7:7)
!---title for user-output files
outtitle = runtitle(1:6)//runtitle(7:7)
!---number of output files
nosetstsr = 1 !MERGE(2,1,iopt_verif.EQ.0)
novarstsr = 4

!nosetsavr = 2
!novarsavr = 7

!nosetsanal = 1
!nofreqsanal = 1
!novarsanal = 4


CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_mod_params

!========================================================================

SUBROUTINE usrdef_grid
!************************************************************************
!
! *usrdef_grid* Define model grid arrays
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description - test case river
!
! Reference -
!
! Calling program - initialise_model
!
!************************************************************************
!
USE grid
USE gridpars
USE iopars
USE time_routines, ONLY: log_timer_in, log_timer_out


procname(pglev+1) = 'usrdef_grid'
CALL log_timer_in()

!
!3. Open boundary locations
!--------------------------
!
!---U-nodes
iobu(1:(nr-1))      = (/(1   ,l=1,nr-1)/)!   1 !(/1,nc/)
iobu(nr: 2*(nr-1))  = (/(nc  ,l=1,nr-1)/)

jobu(1:(nr-1))      = (/(l   ,l=1,nr-1)/)
jobu(nr:2*(nr-1))   = (/(l   ,l=1,nr-1)/)

!---V-nodes
iobv(1:(nc-1))      = (/(l   ,l=1,nc-1)/)!   1 !(/1,nc/)
iobv(nc:2*(nc-1))   = (/(l   ,l=1,nc-1)/)

jobv(1:(nc-1))      = (/(1   ,l=1,nc-1)/)
jobv(nc:2*(nc-1))   = (/(nr  ,l=1,nc-1)/)



CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_grid

!========================================================================

SUBROUTINE usrdef_partition
!************************************************************************
!
! *usrdef_partition* Define domain decomposition
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - domain_decomposition
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_partition

!========================================================================

SUBROUTINE usrdef_phsics
!************************************************************************
!
! *usrdef_phsics* Define arrays initial conditions for physics
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.1.1
!
! Description - test case river
!
! Reference -
!
! Calling program - initialise_model
!
! External calls - read_phsics
!
! Module calls - exchange_mod, file_suffix
!
!************************************************************************
!
USE density
USE grid
USE gridpars
USE iopars
USE modids
USE switches
USE paral_comms, ONLY: exchange_mod
USE time_routines, ONLY: log_timer_in, log_timer_out
USE utility_routines, ONLY : file_suffix

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j
REAL :: s, sl, sr, xc, xl, xr
INTEGER, DIMENSION(3) :: lbounds
INTEGER, DIMENSION(4) :: nhexch


IF (runtitle(6:6).EQ.'0') RETURN

procname(pglev+1) =  'usrdef_phsics'
CALL log_timer_in()

!
!1. Read initial conditions from spin-up phase
!---------------------------------------------
!

!modfiles(io_inicon,1,1)%filename = TRIM(intitle)//'.phsfin.'//&
!     & TRIM(file_suffix(modfiles(io_inicon,1,1)%form))
!iopt_out_anal = 0
!CALL read_phsics
!iopt_out_anal = 1

!
!2. Salinity
!-----------
!

sl = 35.0; sr = 35.0
xl = 25000.0; xr = 45000.0

i_210: DO i=1,ncloc
j_210: DO j=1,nrloc
   IF (maskatc_int(i,j)) THEN
      xc = 0.5*(gxcoord(i,j)+gxcoord(i+1,j))
      IF (xc.LT.xl) THEN
         s = sl
      ELSEIF (xc.GT.xr) THEN
         s = sr
      ELSE
         s = (sr-sl)*(xc-xl)/(xr-xl)+sl
      ENDIF
      sal(i,j,:) = s
   ENDIF
ENDDO j_210
ENDDO i_210

!
!3. Exchange halos
!-----------------
!

IF (iopt_MPI.EQ.1) THEN
   lbounds = (/1-nhalo,1-nhalo,1/); nhexch = nhdens
   CALL exchange_mod(sal,lbounds,nhexch,iarr_sal)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE usrdef_phsics

!========================================================================

SUBROUTINE usrdef_1dsur_spec
!************************************************************************
!
! *usrdef_1dsur_spec* Define specifier arrays for 1-D forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_1dsur_spec
!
!************************************************************************
!

IMPLICIT NONE


RETURN

END SUBROUTINE usrdef_1dsur_spec

!========================================================================

SUBROUTINE usrdef_2dobc_spec(iddesc)
!************************************************************************
!
! *usrdef_2dobc_spec* Define specifier arrays for 2-D open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description - test case river
!
! Reference -
!
! Calling program - define_2dobc_spec
!
!************************************************************************
!
USE gridpars
USE iopars
USE obconds
USE physpars
USE syspars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
!
!------------------------------------------------------------------------------
!
!*Local variables
!
REAL :: amp, crad, phase 


procname(pglev+1) = 'usrdef_2dobc_spec'
CALL log_timer_in()

!---type of conditions at open boundaries
ityp2dobu(1:(nr-1))    = 11
ityp2dobu(nr:2*(nr-1)) = 13

!---location elevation point
iloczobu = 1; iloczobv = 1

!---amplitudes
crad = SQRT(gacc_mean*depmean_cst)
amp = 0.8
phase = -halfpi
udatobu_amp(1:(nr-1),1) = crad*amp
zdatobu_amp(1:(nr-1),1) = amp

!---phases
udatobu_pha(1:(nr-1),1) = phase
zdatobu_pha(1:(nr-1),1) = phase
CALL log_timer_out()

RETURN

END SUBROUTINE usrdef_2dobc_spec

!========================================================================

SUBROUTINE usrdef_profobc_spec(iddesc,itypobux,itypobvy,iprofobux,iprofobvy,&
                             & noprofsd,indexprof,nofiles,nobux,nobvy,novars)
!************************************************************************
!
! *usrdef_profobc_spec* Define specifier arrays for open boundary conditions
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.12
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_profobc_spec
!
!************************************************************************
!

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, nobux, nobvy, nofiles, novars
INTEGER, INTENT(INOUT), DIMENSION(2:nofiles,novars) :: noprofsd
INTEGER, INTENT(INOUT), DIMENSION(nobux,novars) :: iprofobux, itypobux
INTEGER, INTENT(INOUT), DIMENSION(nobvy,novars) :: iprofobvy, itypobvy
INTEGER, INTENT(INOUT), DIMENSION(nobux+nobvy,2:nofiles,novars) :: indexprof

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id
!*itypobux*  INTEGER Type of U- or X-open boundary condition
!*itypobvy*  INTEGER Type of V- or Y-open boundary condition
!*iprofobux* INTEGER Profile numbers at U- or X-open boundaries
!*iprofobvy* INTEGER Profile numbers at V- or Y-open boundaries
!*noprofsd*  INTEGER Number of profiles per data file
!*indexprof* INTEGER Mapping array (for each data variable) of the profile
!                    numbers in the data files to the profile numbers assigned
!                    to the open boundaries. The physical size of the first
!                    dimension equals the number of profiles in a data file.
!*nofiles*   INTEGER Number of data files (+1)
!*nobux*     INTEGER Number of nodes at U- or X-open boundaries
!*nobvy*     INTEGER Number of nodes at V- or Y-open boundaries
!*novars*    INTEGER Number of data variables
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_profobc_spec

!========================================================================

SUBROUTINE usrdef_1dsur_data(ciodatetime,data1d,novars)
!************************************************************************
!
! *usrdef_1dsur_data* Define data for 1-D forcing
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_1dsur_data
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: novars
REAL, INTENT(INOUT), DIMENSION(novars) :: data1d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*ciodatetime* CHAR    Date/time in data file
!*data1d*      REAL    Surface data
!*novars*      INTEGER Number of surface data
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_1dsur_data

!========================================================================

SUBROUTINE usrdef_2dobc_data(iddesc,ifil,ciodatetime,data2d,nodat,novars)
!************************************************************************
!
! *usrdef_2dobc_data* Define open boundary data for 2-D mode
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.7
!
! Description -
!
! Reference -
!
! Calling program - define_2dobc_data
!
!************************************************************************
!
USE syspars

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: data2d

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Data file id (io_2uvobc or io_2xyobc)
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*data2d*      REAL    Input data
!*nodat*       INTEGER Number of input data
!*novars*      INTEGER Number of input parameters
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_2dobc_data

!========================================================================

SUBROUTINE usrdef_profobc_data(iddesc,ifil,ciodatetime,psiprofdat,&
                             & numprofs,nobcvars)
!************************************************************************
!
! *usrdef_profobc_data* Define physical open boundary profiles
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Model.f90  V2.12
!
! Description - empty default file
!
! Reference -
!
! Calling program - define_profobc_data
!
!************************************************************************
!
USE gridpars
USE syspars

IMPLICIT NONE

!
!*  Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nobcvars, numprofs
REAL, INTENT(INOUT), DIMENSION(numprofs,nz,nobcvars) :: psiprofdat

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*psiprofdat*  REAL     Profile arrays
!*numprofs*    INTEGER  Number of profiles in data file
!*nobcvars*    INTEGER  Effective number of data variables for which open
!                       boundary conditions are applied
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_profobc_data
