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

MODULE iopars
!************************************************************************
!
! *iopars* Parameters for input/output
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)iopars.f90  V2.x
!
!
! $Date: 2021-11-28 10:15:48 +0100 (So, 28 Nov 2021) $
!
! $Revision: 1435 $
!
! Description - 
!
!************************************************************************
!
USE datatypes
USE syspars

IMPLICIT NONE


!
!1. General parameters
!---------------------
!

LOGICAL :: cold_start = .FALSE., next_simul = .FALSE.
INTEGER :: isimul = 0, nopenf = 0, nosimul = 0
CHARACTER (LEN=lentitle) :: intitle, outtitle, runtitle
CHARACTER (LEN=lenname), DIMENSION(MaxProgLevels) :: procname

!
!2. Model input/output
!---------------------
!
!2.1 Attributes
!--------------
!

INTEGER :: numglbatts = 7
INTEGER, DIMENSION(MaxIOTypes,2) :: maxdatafiles
CHARACTER (LEN=lendesc) ::  Conventions_CF, history_CF, institution_CF, &
                          & references_CF, source_CF
TYPE (FileParams), DIMENSION(MaxIOTypes,MaxIOFiles,2) :: modfiles
TYPE (GridParams), DIMENSION(MaxGridTypes,MaxGridFiles) :: surfacegrids
TYPE (GlbFileAtts), DIMENSION(MaxGlbAtts) :: glbatts

!
!2.2 File descriptor key ids
!---------------------------
!
!---domain decomposition
INTEGER, PARAMETER :: io_mppmod = 1
!---initial conditions
INTEGER, PARAMETER :: io_inicon = 2, io_fincon = 3
!---model grid and data locations
INTEGER, PARAMETER :: io_modgrd = 4
!---external grids (absolute coordinates)
INTEGER, PARAMETER :: io_metabs = 5, io_sstabs = 6, io_wavabs = 7, &
                    & io_parabs = 8, io_nstgrd = 9
!---external grids (relative coordinates)
INTEGER, PARAMETER :: io_metrel = 10, io_sstrel = 11, io_wavrel = 12, &
                    & io_parrel = 13
!---specifiers for external modules
INTEGER, PARAMETER :: io_sedspc = 14, io_darspc = 15, io_parspc = 16
!---open (surface) boundaries
INTEGER, PARAMETER :: io_1uvsur = 17, io_2uvobc = 18, io_2xyobc = 19, &
                    & io_3uvobc = 20, io_3xyobc = 21, io_salobc = 22, &
                    & io_tmpobc = 23, io_sedobc = 24, io_bioobc = 25
!---nested output
INTEGER, PARAMETER :: io_nstspc = 26, io_2uvnst = 27, io_2xynst = 28, &
                    & io_3uvnst = 29, io_3xynst = 30, io_salnst = 31, &
                    & io_tmpnst = 32, io_sednst = 33, io_bionst = 34
!---surface data
INTEGER, PARAMETER :: io_metsur = 35, io_sstsur = 36, io_wavsur = 37, &
                    & io_parsur = 38
!---structures
INTEGER, PARAMETER :: io_drycel = 39, io_thndam = 40, io_weibar = 41
!---discharges
INTEGER, PARAMETER :: io_disspc = 42, io_disloc = 43, io_disvol = 44 , &
                    & io_discur = 45, io_dissal = 46, io_dissed = 47, &
                    & io_distmp = 48
!--data for particle model
INTEGER, PARAMETER :: io_parcld = 49, io_pargrd = 50, io_parphs = 51

!---key ids for surface data grids
INTEGER, PARAMETER :: igrd_model = 1, igrd_meteo = 2, igrd_sst = 3, &
                    & igrd_waves = 4, igrd_par = 5

!---file descriptor in string format
CHARACTER (LEN=6), DIMENSION(MaxIOTypes) :: modfiles_desc = &
&  (/'mppmod','inicon','fincon','modgrd','metabs','sstabs','wavabs','bioabs',&
&    'nstgrd','metrel','sstrel','wavrel','biorel','sedspc','darspc','parspc',&
&    '1uvsur','2uvobc','2xyobc','3uvobc','3xyobc','salobc','tmpobc','sedobc',&
&    'bioobc','nstspc','2uvnst','2xynst','3uvnst','3xynst','salnst','tmpnst',&
&    'sednst','bionst','metsur','sstsur','wavsur','biosur','drycel','thndam',&
&    'weibar','disspc','disloc','disvol','discur','dissal','dissed','distmp',&
&    'parcld','pargrd','parphs'/)

!---key ids for initial conditions
INTEGER, PARAMETER :: ics_phys = 1, ics_sed = 2, ics_morph = 3, ics_part = 4, &
                    & ics_bio = 5

!
!2.3 CIF parameters
!------------------
!
INTEGER :: ciflinenum
TYPE (FileParams) :: ciffile
LOGICAL, DIMENSION(MaxCIFTypes) :: cifstatus = .FALSE.

!---key ids for CIF
INTEGER, PARAMETER :: icif_anal = 1, icif_avr = 2, icif_bio = 3, &
                    & icif_dar = 4, icif_def = 5, icif_disspc = 6, &
                    & icif_freq = 7, icif_mod = 8, icif_mon = 9, &
                    & icif_morph = 10, icif_nstspc = 11, icif_part = 12, &
                    & icif_partspc = 13, icif_sed = 14, icif_sedspc = 15, &
                    & icif_tspart = 16, icif_tsr = 17

!---CIF block names
CHARACTER (LEN=lencifblock), DIMENSION(MaxCIFTypes) :: cifblocks = &
& (/'anal_params                         ',&
&   'avr_params                          ',&
&   'bio_params                          ',&
&   'dar_params                          ',&
&   'defruns                             ',&
&   'dischr_spec                         ',&
&   'anal_freqs                          ',&
&   'mod_params                          ',&
&   'init_params                         ',&
&   'morph_params                        ',&
&   'nstgrd_spec                         ',&
&   'part_params                         ',&
&   'part_spec                           ',&
&   'sed_params                          ',&
&   'sed_spec                            ',&
&   'pout_params                         ',&
&   'tsr_params                          '/)

!
!3. User output
!--------------
!
!---fill and minimum values
INTEGER :: int_fill
REAL :: float_fill, float_fill_eps, real_fill, real_fill_eps
REAL (KIND=kndrlong) :: double_fill, double_fill_eps

!---ASCII formats
LOGICAL :: freefmt
CHARACTER (Len=lenformat) :: characterfmt, doublefmt, integerfmt, logicalfmt, &
                           & realfmt  

!---time series
INTEGER :: nosetstsr = 0, nostatstsr = 0,  novarstsr = 0
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ivarstsr, lstatstsr
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: weightstatstsr

TYPE (FileParams), ALLOCATABLE, DIMENSION(:) :: tsr0d, tsr2d, tsr3d
TYPE (OutGridParams), ALLOCATABLE, DIMENSION(:) :: tsrgpars
TYPE (StationLocs), ALLOCATABLE, DIMENSION(:) :: tsrstatlocs
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: tsrvars

!---time averages
INTEGER :: nosetsavr = 0, nostatsavr = 0, novarsavr = 0
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ivarsavr, lstatsavr
TYPE (FileParams), ALLOCATABLE, DIMENSION(:) :: avr0d, avr2d, avr3d
TYPE (OutGridParams), ALLOCATABLE, DIMENSION(:) :: avrgpars
TYPE (StationLocs), ALLOCATABLE, DIMENSION(:) :: avrstatlocs
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: avrvars

!---harmonic analysis
INTEGER :: nosetsanal = 0, nofreqsanal = 0, nostatsanal = 0, novarsanal = 0
CHARACTER (LEN=lenfreq), ALLOCATABLE, DIMENSION(:,:) :: harm_freq_names
INTEGER, ALLOCATABLE, DIMENSION(:) :: index_anal, nofreqsharm
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ifreqsharm, ivarsanal, lstatsanal
REAL, ALLOCATABLE, DIMENSION(:) :: fnodal_anal, phase_anal
REAL, ALLOCATABLE, DIMENSION(:,:) :: harm_freqs 
TYPE (FileParams), ALLOCATABLE, DIMENSION(:) :: res0d, res2d, res3d
TYPE (FileParams), ALLOCATABLE, DIMENSION(:,:) :: amp0d, amp2d, amp3d, pha0d, &
                                                & pha2d, pha3d, ell2d, ell3d
TYPE (OutGridParams), ALLOCATABLE, DIMENSION(:) :: analgpars
TYPE (StationLocs), ALLOCATABLE, DIMENSION(:) :: analstatlocs
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: analvars

!---elliptic parameters
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ivarsell, ivecell2d, ivecell3d
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: ellvars

!---output operators
INTEGER, PARAMETER :: oopt_null = 0, oopt_dep = 1, oopt_int = 2, &
                    & oopt_klev = 3, oopt_min = 4, oopt_max = 5, oopt_mean = 6

!
!4. Monitoring
!-------------
!
!4.1 log files
!-------------
!

LOGICAL :: exitlog
CHARACTER (LEN=*), PARAMETER :: logfmt1 = '(A,I1,'':'',A)', &
                              & logfmt2 = '(A,I2,'':'',A)'
CHARACTER (LEN=1), PARAMETER :: logexit = 'R'
CHARACTER (len=leniofile) :: inilog_file, runlog_file
INTEGER :: iolog = 0, loglev1, loglev2, pglev, runlog_count
INTEGER, ALLOCATABLE, DIMENSION(:) :: levprocs_ini, levprocs_run

!
!4.2 Error files
!---------------
!
!---parameters
LOGICAL :: errchk
CHARACTER (len=leniofile) :: errlog_file
INTEGER :: errstat, ioerr = 0, maxerrors, nerrs = 0
INTEGER, ALLOCATABLE, DIMENSION(:) :: levprocs_err

!---key ids
INTEGER, PARAMETER :: ierrno_fopen = 1, ierrno_fclose = 2, ierrno_read = 3,&
                    & ierrno_write = 4, ierrno_fend = 5, ierrno_input = 6,&
                    & ierrno_inival = 7, ierrno_runval = 8, ierrno_alloc = 9,&
                    & ierrno_arg = 10, ierrno_waves = 11, ierrno_comms = 12, &
                    & ierrno_MPI = 13, ierrno_CDF = 14, ierrno_MCT = 15

!---messages
CHARACTER (LEN=lenerrcode), PARAMETER, DIMENSION(MaxErrCodes) :: error_code = &
      & (/'Not possible to open file                            ',&
      &   'Unable to close file                                 ',&
      &   'Read error                                           ',&
      &   'Write error                                          ',&
      &   'End of file condition                                ',&
      &   'Wrong input values                                   ',&
      &   'Invalid initial values for model parameters or arrays',&
      &   'Invalid values for variables at run time             ',&
      &   'Not possible to allocate arrays                      ',&
      &   'Missing or invalid argument in routine call          ',&
      &   'Error in wave model                                  ',&   
      &   'Communication error                                  ',&
      &   'Error in MPI call                                    ',&
      &   'Error in NetCDF call                                 ',&
      &   'Error in MCT call                                    '/)
!
!4.3 Warning file
!----------------
!

LOGICAL :: warnflag, warning
CHARACTER (len=leniofile) :: warlog_file
INTEGER :: iowarn = 0

!
!4.4 Monitoring files
!--------------------
!
!---implicit scheme
LOGICAL :: monflag
CHARACTER (len=leniofile) :: monlog_file
INTEGER :: iomon = 0, monlog = 0

!---sediment modules
LOGICAL :: sedflag, sedlog, sedlog_date
CHARACTER (LEN=leniofile) :: sedlog_file
INTEGER :: iosed 

!
!4.5 Timer
!---------
!
!---parameters
LOGICAL :: timer = .FALSE.
CHARACTER (len=leniofile) :: timing_file
INTEGER :: levtimer = 0, maxwaitsecs = 3600, nowaitsecs = 0, npcc_max, &
         & npcc_rate, timer_format = 1
INTEGER (KIND=kndilong), DIMENSION(MaxTimers) :: nopcc

!---key-ids
INTEGER, PARAMETER :: &
  & itm_hydro = 1, itm_1dmode = 2, itm_2dmode = 3, itm_3dmode = 4, &
  & itm_dens = 5, itm_temp = 6, itm_sal = 7, itm_init = 8, &
  & itm_trans = 9, itm_adv = 10, itm_hdif = 11, itm_vdif = 12, &
  & itm_mg = 13, itm_phgrad = 14, itm_input = 15, itm_output = 16, &
  & itm_inout = 17, itm_com_coll = 18, itm_com_comb = 19, itm_com_copy = 20, &
  & itm_com_dist = 21, itm_com_exch = 22, itm_com_util = 23, itm_coms = 24, &
  & itm_MPI = 25, itm_CDF = 26, itm_MCT = 27, itm_arrint = 28, itm_user = 29, &
  & itm_nest = 30, itm_libs = 31, itm_astro = 32, itm_bconds = 33, &
  & itm_meteo = 34, itm_waves = 35, itm_wavemod = 36, itm_structs = 37, &
  & itm_wait = 38, itm_sed = 39, itm_bio = 40, itm_morph = 41, itm_part = 42


!---timer descriptions
CHARACTER (LEN=20), DIMENSION(MaxTimers) :: desctimer = &
   & (/'Hydrodynamics       ', '1D mode             ', &
     & '2D mode             ', '3D mode             ', &
     & 'Density             ', 'Temperature         ', &
     & 'Salinity            ', 'Initialisation      ', &
     & 'Transport           ', 'Advection           ', &
     & 'Horizontal diffusion', 'Vertical diffusion  ', &
     & 'Multigrid           ', 'Baroclinic pressure ', &
     & 'Input               ', 'Output              ', &
     & 'Input/output        ', 'Collect comms       ', &
     & 'Combine comms       ', 'Copy comms          ', &
     & 'Distribute comms    ', 'Exchange comms      ', &
     & 'Utility comms       ', 'Parallel comms      ', &
     & 'MPI calls           ', 'netCDF calls        ', &
     & 'MCT calls           ', 'Array interpolation ', &
     & 'User calls          ', 'Nesting procedures  ', &
     & 'Library calls       ', 'Astronomical tide   ', &
     & 'Boundary conditions ', 'Meteo               ', &
     & 'Surface waves       ', 'Wave model          ', &
     & 'Structures          ', 'Wait calls          ', &
     & 'Sediment            ', 'Biology             ', &
     & 'Morphology          ', 'Particle model      '/)

!
!4.6  File for dredging and relocation
!-------------------------------------
!

CHARACTER (len=leniofile) :: darlog_file

!
!5. netCDF parameters
!--------------------
!

INTEGER :: char_NF90 = 0, clobber_NF90 = 0, double_NF90 = 0, fill_NF90 = 0, &
         & fill_int_NF90 = 0, global_NF90 = 0, int_NF90 = 0, log_NF90 = 0, &
         & noerr_NF90 = 0,  nofill_NF90 = 0, nowrite_NF90 = 0, &
         & offset_64bit_NF90 = 0, real_NF90 = 0, share_NF90 = 0, &
         & sizehint_default_NF90 = 0, unlimited_NF90 = 0, write_NF90 = 0
REAL :: fill_real_NF90 = 0.0
REAL (KIND=kndrlong):: fill_double_NF90 = 0.0_kndrlong

!
!6. Work space parameters
!------------------------
!

LOGICAL :: inflagdis, inflagobc2d, inflagobc3d, inflagobc4d, &
         & outflagdis, outflagobc2d, outflagobc3d, outflagobc4d, &
         & outflagnst2d, outflagnst3d,  outflagnst4d, outflagnstbio

SAVE

!
! Name            Type     Purpose
!------------------------------------------------------------------------------
! general parameters
!*cold_start*     LOGICAL If .TRUE., main program stops after initialisation
!*next_simul*     LOGICAL .TRUE. to start next simulation, .FALSE. to exit 
!                         program
!*isimul*         INTEGER Simulation number
!*nopenf*         INTEGER Number of connected files
!*nosimul*        INTEGER Total number of simulations defined in the defruns
!                         file
!*intitle*        CHAR    Title for model forcing files
!*outtitle*       CHAR    Title for user output files
!*runtitle*       CHAR    Simulation title
!*procname*       CHAR    Name of subprogram at current and higher levels
!
! model I/O
!*numglbatts*     INTEGER Number of global attributes for files in COHERENS
!                         standard (CF-1.6) format
!*maxdatafiles*   INTEGER Largest file index for an active file of given type
!*ciflinenum*     INTEGER Number of the last input line read from a CIF file 
!*modfiles*       DERIVED Attributes of model I/O files
!                   first dimension  => file descriptor
!                   second dinebsion => file number
!                   third dimension  => input (1), output (2)
!*surfacegrids*   DERIVED Attributes of surface data grids
!*glbatts*        DERIVED Global attributes for files in COHERENS standard
!                         (CF-1.6) format
!*ciffile         DERIVED Attributes of the CIF file
!*ciflinenum*     INTEGER CIF line number
!*cifstatus*      LOGICAL .TRUE./.FALSE. if a CIF block is present in the CIF
!                         file
!
! key ids and descriptors
! io_?            INTEGER File descriptor key ids
!*io_mppmod*      INTEGER Domain decomposition
!*io_inicon*      INTEGER Initial condition files
!*io_fincon*      INTEGER Final condition files
!*io_modgrd*      INTEGER Model bathymetry
!*io_metabs*      INTEGER Meteo surface grid in absolute coordinates
!*io_sstabs*      INTEGER Surface SST grid in absolute coordinates
!*io_wavabs*      INTEGER Surface grid for wave data in absolute coordinates
!*io_parabs*      INTEGER Surface grid for PAR data in absolute
!                         coordinates
!*io_nstgrd*      INTEGER Locations of nested open boundaries
!*io_metrel*      INTEGER Meteo surface grid in relative coordinates
!*io_sstrel*      INTEGER Surface SST grid in relative coordinates
!*io_wavrel*      INTEGER Surface grid for wave data in relative coordinates
!*io_parrel*      INTEGER Surface grid for biological variables in relative
!                         coordinates
!*io_sedspc*      INTEGER Specifiers arrays for sediment model
!*io_darspc*      INTEGER Specifiers for dredging/relocation module
!*io_parspc*      INTEGER Specifiers for particle model
!*io_1uvsur*      INTEGER Surface boundary conditions/data for 1-D mode
!*io_2uvobc*      INTEGER Open boundary conditions/data for 2-D mode
!*io_2xyobc*      INTEGER Tangential open boundary conditions/data for 2-D mode
!*io_3uvobc*      INTEGER Open boundary conditions/data for 3-D mode
!*io_3xyobc*      INTEGER Tangential open boundary conditions/data for 3-D mode
!*io_salobc*      INTEGER Open boundary conditions/data for salinity
!*io_tmpobc*      INTEGER Open boundary conditions/data for temperature
!*io_sedobc*      INTEGER Open boundary conditions/data for sediment
!                         variables
!*io_bioobc*      INTEGER Open boundary conditions/data for biological
!                         variables
!*io_nstspc*      INTEGER Specifiers for nesting
!*io_2uvnst*      INTEGER 2-D mode output for nesting
!*io_2xynst*      INTEGER 2-D mode output for nesting for tangential
!*io_3uvnst*      INTEGER 3-D mode output for nesting
!*io_3xynst*      INTEGER 3-D mode output for nesting for tangential
!*io_salnst*      INTEGER Salinity output for nesting
!*io_tmpnst*      INTEGER Temperature output for nesting
!*io_sednst*      INTEGER Sediment output for nesting
!*io_bionst*      INTEGER Biological output for nesting
!*io_metsur*      INTEGER Surface meteo data
!*io_sstsur*      INTEGER Surface SST data
!*io_wavsur*      INTEGER Wave data
!*io_parsur*      INTEGER 2-D biological data
!*io_drycel*      INTEGER Dry cells
!*io_thndam*      INTEGER Thin dams
!*io_weibar*      INTEGER Weirs/barriers
!*io_disspc*      INTEGER Discharge specifier arrays
!*io_disloc*      INTEGER Discharge locations
!*io_disvol*      INTEGER Volume discharge data
!*io_discur*      INTEGER Discharge area and direction
!*io_dissal*      INTEGER Salinity discharge data
!*io_distmp*      INTEGER Temperature discharge data
!*io_dissed*      INTEGER Sediment discharge data
!*io_parcld*      INTEGER Particle cloud discharge locations
!*io_pargrd*      INTEGER Model grid arrays used in particle module
!*io_parphs*      INTEGER Physical data for the particle model
!
!*icif_?*         INTEGER Type of CIF blocks
!*icif_anal*      INTEGER parameters for harmonic analysis
!*icif_avr*       INTEGER time averaging output parameters
!*icif_bio*       INTEGER biological model setup parameters
!*icif_dar*       INTEGER dredging/relocation model setup parameters
!*icif_def*       INTEGER 'defruns' file
!*icif_disspc*    INTEGER specifiers for discharge module
!*icif_freq*      INTEGER frequencies for harmonic analysis
!*icif_mod*       INTEGER model setup parameters 
!*icif_mon*       INTEGER monitoring parameters
!*icif_morph*     INTEGER morphology model setup parameters
!*icif_nstspc*    INTEGER specifiers for nesting
!*icif_part*      INTEGER particle model setup parameters
!*icif_partspc*   INTEGER particle model specifiers
!*icif_sed*       INTEGER sediment model setup parameters
!*icif_sedspc*    INTEGER sediment model specifiers
!*icif_tspart*    INTEGER Particle trajectory output
!*icif_tsr*       INTEGER time series output parameters


! ics_?           INTEGER Initial condition files
!*ics_phys*       INTEGER Physical model
!*ics_bio*        INTEGER Biological model
!*ics_sed*        INTEGER Sediment (multi-fraction) transport model
!*ics_morph*      INTEGER Morphology model
!*ics_part*       INTEGER Particle model
!
! igrd_?          INTEGER Surface grid key ids
!*igrd_model*     INTEGER Model grid
!*igrd_meteo*     INTEGER Meteo grid
!*igrd_sst*       INTEGER SST grid
!*igrd_waves*     INTEGER Surface grid for wave data
!*igrd_par*       INTEGER Surface grid for surface PAR (biology)
!
!*double_fill*    DOUBLE  Flag for undefined/invalid/missing double real
!                         data or model variables
!*double_fill_eps*DOUBLE  =0.00001*ABS(real_double)
!*int_fill*       INTEGER Flag for undefined/invalid/missing integer
!                         data or model variables
!*longint_fill*   LONGINT Flag for undefined/invalid/missing long (64B) integer
!                         data or model variables
!*float_fill*      REAL/DOUBLE Flag for undefined/invalid/missing real
!                          data or model variables
!*float_fill_eps*  REAL/DOUBLE =0.00001*ABS(float_fill)
! ASCII formats
!*freefmt*         LOGICAL All FORTRAN output is written in free format
!                          if .TRUE.
!*characterfmt*    CHAR    FORTRAN output format for character data
!*doublefmt*       CHAR    FORTRAN output format for formatted double precision
!                          real data
!*integerfmt*      CHAR    FORTRAN output format for formatted integer data
!*logicalfmt*      CHAR    FORTRAN output format for formatted logical data
!*realfmt*         CHAR    FORTRAN output format for formatted simple precision
!                          real data
! time series output
!*nosetstsr*      INTEGER Number of output file series
!*nostatstsr*     INTEGER Number of stations for irregular output
!*novarstsr*      INTEGER Number of output variables
!*ivarstsr*       INTEGER Output variable indices per file index
!*lstatstsr*      INTEGER Station label per file index
!*tsr0d*          DERIVED Attributes 0-D output files
!*tsr2d*          DERIVED Attributes 2-D output files
!*tsr3d*          DERIVED Attributes 3-D output files
!*tsrgpars*       DERIVED Attributes of output grid
!*tsrstatlocs*    DERIVED Index positions of output stations
!*tsrvars*        DERIVED Output variable attributes
!
! time averaged output
!*nosetsavr*      INTEGER Number of output file series
!*nostatsavr*     INTEGER Number of stations for irregular output
!*novarsavr*      INTEGER Number of output variables
!*ivarsavr*       INTEGER Output variable indices per file index
!*lstatavr*       INTEGER Station label per file index
!*avr0d*          DERIVED Attributes 0-D output files
!*avr2d*          DERIVED Attributes 2-D output files
!*avr3d*          DERIVED Attributes 3-D output files
!*avrgpars*       DERIVED Attributes of output grid
!*avrstatlocs*    DERIVED Index positions of output stations
!*avrvars*        DERIVED Output variable attributes
!
! harmonic analysis
!*nosetsanal*     INTEGER Number of files for each specific harmonic output
!*nofreqsanal*    INTEGER Number of frequencies for harmonic analysis
!*nostatsanal*    INTEGER Number of stations for irregular output
!*novarsanal*     INTEGER Number of output variables
!*harm_freq*      REAL    Frequencies of the harmonic components     [radian/s]
!*harm_freq_names*CHAR    Names of the frequencies used in harmonic analysis
!*phase_anal*     REAL    Nodal factors of the tidal constituents used in the
!                         harmonic analysis procedure 
!*phase_anal*     REAL    Astronomical phases of the tidal cobnstituents used
!                         in the harmonic analysis procedure 
!*index_anal*     INTEGER Tidal constituent frequency indices (undefined if 0)
!*ifreqsharm*     INTEGER Tidal frequency key index used for each set
!*ivarsanal*      INTEGER Output variable indices per set
!*lstatanal*      INTEGER Station label per set
!*nofreqsharm*    INTEGER Number of frequencies per file series
!*res0d*          DERIVED Attributes 0-D output residual files
!*res2d*          DERIVED Attributes 2-D output residual files
!*res3d*          DERIVED Attributes 3-D output residual files
!*amp0d*          DERIVED Attributes 0-D output amplitude files
!*amp2d*          DERIVED Attributes 2-D output amplitude files
!*amp3d*          DERIVED Attributes 3-D output amplitude files
!*pha0d*          DERIVED Attributes 0-D output phase files
!*pha2d*          DERIVED Attributes 2-D output phase files
!*pha3d*          DERIVED Attributes 3-D output phase files
!*ell2d*          DERIVED Attributes 2-D elliptic output files
!*ell3d*          DERIVED Attributes 3-D elliptic output files
!*analgpars*      DERIVED Attributes of output grid
!*analstatlocs*   DERIVED Index positions of output stations
!*analvars*       DERIVED Attributes of analysed variables
!*ivarsell*       INTEGER Output variable indices for elliptic parameters per
!                         file index
!*ivecell2d*      INTEGER Variable indices of 2-D elliptic vector components
!*ivecell3d*      INTEGER Variable indices of 3-D elliptic vector components
!*ellvars*        DERIVED Attributes of elliptic variables
!
! operators for user output data
!*oopt_null*      INTEGER Null operator (i.e. no operator applied)
!*oopt_dep*       INTEGER Value at a prescribed vertical water depth
!*oopt_int*       INTEGER Integrated value
!*oopt_klev*      INTEGER Value at a prescribed vertical (C-node) grid level
!*oopt_min*       INTEGER Minimum value
!*oopt_max*       INTEGER Maximum value
!*oopt_mean*      INTEGER Mean value
!
! log files
!*exitlog*        LOGICAL Enables/disables an 'exit' message in log file just
!                         before RETURN statement
!*logfmt1*        CHAR    Format of log message
!*logfmt2*        CHAR    Format of log message
!*logexit*        CHAR    Exit code
!*inilog_file*    CHAR    Name of (generic) 'inilog' file
!*runlog_file*    CHAR    Name of (generic) 'runlog' file
!*loglev_def*     INTEGER Default number of program leveling in log files
!*iolog*          INTEGER File unit of log file
!*loglev1*        INTEGER Program leveling in log files
!*loglev2*        INTEGER Program leveling in log files (exit code)
!*pglev*          INTEGER Current program level
!*runlog_count*   INTEGER Determines the number of time steps after which the
!                         log-file will be re-written
!*levprocs_ini*   INTEGER Program leveling in 'inilog' files for each process
!*levprocs_run*   INTEGER Program leveling in 'runlog' files for each process
!
! error files
!*errchk*         LOGICAL Enables/disables error checking
!*errlog_file*    CHAR    Name of (generic) error file
!*errstat*        INTEGER Error status number as returned by FORTRAN, MPI and
!                         netcdf calls
!*ioerr*          INTEGER File unit of error file
!*maxerrors*      INTEGER Maximum allowed number of error messages
!*nerrs*          INTEGER Number of detected errors
!*levprocs_err*   INTEGER Level of error coding in 'errlog' files for each
!                         process
!                    = 0 => no error messages written
!                    = 1 => errors are writen during initialisation only
!                    = 2 => errors are written at any time during the run
! ierrno_?        INTEGER Error code key ids
!*ierrno_fopen*   INTEGER File opening error
!*ierrno_fclose*  INTEGER File closing error
!*ierrno_read*    INTEGER Read error
!*ierrno_write*   INTEGER Write error
!*ierrno_fend*    INTEGER End of file condition
!*ierrno_input*   INTEGER Input error
!*ierrno_inival*  INTEGER Invalid parameter or initial value
!*ierrno_runval*  INTEGER Invalid run time value
!*ierrno_alloc*   INTEGER Allocation error
!*ierrno_arg*     INTEGER Missing argument or invalid argument value
!*ierrno_comms*   INTEGER Communication error
!*ierrno_MPI*     INTEGER MPI error
!*ierrno_CDF*     INTEGER Netcdf error
!*error_code*     CHAR    List of error message codes
!
! warning file
!*warning*        LOGICAL Enables issuing of warning messages by master
!*warlog_file*    CHAR    Name of (generic) 'warning' file
!*iowarn*         INTEGER File unit of warning file
!
! monitoring for the implicit/multigrid scheme
!*monflag*        LOGICAL Enables issuing of monitoring messages for the
!                         implicit scheme 
!*monlog_file*    CHAR    Name of (generic) 'monitoring' file
!*monlog*         INTEGER Parameter to control level of monitoring
!                    = 0 => no monitoring
!                    = 1 => outer loop only
!                    = 2 => outer and inner loop (without smoothing iterations)
!                    = 3 => outer and inner loop (including smoothing
!                           iterations)
!*iomon*          INTEGER File unit of monitoring file
!
! sediment monitoring
!*sedmon*         LOGICAL  Enables issuing monitoring messages for the sediment
!                          module
!*sedlog_file*    CHAR     Name of the sediment monitoring file
!*iosed*          INTEGER  File unit of the sediment monitoring file
!
! timer
!*timer*          LOGICAL .TRUE. if levtimer > 0
!*timing_file*    CHAR    Name of (generic) timing file
!*levtimer*       INTEGER Parameter to control level of info in timing file
!                     = 0 => no timing
!                     = 1 => execution time only
!                     = 2 => as 1, including timing (in %) of different model
!                            components, in parallel case given by maximum,
!                            minimum and mean over all processes
!                     = 3 => as 2, in parallel case the timing (in %) is listed
!                            for all processes
!*maxwaitsecs*    INTEGER Maximum allowed time (in seconds) for program 
!                         suspension
!*nowaitsecs*     INTEGER Number of seconds to suspend process in wait call
!*npcc_max*       INTEGER Maximum clock count on process clock
!*npcc_rate*      INTEGER Number of process clock counts per second
!*timer_format*   INTEGER Switch to select for time formatting
!*nopcc*          LONGINT Number of process clock counts per timer
! itm_?           INTEGER Key ids of timer processes
!*itm_hydro*      INTEGER hydrodynamics
!*itm_1dmode*     INTEGER 1-D mode calculations
!*itm_2dmode*     INTEGER 2-D mode calculations
!*itm_3dmode*     INTEGER 3-D mode calculations
!*itm_dens*       INTEGER Density calculations
!*itm_temp*       INTEGER Temperature
!*itm_sal*        INTEGER Salinity
!*itm_init*       INTEGER Model initialisation
!*itm_trans*      INTEGER Transport modules
!*itm_adv*        INTEGER Advection
!*itm_hdif*       INTEGER Horizontal diffusion
!*itm_vdif*       INTEGER Vertical diffusion
!*itm_phgrad*     INTEGER Baroclonic pressure gradient
!*itm_input*      INTEGER Reading
!*itm_output*     INTEGER Writing
!*itm_inout*      INTEGER Reading and writing
!*itm_com_coll*   INTEGER Collect communications
!*itm_com_comb*   INTEGER Combine communications
!*itm_com_copy*   INTEGER Copy communications
!*itm_com_dist*   INTEGER Distribute communications
!*itm_com_exch*   INTEGER Exchange communications
!*itm_com_util*   INTEGER Utility communications
!*itm_coms*       INTEGER All comminications
!*itm_MPI*        INTEGER MPI routines
!*itm_CDF*        INTEGER Netcdf routines
!*itm_arrint*     INTEGER Array interpolation
!*itm_user*       INTEGER User output
!*itm_nest*       INTEGER Nesting
!*itm_libs*       INTEGER Libraries
!*itm_astro*      INTEGER Astronomical forcing
!*itm_bconds*     INTEGER Boundary conditions
!*itm_meteo*      INTEGER Meteo
!*itm_structs*    INTEGER Structures
!*itm_wait*       INTEGER Wait calls
!*itm_biolgy*     INTEGER Biology
!*itm_sed*        INTEGER Sediment
!*desctimer*      INTEGER Desciptions of timers
!
!************************************************************************
!

END MODULE iopars
