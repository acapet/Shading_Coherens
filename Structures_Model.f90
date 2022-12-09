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
!
!
!*****************************************************************
!
! *Structures_Model* Routines for the structures (dry cells, thin dams, weirs,
!                    barriers) and discharge model units
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90  V2.12.1
!
! $Date: 2021-11-28 10:15:48 +0100 (So, 28 Nov 2021) $
!
! $Revision: 1435 $
!
! Description -
!
! Reference -
!
! Routines - allocate_struc_arrays, define_dischr_data, discharges, dry_cells,
!            local_struc_arrays, momentum_discharge_2d, momentum_discharge_3d,
!            read_dischr_data, read_dry_cells, read_thin_dams, read_weirs,
!            scalar_discharge, surface_discharge, thin_dams,
!            update_dischr_data, weirs_bstress, weirs_depth, weirs_loss,
!            weirs_mask, weirs_sink, write_dischr_data, write_dischr_spec,
!            write_dry_cells, write_thin_dams, write_weirs, read_mpv, write_mpv
!            
!
!*****************************************************************
!

!=================================================================

SUBROUTINE allocate_struc_arrays
!*****************************************************************
!
! *allocate_struc_arrays* Allocate structure/discharge arrays
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.11
!
! Description - allocation of the structures arrays
!
! Reference -
!
! Calling program - initialise_model
!
! Module calls - error_alloc, error_alloc_struc, num_halo
!
!*****************************************************************
!
USE iopars
USE paralpars
USE structures
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE time_routines, ONLY: log_timer_in, log_timer_out

!AC 16.11.2022
USE grid
USE gridpars

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: npcc


procname(pglev+1) = 'allocate_struc_arrays'
CALL log_timer_in(npcc)

!
!1. Dry cells
!------------
!

IF (iopt_drycel.EQ.1) THEN

   ALLOCATE (idry(numdry),STAT=errstat)
   CALL error_alloc('idry',1,(/numdry/),kndint)
   IF (numdry.GT.0) idry = 0
   ALLOCATE (jdry(numdry),STAT=errstat)
   CALL error_alloc('jdry',1,(/numdry/),kndint)
   IF (numdry.GT.0) jdry = 0

ENDIF

!
!2. Thin dams
!------------
!

IF (iopt_thndam.EQ.1) THEN

   ALLOCATE (ithinu(numthinu),STAT=errstat)
   CALL error_alloc('ithinu',1,(/numthinu/),kndint)
   IF (numthinu.GT.0) ithinu = 0
   ALLOCATE (jthinu(numthinu),STAT=errstat)
   CALL error_alloc('jthinu',1,(/numthinu/),kndint)
   IF (numthinu.GT.0) jthinu = 0
   ALLOCATE (ithinv(numthinv),STAT=errstat)
   CALL error_alloc('ithinv',1,(/numthinv/),kndint)
   IF (numthinv.GT.0) ithinv = 0
   ALLOCATE (jthinv(numthinv),STAT=errstat)
   CALL error_alloc('jthinv',1,(/numthinv/),kndint)
   IF (numthinv.GT.0) jthinv = 0

ENDIF

!
!3. Weirs/barriers
!-----------------
!

IF (iopt_weibar.EQ.1) THEN

   ALLOCATE (iwbaru(numwbaru),STAT=errstat)
   CALL error_alloc('iwbaru',1,(/numwbaru/),kndint)
   IF (numwbaru.GT.0) iwbaru = 0

   ALLOCATE (iwbarv(numwbarv),STAT=errstat)
   CALL error_alloc('iwbarv',1,(/numwbarv/),kndint)
   IF (numwbarv.GT.0) iwbarv = 0

   ALLOCATE (jwbaru(numwbaru),STAT=errstat)
   CALL error_alloc('jwbaru',1,(/numwbaru/),kndint)
   IF (numwbaru.GT.0) iwbarv = 0

   ALLOCATE (jwbarv(numwbarv),STAT=errstat)
   CALL error_alloc('jwbarv',1,(/numwbarv/),kndint)
   IF (numwbarv.GT.0) jwbarv = 0

   ALLOCATE (oricoefu(numwbaru),STAT=errstat)
   CALL error_alloc('oricoefu',1,(/numwbaru/),kndrtype)
   IF (numwbaru.GT.0) oricoefu = 0.0

   ALLOCATE (oricoefv(numwbarv),STAT=errstat)
   CALL error_alloc('oricoefv',1,(/numwbarv/),kndrtype)
   IF (numwbarv.GT.0) oricoefv = 0.0

   ALLOCATE (oriwidthu(numwbaru),STAT=errstat)
   CALL error_alloc('oriwidthu',1,(/numwbaru/),kndrtype)
   IF (numwbaru.GT.0) oriwidthu = 0.0

   ALLOCATE (oriwidthv(numwbarv),STAT=errstat)
   CALL error_alloc('oriwidthv',1,(/numwbarv/),kndrtype)
   IF (numwbarv.GT.0) oriwidthv = 0.0

   ALLOCATE (oriheightu(numwbaru),STAT=errstat)
   CALL error_alloc('oriheightu',1,(/numwbaru/),kndrtype)
   IF (numwbaru.GT.0) oriheightu = 0.0

   ALLOCATE (oriheightv(numwbarv),STAT=errstat)
   CALL error_alloc('oriheightv',1,(/numwbarv/),kndrtype)
   IF (numwbarv.GT.0) oriheightv = 0.0

   ALLOCATE (wbarcoefu(numwbaru),STAT=errstat)
   CALL error_alloc('wbarcoefu',1,(/numwbaru/),kndrtype)
   IF (numwbaru.GT.0) wbarcoefu = 0.0

   ALLOCATE (wbarcoefv(numwbarv),STAT=errstat)
   CALL error_alloc('wbarcoefv',1,(/numwbarv/),kndrtype)
   IF (numwbarv.GT.0) wbarcoefv = 0.0

   ALLOCATE (wbarcrestu(numwbaru),STAT=errstat)
   CALL error_alloc('wbarcrestu',1,(/numwbaru/),kndrtype)
   IF (numwbaru.GT.0) wbarcrestu = 0.0

   ALLOCATE (wbarcrestv(numwbarv),STAT=errstat)
   CALL error_alloc('wbarcrestv',1,(/numwbarv/),kndrtype)
   IF (numwbarv.GT.0) wbarcrestv = 0.0

   ALLOCATE (wbarmodlu(numwbaru),STAT=errstat)
   CALL error_alloc('wbarmodlu',1,(/numwbaru/),kndrtype)
   IF (numwbaru.GT.0) wbarmodlu = 0.0

   ALLOCATE (wbarmodlv(numwbarv),STAT=errstat)
   CALL error_alloc('wbarmodlv',1,(/numwbarv/),kndrtype)
   IF (numwbarv.GT.0) wbarmodlv = 0.0

   ALLOCATE (wbarelossu(numwbaru),STAT=errstat)
   CALL error_alloc('wbarelossu',1,(/numwbaru/),kndrtype)
   IF (numwbaru.GT.0) wbarelossu = 0.0

   ALLOCATE (wbarelossv(numwbarv),STAT=errstat)
   CALL error_alloc('wbarelossv',1,(/numwbarv/),kndrtype)
   IF (numwbarv.GT.0) wbarelossv = 0.0

   ALLOCATE (nowbaruprocs(nprocs),STAT=errstat)
   CALL error_alloc('nowbaruprocs',1,(/nprocs/),kndint)
   nowbaruprocs = 0

   ALLOCATE (nowbarvprocs(nprocs),STAT=errstat)
   CALL error_alloc('nowbarvprocs',1,(/nprocs/),kndint)
   nowbarvprocs = 0

   ALLOCATE (indexwbaruprocs(numwbaru,nprocs),STAT=errstat)
   CALL error_alloc('indexwbaruprocs',2,(/numwbaru,nprocs/),kndint)
   indexwbaruprocs = 0

   ALLOCATE (indexwbarvprocs(numwbarv,nprocs),STAT=errstat)
   CALL error_alloc('indexwbarvprocs',2,(/numwbarv,nprocs/),kndint)
   indexwbarvprocs = 0

ENDIF

!
!4. Discharges
!-------------
!

IF (iopt_dischr.EQ.1) THEN

   ALLOCATE (idis(numdis),STAT=errstat)
   CALL error_alloc('idis',1,(/numdis/),kndint)
   idis = 0

   ALLOCATE (jdis(numdis),STAT=errstat)
   CALL error_alloc('jdis',1,(/numdis/),kndint)
   jdis = 0

   ALLOCATE (kdis(numdis),STAT=errstat)
   CALL error_alloc('kdis',1,(/numdis/),kndint)
   kdis = 0

   ALLOCATE (kdistype(numdis),STAT=errstat)
   CALL error_alloc('kdistype',1,(/numdis/),kndint)
   kdistype = 0

   ALLOCATE(indexdisloc(numdis),STAT=errstat)
   CALL error_alloc('indexdisloc',1,(/numdis/),kndint)
   indexdisloc = 0

   ALLOCATE (disarea(numdis),STAT=errstat)
   CALL error_alloc('disarea',1,(/numdis/),kndrtype)
   disarea = 0.0

   ALLOCATE (disdir(numdis),STAT=errstat)
   CALL error_alloc('disdir',1,(/numdis/),kndrtype)
   disdir = 0.0

   ALLOCATE (disflag(numdis),STAT=errstat)
   CALL error_alloc('disflag',1,(/numdis/),kndlog)
   disflag = .FALSE.

   ALLOCATE (disspeed(numdis),STAT=errstat)
   CALL error_alloc('disspeed',1,(/numdis/),kndrtype)
   disspeed = 0.0

   ALLOCATE (disvol(numdis),STAT=errstat)
   CALL error_alloc('disvol',1,(/numdis/),kndrtype)
   disvol = 0.0

   ALLOCATE (xdiscoord(numdis),STAT=errstat)
   CALL error_alloc('xdiscoord',1,(/numdis/),kndrtype)
   xdiscoord = 0.0

   ALLOCATE (ydiscoord(numdis),STAT=errstat)
   CALL error_alloc('ydiscoord',1,(/numdis/),kndrtype)
   ydiscoord = 0.0

   ALLOCATE (zdiscoord(numdis),STAT=errstat)
   CALL error_alloc('zdiscoord',1,(/numdis/),kndrtype)
   zdiscoord = 0.0

   ALLOCATE (nodisprocs(nprocs),STAT=errstat)
   CALL error_alloc('nodisprocs',1,(/nprocs/),kndint)
   nodisprocs = 0

   ALLOCATE (indexdisprocs(numdis,nprocs),STAT=errstat)
   CALL error_alloc('indexdisprocs',2,(/numdis,nprocs/),kndint)
   indexdisprocs = 0

ENDIF

!
!5. MPV
!-------------
!
IF (iopt_mpvcov.EQ.1) THEN
   ALLOCATE (mpvcov(ncloc,nrloc),STAT=errstat)
   CALL error_alloc('mpvcov',2,(/ncloc,nrloc/),kndrtype)
   mpvcov = 0.0
ENDIF

CALL log_timer_out(npcc,itm_structs)

RETURN

END SUBROUTINE allocate_struc_arrays

!========================================================================

SUBROUTINE define_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)
!************************************************************************
!
! *define_dischr_data* Obtain discharge data
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.11
!
! Description -
!
! Reference -
!
! Calling program - update_dischr_data
!
! External calls - read_dischr_data, usrdef_dischr_data, write_dischr_data
!
! Module calls - error_file, suspend_procs
!
!************************************************************************
!
USE iopars
USE syspars
USE error_routines, ONLY: error_file
USE time_routines, ONLY: log_timer_in, log_timer_out, suspend_proc

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: disdata

!
! Name         Type    Purpose
!----------------------------------------------------------------------------
!*iddesc*      INTEGER Data file id
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!*nodat*       INTEGER Number of discharge locations
!*novars*      INTEGER Number of input variables
!
!----------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: define_data
INTEGER :: endfile, nosecs


procname(pglev+1) = 'define_dischr_data'
CALL log_timer_in()

!
!1. Initialise
!-------------
!

endfile = modfiles(iddesc,ifil,1)%endfile
nosecs = 0
disdata = 0.0

!
!2. Open file on first call
!--------------------------
!

DO WHILE (modfiles(iddesc,ifil,1)%iostat.EQ.0)
!  ---standard
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)

!  ---user-dfined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)
   ENDIF

!   ---suspend/abort if needed
   IF (modfiles(iddesc,ifil,1)%iostat.EQ.-1) THEN
      SELECT CASE (endfile)
         CASE (0:1)
            CALL error_file(ierrno_fopen,filepars=modfiles(iddesc,ifil,1))
         CASE (2)
            CALL suspend_proc(nowaitsecs)
            nosecs = nosecs + nowaitsecs
            IF (nosecs.GT.maxwaitsecs) THEN
               CALL error_file(ierrno_fopen,filepars=modfiles(iddesc,ifil,1))
            ENDIF
      END SELECT
   ENDIF
         
ENDDO

!
!3. Obtain data
!-------------
!

define_data = .TRUE.

DO WHILE (define_data)

!  ---read
   IF (modfiles(iddesc,ifil,1)%status.EQ.'R') THEN
      CALL read_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)
!   ---user-defined
   ELSEIF (modfiles(iddesc,ifil,1)%status.EQ.'N') THEN
      CALL usrdef_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)
   ENDIF

!  ---Suspend/abort if needed
   IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
      SELECT CASE (endfile)
         CASE (0)
            CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
         CASE (1)
            define_data = .FALSE.
         CASE (2)
            CALL suspend_proc(nowaitsecs)
            nosecs = nosecs + nowaitsecs
            IF (nosecs.GT.maxwaitsecs) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
            ENDIF
      END SELECT
   ELSE
      define_data = .FALSE.
   ENDIF

ENDDO
 
!
!4. Write data
!-------------
!

IF (modfiles(iddesc,ifil,2)%defined.AND.&
  & modfiles(iddesc,ifil,1)%iostat.LT.3) THEN
   CALL write_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)
ENDIF

CALL log_timer_out()


RETURN

END SUBROUTINE define_dischr_data

!=================================================================

SUBROUTINE discharges
!*****************************************************************
!
! *discharges* define discharge locations, rates and directions
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference - 
!
! Calling program - coherens_main, initialise_model
!
! External calls - update_dischr_data
!
! Module calls - combine_stats_glb, data_to_model_hcoords, error_alloc,
!                local_proc
!
!*****************************************************************
!
USE datatypes
USE depths
USE grid
USE gridpars
USE iopars
USE modids
USE paralpars
USE structures
USE switches
USE syspars
USE timepars
USE datatypes_init, ONLY: hrelativecoords_init
USE error_routines, ONLY: error_alloc
USE grid_interp, ONLY: data_to_model_hcoords
USE grid_routines, ONLY: local_proc
USE paral_comms, ONLY: combine_stats_glb
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
LOGICAL :: update
INTEGER :: i, iproc, j, k, lloc, n, nloc, npcc
INTEGER (KIND=kndilong), SAVE, DIMENSION(3,2) :: nosecsdat
REAL :: gsig
REAL, DIMENSION(numdis,3) :: disvals
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: disdata
INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ind, jnd 
TYPE (HRelativeCoords), DIMENSION(numdis,1) :: hdiscoords


procname(pglev+1) = 'discharges'
CALL log_timer_in(npcc)

!
!1. Allocate on first call
!-------------------------
!

IF (nt.EQ.0) THEN
   ALLOCATE (disdata(numdis,6,2),STAT=errstat)
   CALL error_alloc('disdata',3,(/numdis,6,2/),kndrtype)
ENDIF

!
!2. Update horizontal locations
!------------------------------
!
!2.1 Update geographical locations
!---------------------------------
!

CALL update_dischr_data(io_disloc,1,disdata(:,1:3,:),disvals(:,1:3),numdis,3,&
                      & nosecsdat(1,:),update)
IF (.NOT.update) GOTO 1000
xdiscoord = disvals(:,1); ydiscoord = disvals(:,2); zdiscoord = disvals(:,3)

! 
!2.2 Discharge locations within the model grid
!---------------------------------------------
!

CALL hrelativecoords_init(hdiscoords,.TRUE.)
CALL data_to_model_hcoords(hdiscoords,'UV ',numdis,1,xdiscoord,ydiscoord,&
                         & land=iopt_dischr_land.EQ.2)
idis = hdiscoords(:,1)%icoord; jdis = hdiscoords(:,1)%jcoord

!
!2.3 Apply masks
!----------------
!


n_230: DO n=1,numdis
   i = idis(n); j = jdis(n)
   disflag(n) = MERGE(.FALSE.,.TRUE.,idis(n).EQ.int_fill)
ENDDO n_230

!
!2.4 Local coordinates
!---------------------
!
!2.4.1 Allocate
!--------------
!

ALLOCATE (ind(numdis),STAT=errstat)
CALL error_alloc('ind',1,(/numdis/),kndint)
ALLOCATE (jnd(numdis),STAT=errstat)
CALL error_alloc('jnd',1,(/numdis/),kndint)

!
!2.4.2 Evaluate local arrays
!---------------------------
!

iproc_242: DO iproc=1,nprocs

   lloc = 0; nodisprocs(iproc) = 0

!  ---interior points
   n_2421: DO n=1,numdis
      i = idis(n); j = jdis(n)
      IF (local_proc(i,j,iproc=iproc)) THEN
         lloc = lloc + 1
         ind(lloc) = i-nc1procs(iproc)+1
         jnd(lloc) = j-nr1procs(iproc)+1
         indexdisprocs(lloc,iproc) = n
         nodisprocs(iproc) = nodisprocs(iproc) + 1
      ENDIF
   ENDDO n_2421
   IF (idloc.EQ.idprocs(iproc)) numdisloc = lloc

!  ---locations in the western/southern halo
   n_2422: DO n=1,numdis
      i = idis(n); j = jdis(n)
      IF (i.EQ.(nc1procs(iproc)-1).AND.&
        & j.GE.nr1procs(iproc).AND.j.LE.nr2procs(iproc)) THEN
         lloc = lloc + 1
         ind(lloc) = 0
         jnd(lloc) = j-nr1procs(iproc)+1
         indexdisprocs(lloc,iproc) = n
      ENDIF
      IF (j.EQ.(nr1procs(iproc)-1).AND.&
        & i.GE.nc1procs(iproc).AND.i.LE.nc2procs(iproc)) THEN
         lloc = lloc + 1
         ind(lloc) = i-nc1procs(iproc)+1
         jnd(lloc) = 0
         indexdisprocs(lloc,iproc) = n
      ENDIF
   ENDDO n_2422
   IF (idloc.EQ.idprocs(iproc)) THEN
      numdisloc_ext = lloc   
      IF (ALLOCATED(idisloc)) DEALLOCATE (idisloc,jdisloc)
      ALLOCATE(idisloc(numdisloc_ext),STAT=errstat)
      CALL error_alloc('idisloc',1,(/numdisloc_ext/),kndint)
      ALLOCATE(jdisloc(numdisloc_ext),STAT=errstat)
      CALL error_alloc('jdisloc',1,(/numdisloc_ext/),kndint)
      IF (numdisloc_ext.GT.0) THEN
         idisloc = ind(1:numdisloc_ext)
         jdisloc = jnd(1:numdisloc_ext)
         indexdisloc(1:numdisloc_ext) = indexdisprocs(1:numdisloc_ext,iproc)
      ENDIF
   ENDIF

ENDDO iproc_242

!
!2.5 Deallocate
!--------------
!

DEALLOCATE (ind,jnd)

1000 CONTINUE

!
!3. Vertical locations
!---------------------
!

IF (predstep.AND.iopt_grid_nodim.EQ.3) THEN

   nloc_310: DO nloc=1,numdisloc_ext
      i = idisloc(nloc); j = jdisloc(nloc)
      n = indexdisloc(nloc)
      IF (disflag(n)) THEN
         SELECT CASE(kdistype(n))
            CASE (0); kdis(n) = 0
            CASE (1); kdis(n) = NINT(zdiscoord(n))
            CASE (2)
               gsig = zdiscoord(n)/deptotatc(i,j)
               IF (gsig.LE.1.0) THEN
                  k = 1
                  DO WHILE (k.LE.nz.AND.gscoordatw(i,j,k+1).LE.gsig)
                     k = k + 1
                  END DO
                  kdis(n) = k
               ELSE
                  kdis(n) = -1
               ENDIF
            CASE (3)
               gsig = (1.0-zdiscoord(n)/deptotatc(i,j))
               IF (gsig.GE.0.0) THEN
                  k = 1
                  DO WHILE (k.LE.nz.AND.gscoordatw(i,j,k+1).LE.gsig)
                     k = k + 1
                  END DO
                  kdis(n) = k
               ELSE
                  kdis(n) = -1
               ENDIF
         END SELECT
      ENDIF
      IF (kdis(n).EQ.-1) disflag(n) = .FALSE.
   ENDDO nloc_310

   CALL combine_stats_glb(kdis,numdis,nodisprocs,indexdisprocs,&
                        & iarr_kdis,commall=.TRUE.)
   CALL combine_stats_glb(disflag,numdis,nodisprocs,indexdisprocs,&
                        & iarr_disflag,commall=.TRUE.)

ENDIF

!
!4. Volume discharge data
!------------------------
!

CALL update_dischr_data(io_disvol,1,disdata(:,4:4,:),disvol,numdis,1,&
                      & nosecsdat(2,:),update)

IF (update) THEN
   n_410: DO nloc=1,numdisloc_ext
      n = indexdisloc(nloc)
      IF (disflag(n)) THEN
         i = idisloc(nloc); j = jdisloc(nloc)
         disspeed(n) = disvol(n)/garea(i,j)
      ENDIF
   ENDDO n_410
   CALL combine_stats_glb(disspeed,numdis,nodisprocs,indexdisprocs,&
                        & iarr_disspeed,commall=.TRUE.)
ENDIF

!
!5. Momentum discharge data
!--------------------------
!

CALL update_dischr_data(io_discur,1,disdata(:,5:6,:),disvals(:,1:2),numdis,2,&
                      & nosecsdat(3,:),update)

IF (update) THEN
   disarea = disvals(:,1)
   disdir = disvals(:,2)

!  ---correct direction for rotated or curvilinear grid   
   IF (iopt_grid_htype.EQ.3.OR.surfacegrids(igrd_model,1)%rotated) THEN
      nloc_510: DO nloc=1,numdisloc
         i = idisloc(nloc); j = jdisloc(nloc)
         n = indexdisloc(nloc)
         disdir(n) = disdir(n) - gangleatc(i,j)
      ENDDO nloc_510
      CALL combine_stats_glb(disdir,numdis,nodisprocs,indexdisprocs,&
                           & iarr_disdir,commall=.TRUE.)
   ENDIF

ENDIF

!
!6. Deallocate at final step
!---------------------------
!

IF (cold_start.OR.nt.EQ.nstep) DEALLOCATE (disdata)

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE discharges

!=================================================================

SUBROUTINE dry_cells
!*****************************************************************
!
! *dry_cells* Dry cells model unit
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description - dry cells model unit 
!
! Reference - 
!
! Calling program - mask_function, pointer_arrays
!
! Module calls - error_limits_arr 
!
!*****************************************************************
!
USE depths
USE gridpars
USE iopars
USE physpars
USE structures
USE error_routines, ONLY: error_limits_arr
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: i, j, l, npcc


procname(pglev+1) = 'dry_cells'
CALL log_timer_in(npcc)

!
!1. Check dry cell locations
!---------------------------
!

l_110: DO l=1,numdry
   i = idry(l); j = jdry(l)
   CALL error_limits_arr(i,'idry',1,nc-1,1,indx=(/l/))
   CALL error_limits_arr(j,'jdry',1,nr-1,1,indx=(/l/))
ENDDO l_110

!
!2. Set bathymetry flag
!----------------------
!

l_210: DO l=1,numdry
   i = idry(l); j = jdry(l)
   depmeanglb(i,j) = depmean_flag
ENDDO l_210

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE dry_cells

!=================================================================

SUBROUTINE local_struc_arrays
!*****************************************************************
!
! *local_struc_arrays* local arrays for structures module
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Reference - 
!
! Description - defines local values for the structures arrays
!               parameters
!
! Calling program - initialise_model
!
! Module calls - error_alloc, local_proc
!
!*****************************************************************
!
USE grid
USE iopars
USE paralpars
USE structures
USE switches
USE syspars
USE error_routines, ONLY: error_alloc
USE grid_routines, ONLY: local_proc
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: i, ii, iiloc, iproc, j, jj, jjloc, nolocsu, nolocsv, nolocswb, npcc
INTEGER, ALLOCATABLE, DIMENSION(:) :: indexarr, indu, indv, jndu, jndv


procname(pglev+1) = 'local_struc_arrays'
CALL log_timer_in(npcc)

!
!1. Allocate work space arrays
!-----------------------------
!

nolocsu = MAX(numthinu,numwbaru)
nolocsv = MAX(numthinv,numwbarv)
nolocswb = MAX(numwbaru,numwbarv)

ALLOCATE (indu(nolocsu),STAT=errstat)
CALL error_alloc('indu',1,(/nolocsu/),kndint)
IF (nolocsu.GT.0) indu = 0

ALLOCATE (jndu(nolocsu),STAT=errstat)
CALL error_alloc('jndu',1,(/nolocsu/),kndint)
IF (nolocsu.GT.0) jndu = 0

ALLOCATE (indv(nolocsv),STAT=errstat)
CALL error_alloc('indv',1,(/nolocsv/),kndint)
IF (nolocsv.GT.0) indv = 0

ALLOCATE (jndv(nolocsv),STAT=errstat)
CALL error_alloc('jndv',1,(/nolocsv/),kndint)
IF (nolocsv.GT.0) jndv = 0

IF (iopt_weibar.EQ.1) THEN
   ALLOCATE (indexarr(nolocswb),STAT=errstat)
   CALL error_alloc('indexarr',1,(/nolocswb/),kndint)
   IF (nolocswb.GT.0) indexarr = 0
ENDIF

!
!2. Thin dams (local)
!--------------------
!

IF (iopt_thndam.EQ.1) THEN

!
!2.1 U-nodes
!-----------
!

   iproc_210: DO iproc=1,nprocs

!     ---evaluate
      iiloc = 0
      ii_211: DO ii=1,numthinu
         i = ithinu(ii); j = jthinu(ii)
         IF (local_proc(i,j,iproc=iproc)) THEN
            iiloc = iiloc + 1
            indu(iiloc) = i-nc1procs(iproc)+1
            jndu(iiloc) = j-nr1procs(iproc)+1
         ENDIF
      ENDDO ii_211
      IF (idloc.EQ.idprocs(iproc)) THEN
         numthinuloc = iiloc
         ALLOCATE(ithinuloc(numthinuloc),STAT=errstat)
         CALL error_alloc('ithinuloc',1,(/numthinuloc/),kndint)
         ALLOCATE(jthinuloc(numthinuloc),STAT=errstat)
         CALL error_alloc('jthinuloc',1,(/numthinuloc/),kndint)
         IF (numthinuloc.GT.0) THEN
            ithinuloc = indu(1:numthinuloc)
            jthinuloc = jndu(1:numthinuloc)
         ENDIF
      ENDIF

   ENDDO iproc_210

!
!2.2 V-nodes
!-----------
!

   iproc_220: DO iproc=1,nprocs

!     ---evaluate
      jjloc = 0  
      jj_221: DO jj=1,numthinv
         i = ithinv(jj); j = jthinv(jj)
         IF (local_proc(i,j,iproc=iproc)) THEN
            jjloc = jjloc + 1
            indv(jjloc) = i-nc1procs(iproc)+1
            jndv(jjloc) = j-nr1procs(iproc)+1
         ENDIF
      ENDDO jj_221
      IF (idloc.EQ.idprocs(iproc)) THEN
         numthinvloc = jjloc
         ALLOCATE(ithinvloc(numthinvloc),STAT=errstat)
         CALL error_alloc('ithinvloc',1,(/numthinvloc/),kndint)
         ALLOCATE(jthinvloc(numthinvloc),STAT=errstat)
         CALL error_alloc('jthinvloc',1,(/numthinvloc/),kndint)
         IF (numthinvloc.GT.0) THEN
            ithinvloc = indv(1:numthinvloc)
            jthinvloc = jndv(1:numthinvloc)
         ENDIF
      ENDIF
      
   ENDDO iproc_220

ENDIF

!
!4. Weirs/barriers (local)
!-------------------------
!
!4.1 U-nodes
!-----------
!

IF (iopt_weibar.EQ.1) THEN

   iproc_410: DO iproc=1,nprocs

!     ---evaluate
      iiloc = 0; nowbaruprocs(iproc) = 0
      ii_411: DO ii=1,numwbaru
         i = iwbaru(ii); j = jwbaru(ii)
         IF (local_proc(i,j,iproc=iproc)) THEN
            iiloc = iiloc + 1
            indu(iiloc) = i-nc1procs(iproc)+1
            jndu(iiloc) = j-nr1procs(iproc)+1
            indexwbaruprocs(iiloc,iproc) = ii
            nowbaruprocs(iproc) = nowbaruprocs(iproc) + 1
         ENDIF
      ENDDO ii_411
      IF (idloc.EQ.idprocs(iproc)) THEN
         numwbaruloc = iiloc
         ALLOCATE(iwbaruloc(numwbaruloc),STAT=errstat)
         CALL error_alloc('iwbaruloc',1,(/numwbaruloc/),kndint)
         ALLOCATE(jwbaruloc(numwbaruloc),STAT=errstat)
         CALL error_alloc('jwbaruloc',1,(/numwbaruloc/),kndint)
         ALLOCATE (indexwbaru(numwbaruloc),STAT=errstat)
         CALL error_alloc('indexwbaru',1,(/numwbaruloc/),kndreal)
         IF (numwbaruloc.GT.0) THEN
            iwbaruloc = indu(1:numwbaruloc)
            jwbaruloc = jndu(1:numwbaruloc)
            indexwbaru = indexwbaruprocs(1:numwbaruloc,iproc)
         ENDIF
      ENDIF

   ENDDO iproc_410

!
!4.2 V-nodes
!-----------
!

   iproc_420: DO iproc=1,nprocs

!     ---evaluate
      jjloc = 0; nowbarvprocs(iproc) = 0
      jj_421: DO jj=1,numwbarv
         i = iwbarv(jj); j = jwbarv(jj)
         IF (local_proc(i,j,iproc=iproc)) THEN
            jjloc = jjloc + 1
            indv(jjloc) = i-nc1procs(iproc)+1
            jndv(jjloc) = j-nr1procs(iproc)+1
            indexwbarvprocs(jjloc,iproc) = ii
            nowbarvprocs(iproc) = nowbarvprocs(iproc) + 1
         ENDIF
      ENDDO jj_421
      IF (idloc.EQ.idprocs(iproc)) THEN
         numwbarvloc = jjloc
         ALLOCATE(iwbarvloc(numwbarvloc),STAT=errstat)
         CALL error_alloc('iwbarvloc',1,(/numwbarvloc/),kndint)
         ALLOCATE(jwbarvloc(numwbarvloc),STAT=errstat)
         CALL error_alloc('jwbarvloc',1,(/numwbarvloc/),kndint)
         ALLOCATE (indexwbarv(numwbarvloc),STAT=errstat)
         CALL error_alloc('indexwbarv',1,(/numwbarvloc/),kndrtype)
         IF (numwbarvloc.GT.0) THEN
            iwbarvloc = indv(1:numwbarvloc)
            jwbarvloc = jndv(1:numwbarvloc)
            indexwbarv = indexwbarvprocs(1:numwbarvloc,iproc)
         ENDIF
      ENDIF

   ENDDO iproc_420

ENDIF

!
!5. Deallocate work space arrays
!-------------------------------
!

DEALLOCATE (indu,indv,jndu,jndv)

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE local_struc_arrays

!=================================================================

SUBROUTINE momentum_discharge_2d(sourceatu,sourceatv)
!*****************************************************************
!
! *momentum_discharge_2d* source term for 2-D momentum discharge
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling program - current_2d
!
! Module calls -
!
!*****************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc) :: sourceatu, sourceatv

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*sourceatu* REAL    Source terms in the U- or u-momentum equation     [m^2/s^2]
!*sourceatv* REAL    Source terms in the V- or v-momentum equation     [m^2/s^2]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, j, jj, k, l, n, nloc, npcc
REAL :: udis, vdis
REAL, DIMENSION(2) :: weights


procname(pglev+1) = 'momentum_discharge_2d'
CALL log_timer_in(npcc)

!
!1. U-nodes
!----------
!

nloc_110: DO nloc=1,numdisloc_ext
   n = indexdisloc(nloc); k = kdis(n)
   IF (disflag(n).AND.disarea(n).GT.0.0) THEN
      i = idisloc(nloc); j = jdisloc(nloc)
      IF (node2du(i,j).EQ.2) THEN
         IF (node2du(i+1,j).EQ.2) THEN
            weights = 0.5
         ELSE
            weights = (/1.0,0.0/)
         ENDIF
      ELSE
         IF (node2du(i+1,j).EQ.2) THEN
            weights = (/0.0,1.0/)
         ELSE
            weights = 0.0
         ENDIF
      ENDIF
      IF (SUM(weights).GT.0) THEN
         ii_111: DO ii=i,i+1
            IF (ii.GE.1.AND.ii.LE.ncloc) THEN
               l = MERGE(1,2,ii.EQ.i)
               udis = disvol(n)*COS(disdir(n))/disarea(n)
               IF (k.EQ.0) THEN
                  sourceatu(ii,j) = sourceatu(ii,j) + weights(l)*disspeed(n)*&
                                  & (udis-umvel(i,j))
               ELSE
                  sourceatu(ii,j) = sourceatu(ii,j) + weights(l)*disspeed(n)*&
                                  & (udis-uvel(i,j,k))
               ENDIF
            ENDIF
         ENDDO ii_111
      ENDIF
   ENDIF
ENDDO nloc_110

!
!2. V-nodes
!----------
!

nloc_210: DO nloc=1,numdisloc_ext
   n = indexdisloc(nloc); k = kdis(n)
   IF (disflag(n).AND.disarea(n).GT.0.0) THEN
      i = idisloc(nloc); j = jdisloc(nloc)
      IF (node2dv(i,j).EQ.2) THEN
         IF (node2dv(i,j+1).EQ.2) THEN
            weights = 0.5
         ELSE
            weights = (/1.0,0.0/)
         ENDIF
      ELSE
         IF (node2dv(i,j+1).EQ.2) THEN
            weights = (/0.0,1.0/)
         ELSE
            weights = 0.0
         ENDIF
      ENDIF
      IF (SUM(weights).GT.0.0) THEN
         jj_211: DO jj=j,j+1
            IF (jj.GE.1.AND.jj.LE.nrloc) THEN
               l = MERGE(1,2,jj.EQ.j)
               vdis = disvol(n)*SIN(disdir(n))/disarea(n)
               IF (k.EQ.0) THEN
                  sourceatv(i,jj) = sourceatv(i,jj) + weights(l)*disspeed(n)*&
                               & (vdis-vmvel(i,j))
               ELSE
                  sourceatv(i,jj) = sourceatv(i,jj) + weights(l)*disspeed(n)*&
                               & (vdis-vvel(i,j,k))
               ENDIF
            ENDIF
         ENDDO jj_211
      ENDIF
   ENDIF
ENDDO nloc_210

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE momentum_discharge_2d

!=================================================================

SUBROUTINE momentum_discharge_3d(sourceatu,sourceatv)
!*****************************************************************
!
! *momentum_discharge_3d* source term for 3-D momentum discharge
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling program - current_pred
!
! Module calls -
!
!*****************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz) :: sourceatu, sourceatv

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*sourceatu* REAL    Source terms in the U- or u-momentum equation       [m/s^2]
!*sourceatv* REAL    Source terms in the V- or v-momentum equation       [m/s^2]
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ii, j, jj, k, l, n, nloc, npcc
REAL :: udis, vdis
REAL, DIMENSION(2) :: weights


procname(pglev+1) = 'momentum_discharge_3d'
CALL log_timer_in(npcc)

!
!1. U-nodes
!----------
!

nloc_110: DO nloc=1,numdisloc_ext
   n = indexdisloc(nloc); k = kdis(n)
   IF (disflag(n).AND.disarea(n).GT.0.0) THEN
      i = idisloc(nloc); j = jdisloc(nloc)
      IF (node2du(i,j).EQ.2) THEN
         IF (node2du(i+1,j).EQ.2) THEN
            weights = 0.5
         ELSE
            weights = (/1.0,0.0/)
         ENDIF
      ELSE
         IF (node2du(i+1,j).EQ.2) THEN
            weights = (/0.0,1.0/)
         ELSE
            weights = 0.0
         ENDIF
      ENDIF
      IF (SUM(weights).GT.0) THEN
         ii_111: DO ii=i,i+1
            IF (ii.GE.1.AND.ii.LE.ncloc) THEN
               l = MERGE(1,2,ii.EQ.i)
               udis = disvol(n)*COS(disdir(n))/disarea(n)
               IF (k.EQ.0) THEN
                  sourceatu(ii,j,:) = sourceatu(ii,j,:) + &
                                    & weights(l)*disspeed(n)*&
                                    & (udis-uvel(ii,j,:))/deptotatc(i,j)
               ELSE
                  sourceatu(ii,j,k) = sourceatu(ii,j,k) + &
                                    & weights(l)*disspeed(n)*&
                                    & (udis-uvel(ii,j,k))/delzatc(i,j,k)
               ENDIF
            ENDIF
         ENDDO ii_111
      ENDIF
   ENDIF
ENDDO nloc_110

!
!2. V-nodes
!----------
!

nloc_210: DO nloc=1,numdisloc
   n = indexdisloc(nloc); k = kdis(n)
   IF (disflag(n).AND.disarea(n).GT.0.0) THEN
      i = idisloc(nloc); j = jdisloc(nloc)
      IF (node2dv(i,j).EQ.2) THEN
         IF (node2dv(i,j+1).EQ.2) THEN
            weights = 0.5
         ELSE
            weights = (/1.0,0.0/)
         ENDIF
      ELSE
         IF (node2dv(i,j+1).EQ.2) THEN
            weights = (/0.0,1.0/)
         ELSE
            weights = 0.0
         ENDIF
      ENDIF
      IF (SUM(weights).GT.0.0) THEN
         jj_211: DO jj=j,j+1
            IF (jj.GE.1.AND.jj.LE.nrloc) THEN
               l = MERGE(1,2,jj.EQ.j)
               vdis = disvol(n)*SIN(disdir(n))/disarea(n)
               IF (k.EQ.0) THEN
                  sourceatv(i,jj,:) = sourceatv(i,jj,:) + &
                                    & weights(l)*disspeed(n)*&
                                    & (vdis-vvel(i,jj,:))/deptotatc(i,j)
               ELSE
                  sourceatv(i,jj,k) = sourceatv(i,jj,k) + &
                                    & weights(l)*disspeed(n)*&
                                   & (vdis-vvel(i,jj,k))/delzatc(i,j,k)
               ENDIF
            ENDIF
         ENDDO jj_211
      ENDIF
   ENDIF
ENDDO nloc_210

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE momentum_discharge_3d

!=================================================================

SUBROUTINE read_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novarsdat)
!*****************************************************************
!
! *read_dischr_data* Read discharge data in standard format
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12.1
!
! Description -
!
! Reference -
!
! Calling program - define_dischr_data
!
! Module calls - check_value_varatts, close_filepars, error_abort,
!                error_alloc_struc, error_value_arr_struc,
!                open_filepars, read_glbatts_mod, read_time,
!                read_varatts_mod, read_vars, set_modfiles_atts,
!                set_modvars_atts, varatts_init
!
!*****************************************************************
!
USE datatypes
USE iopars
USE syspars
USE check_model, ONLY: check_value_varatts
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_abort, error_alloc_struc, &
                        & error_value_arr_struc
USE inout_routines, ONLY: open_filepars, read_glbatts_mod, read_time, &
                        & read_varatts_mod, read_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(OUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novarsdat
REAL, INTENT(INOUT), DIMENSION(nodat,novarsdat) :: disdata

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!*nodat*       INTEGER Number of discharge locations
!*novarsdat*   INTEGER Number of input variables
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: header
INTEGER, PARAMETER :: nodesc = 5
INTEGER, SAVE :: maxvars = 4
INTEGER :: idesc, ivar, nocoords1, nocoords2, novars, novars1, novars2, &
         & numvars1, numvars2
INTEGER, SAVE, DIMENSION(nodesc) :: numvars
TYPE (FileParams) :: filepars, filepars1, filepars2
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:,:) :: varatts
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts1, varatts2


procname(pglev+1) = 'read_dischr_data'
CALL log_timer_in()

filepars = modfiles(iddesc,ifil,1)

SELECT CASE (iddesc)
   CASE (io_disloc); idesc = 1
   CASE (io_disvol); idesc = 2
   CASE (io_discur); idesc = 3
   CASE (io_dissal); idesc = 4
   CASE (io_distmp); idesc = 5
END SELECT

!
!1. Initialise
!-------------
!

IF (filepars%iostat.EQ.0) THEN
   header = modfiles(iddesc,ifil,1)%header

!   
!1.1 Data attributes
!-------------------
!

   IF (header) THEN

!     ---open data file
      filepars1 = modfiles(iddesc,ifil,1)
      CALL open_filepars(filepars1)
   
!     ---file attributes
      CALL read_glbatts_mod(filepars1)

!     ---variable attributes
      nocoords1 = filepars1%nocoords
      novars1 = filepars1%novars
      numvars1 = nocoords1 + novars1
      ALLOCATE (varatts1(numvars1),STAT=errstat)
      CALL error_alloc_struc('varatts1',1,(/numvars1/),'VariableAtts')
      CALL varatts_init(varatts1)
      CALL read_varatts_mod(filepars1,varatts1)

   ENDIF
   
!
!1.2 Model attributes
!--------------------
!
   
!  ---file attributes
   filepars2 = modfiles(iddesc,ifil,1)
   CALL set_modfiles_atts(iddesc,ifil,1,filepars=filepars2)

!  ---data attributes
   nocoords2 = filepars2%nocoords
   novars2 = filepars2%novars
   numvars2 = nocoords2 + novars2
   ALLOCATE (varatts2(numvars2),STAT=errstat)
   CALL error_alloc_struc('varatts2',1,(/numvars2/),'VariableAtts')
   CALL varatts_init(varatts2)
   CALL set_modvars_atts(iddesc,ifil,1,filepars2,numvars2,varatts2)

!
!1.3 Check attributes
!--------------------
!

   IF (header) THEN
      CALL error_value_arr_struc(nocoords1,'modfiles','nocoords',nocoords2,3,&
                              & (/iddesc,ifil,1/))
      CALL check_value_varatts(iddesc,ifil,varatts1,varatts2,filepars1,&
                             & filepars2)
      CALL error_abort('read_dischr_data',ierrno_read)
   ENDIF
 
!
!1.4 Variables names used for idenfication
!-----------------------------------------
!

   IF (header) THEN
      novars = novars1; numvars(idesc) = numvars1
      IF (inflagdis) THEN
         IF (ALLOCATED(varatts)) DEALLOCATE(varatts)
         ALLOCATE (varatts(maxvars,nodesc),STAT=errstat)
         CALL error_alloc_struc('varatts',2,(/maxvars,nodesc/),'VariableAtts')
         CALL varatts_init(varatts)
         inflagdis = .FALSE.
      ENDIF
      varatts(1:numvars1,idesc) = varatts1
      filepars = filepars1
      DEALLOCATE (varatts1)
   ELSE
      novars = novars2; numvars(idesc) = numvars2
      IF (inflagdis) THEN
         IF (ALLOCATED(varatts)) DEALLOCATE(varatts)
         ALLOCATE (varatts(maxvars,nodesc),STAT=errstat)
         CALL error_alloc_struc('varatts',2,(/maxvars,nodesc/),'VariableAtts')
         CALL varatts_init(varatts)
         inflagdis = .FALSE.         
      ENDIF
      varatts(1:numvars2,idesc) = varatts2
      filepars = filepars2
      CALL open_filepars(filepars)

   ENDIF

!
!1.5 Deallocate
!--------------
!

   DEALLOCATE (varatts2)
   
   GOTO 1000
   
ENDIF

!
!2. Read data
!------------
!
!2.1 Date/time
!-------------
!

CALL read_time(ciodatetime,filepars,(/varatts(1,idesc)/))

!
!2.2 Discharge data
!------------------
!

IF (filepars%iostat.LT.3) THEN
   ivar_220: DO ivar=2,numvars(idesc)
      SELECT CASE (TRIM(varatts(ivar,idesc)%f90_name))
         CASE ('xdiscoord')
            CALL read_vars(disdata(:,ivar-1),filepars,ivar,&
                        & (/varatts(ivar,idesc)/))
         CASE ('ydiscoord')
            CALL read_vars(disdata(:,ivar-1),filepars,ivar,&
                        & (/varatts(ivar,idesc)/))
         CASE ('zdiscoord')
            CALL read_vars(disdata(:,ivar-1),filepars,ivar,&
                        & (/varatts(ivar,idesc)/))
         CASE ('disvol')
            CALL read_vars(disdata(:,ivar-1),filepars,ivar,&
                        & (/varatts(ivar,idesc)/))
         CASE ('disarea')
            CALL read_vars(disdata(:,ivar-1),filepars,ivar,&
                        & (/varatts(ivar,idesc)/))
         CASE ('disdir')
            CALL read_vars(disdata(:,ivar-1),filepars,ivar,&
                        & (/varatts(ivar,idesc)/))
         CASE ('dissal')
            CALL read_vars(disdata(:,ivar-1),filepars,ivar,&
                        & (/varatts(ivar,idesc)/))
         CASE ('distmp')
            CALL read_vars(disdata(:,ivar-1),filepars,ivar,&
                        & (/varatts(ivar,idesc)/))
      END SELECT
   ENDDO ivar_220
ENDIF

1000 modfiles(iddesc,ifil,1) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE read_dischr_data

!=================================================================

SUBROUTINE read_dry_cells
!*****************************************************************
!
! *read_dry_cells* Read external file for dry cells
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12.1
!
! Description - reads external input files for dry cells 
!
! Reference -
!
! Calling process - define_global_grid
!
! Module calls - check_value_varatts, close_filepars, error_abort,
!                error_alloc_struc, error_value_arr_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars,
!                set_modfiles_atts, set_modvars_atts, varatts_init  
!
!*****************************************************************
!
USE datatypes
USE iopars
USE structures
USE check_model, ONLY: check_value_varatts
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_abort, error_alloc_struc, error_value_arr_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: header
INTEGER :: ivar, nocoords1, nocoords2, novars, novars1, novars2, npcc
TYPE (FileParams) :: filepars, filepars1, filepars2
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts, varatts1, varatts2


procname(pglev+1) = 'read_dry_cells'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!1.1 Data attributes
!-------------------
!

header = modfiles(io_drycel,1,1)%header

IF (header) THEN

!  ---open data file
   filepars1 = modfiles(io_drycel,1,1)
   CALL open_filepars(filepars1)
   
!  ---file attributes
   CALL read_glbatts_mod(filepars1)

!  ---variable attributes
   nocoords1 = filepars1%nocoords
   novars1 = filepars1%novars
   ALLOCATE (varatts1(novars1),STAT=errstat)
   CALL error_alloc_struc('varatts1',1,(/novars1/),'VariableAtts')
   CALL varatts_init(varatts1)
   CALL read_varatts_mod(filepars1,varatts1)
   
ENDIF
   
!
!1.2 Model attributes
!--------------------
!
!---file attributes
filepars2 = modfiles(io_drycel,1,1)
CALL set_modfiles_atts(io_drycel,1,1,filepars=filepars2)

!---data attributes
nocoords2 = filepars2%nocoords
novars2 = filepars2%novars
ALLOCATE (varatts2(novars2),STAT=errstat)
CALL error_alloc_struc('varatts2',1,(/novars2/),'VariableAtts')
CALL varatts_init(varatts2)
CALL set_modvars_atts(io_drycel,1,1,filepars2,novars2,varatts2)

!
!1.3 Check attributes
!--------------------
!

IF (header) THEN
   CALL error_value_arr_struc(nocoords1,'modfiles','nocoords',nocoords2,3,&
                           & (/io_drycel,1,1/))
   CALL check_value_varatts(io_drycel,1,varatts1,varatts2,filepars1,filepars2)
   CALL error_abort('read_dry_cells',ierrno_read)
ENDIF

!
!1.4 Variables names used for idenfication
!-----------------------------------------
!

IF (header) THEN
   novars = novars1
   ALLOCATE (varatts(novars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/novars/),'VariableAtts')
   CALL varatts_init(varatts)
   varatts = varatts1
   DEALLOCATE (varatts1)
   filepars = filepars1
ELSE
   novars = novars2
   ALLOCATE (varatts(novars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/novars/),'VariableAtts')
   CALL varatts_init(varatts)
   varatts = varatts2
   filepars = filepars2
   CALL open_filepars(filepars)
ENDIF
DEALLOCATE (varatts2)

!
!2. Read dry cell locations
!--------------------------
!

ivar_210: DO ivar=1,novars

   SELECT CASE (TRIM(varatts(ivar)%f90_name))
      CASE ('idry')
         CALL read_vars(idry,filepars,ivar,(/varatts(ivar)/))
      CASE ('jdry')
         CALL read_vars(jdry,filepars,ivar,(/varatts(ivar)/))
   END SELECT

ENDDO ivar_210

!
!3. Finalise
!-----------
!
!---close file 
CALL close_filepars(filepars)
modfiles(io_drycel,1,1) = filepars

!---deallocate
DEALLOCATE(varatts)

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE read_dry_cells

!=================================================================

SUBROUTINE read_thin_dams
!*****************************************************************
!
! *read_thin_dams* Read external file for thin dams
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12.1
!
! Description - read external files for thin dams 
!
! Reference -
!
! Calling process - initialise_model
!
! Module calls - check_value_varatts, close_filepars, error_abort,
!                error_alloc_struc, error_value_arr_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars,
!                set_modfiles_atts, set_modvars_atts, varatts_init  
!
!*****************************************************************
!
USE datatypes
USE iopars
USE structures
USE check_model, ONLY: check_value_varatts
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_abort, error_alloc_struc,  error_value_arr_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                         & read_varatts_mod, read_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
LOGICAL :: header
INTEGER :: ivar, nocoords1, nocoords2, novars, novars1, novars2, npcc
TYPE (FileParams) :: filepars, filepars1, filepars2
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts, varatts1, varatts2


procname(pglev+1) = 'read_thin_dams'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!1.1 Data attributes
!-------------------
!

header = modfiles(io_thndam,1,1)%header

IF (header) THEN

!  ---open data file
   filepars1 = modfiles(io_thndam,1,1)
   CALL open_filepars(filepars1)
   
!  ---file attributes
   CALL read_glbatts_mod(filepars1)

!  ---variable attributes
   nocoords1 = filepars1%nocoords
   novars1 = filepars1%novars
   ALLOCATE (varatts1(novars1),STAT=errstat)
   CALL error_alloc_struc('varatts1',1,(/novars1/),'VariableAtts')
   CALL varatts_init(varatts1)
   CALL read_varatts_mod(filepars1,varatts1)
   
ENDIF
   
!
!1.2 Model attributes
!--------------------
!
!---file attributes
filepars2 = modfiles(io_thndam,1,1)
CALL set_modfiles_atts(io_thndam,1,1,filepars=filepars2)

!---data attributes
nocoords2 = filepars2%nocoords
novars2 = filepars2%novars
ALLOCATE (varatts2(novars2),STAT=errstat)
CALL error_alloc_struc('varatts2',1,(/novars2/),'VariableAtts')
CALL varatts_init(varatts2)
CALL set_modvars_atts(io_thndam,1,1,filepars2,novars2,varatts2)

!
!1.3 Check attributes
!--------------------
!

IF (header) THEN
   CALL error_value_arr_struc(nocoords1,'modfiles','nocoords',nocoords2,3,&
                           & (/io_thndam,1,1/))
   CALL check_value_varatts(io_thndam,1,varatts1,varatts2,filepars1,filepars2)
   CALL error_abort('read_thin_dams',ierrno_read)
ENDIF

!
!1.4 Variables names used for idenfication
!-----------------------------------------
!

IF (header) THEN
   novars = novars1
   ALLOCATE (varatts(novars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/novars/),'VariableAtts')
   CALL varatts_init(varatts)
   varatts = varatts1
   DEALLOCATE (varatts1)
   filepars = filepars1
ELSE
   novars = novars2
   ALLOCATE (varatts(novars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/novars/),'VariableAtts')
   CALL varatts_init(varatts)
   varatts = varatts2
   filepars = filepars2
   CALL open_filepars(filepars)
ENDIF
DEALLOCATE (varatts2)

!
!2. Read thin dams locations
!---------------------------
!

ivar_210: DO ivar=1,novars

   SELECT CASE (TRIM(varatts(ivar)%f90_name))
      CASE ('ithinu')
         CALL read_vars(ithinu,filepars,ivar,(/varatts(ivar)/))
      CASE ('jthinu')
         CALL read_vars(jthinu,filepars,ivar,(/varatts(ivar)/))
      CASE ('ithinv')
         CALL read_vars(ithinv,filepars,ivar,(/varatts(ivar)/))
      CASE ('jthinv')
         CALL read_vars(jthinv,filepars,ivar,(/varatts(ivar)/))
   END SELECT
ENDDO ivar_210

!
!3. Finalise
!-----------
!
!---close file
CALL close_filepars(filepars)
modfiles(io_thndam,1,1) = filepars

!---deallocate
DEALLOCATE(varatts)

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE read_thin_dams

!=================================================================

SUBROUTINE read_weirs
!*****************************************************************
!
! *read_weirs* Read external file for weirs/barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12.1
!
! Description - read external files for weirs/barriers 
!
! Reference -
!
! Calling process - initialise_model
!
! Module calls - check_value_varatts, close_filepars, error_abort,
!                error_alloc_struc, error_value_arr_struc, open_filepars,
!                read_glbatts_mod, read_varatts_mod, read_vars,
!                set_modfiles_atts, set_modvars_atts, varatts_init
!
!*****************************************************************
!
USE datatypes
USE iopars
USE structures
USE check_model, ONLY: check_value_varatts
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_abort, error_alloc_struc, &
                        & error_value_arr_struc
USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                        & read_varatts_mod, read_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
LOGICAL :: header
INTEGER :: ivar, nocoords1, nocoords2, novars, novars1, novars2, npcc
TYPE (FileParams) :: filepars, filepars1, filepars2
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts, varatts1, varatts2


procname(pglev+1) = 'read_weirs'
CALL log_timer_in(npcc)

!
!1. Initialise
!-------------
!
!1.1 Data attributes
!-------------------
!

header = modfiles(io_weibar,1,1)%header

IF (header) THEN

!  ---open data file
   filepars1 = modfiles(io_weibar,1,1)
   CALL open_filepars(filepars1)
   
!  ---file attributes
   CALL read_glbatts_mod(filepars1)

!  ---variable attributes
   nocoords1 = filepars1%nocoords
   novars1 = filepars1%novars
   ALLOCATE (varatts1(novars1),STAT=errstat)
   CALL error_alloc_struc('varatts1',1,(/novars1/),'VariableAtts')
   CALL varatts_init(varatts1)
   CALL read_varatts_mod(filepars1,varatts1)
   
ENDIF
   
!
!1.2 Model attributes
!--------------------
!
!---file attributes
filepars2 = modfiles(io_weibar,1,1)
CALL set_modfiles_atts(io_weibar,1,1,filepars2)

!---data attributes
nocoords2 = filepars2%nocoords
novars2 = filepars2%novars
ALLOCATE (varatts2(novars2),STAT=errstat)
CALL error_alloc_struc('varatts2',1,(/novars2/),'VariableAtts')
CALL varatts_init(varatts2)
CALL set_modvars_atts(io_weibar,1,1,filepars2,novars2,varatts2)

!
!1.3 Check attributes
!--------------------
!

IF (header) THEN
   CALL error_value_arr_struc(nocoords1,'modfiles','nocoords',nocoords2,3,&
                           & (/io_weibar,1,1/))
   CALL check_value_varatts(io_weibar,1,varatts1,varatts2,filepars1,filepars2)
   CALL error_abort('read_weirs',ierrno_read)

ENDIF

!
!1.4 Variables names used for idenfication
!-----------------------------------------
!

IF (header) THEN
   novars = novars1
   ALLOCATE (varatts(novars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/novars/),'VariableAtts')
   CALL varatts_init(varatts)
   varatts = varatts1
   DEALLOCATE (varatts1)
   filepars = filepars1
ELSE
   novars = novars2
   ALLOCATE (varatts(novars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/novars/),'VariableAtts')
   CALL varatts_init(varatts)
   varatts = varatts2
   filepars = filepars2
   CALL open_filepars(filepars)
ENDIF

DEALLOCATE (varatts2)

!
!2. Reading input data of weirs/barriers
!---------------------------------------
!

ivar_210: DO ivar=1,novars

   SELECT CASE (TRIM(varatts(ivar)%f90_name))

!     ---U-nodes
      CASE ('iwbaru')
         CALL read_vars(iwbaru,filepars,ivar,(/varatts(ivar)/))
      CASE ('jwbaru')
         CALL read_vars(jwbaru,filepars,ivar,(/varatts(ivar)/))
      CASE ('oricoefu')
         CALL read_vars(oricoefu,filepars,ivar,(/varatts(ivar)/))
      CASE ('oriheightu')
         CALL read_vars(oriheightu,filepars,ivar,(/varatts(ivar)/))
      CASE ('oriwidthu')
         CALL read_vars(oriwidthu,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarcoefu')
         CALL read_vars(wbarcoefu,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarcrestu')
         CALL read_vars(wbarcrestu,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarmodlu')
         CALL read_vars(wbarmodlu,filepars,ivar,(/varatts(ivar)/))

!-     ---V-nodes
      CASE ('iwbarv')
         CALL read_vars(iwbarv,filepars,ivar,(/varatts(ivar)/))
      CASE ('jwbarv')
         CALL read_vars(jwbarv,filepars,ivar,(/varatts(ivar)/))
      CASE ('oricoefv')
         CALL read_vars(oricoefv,filepars,ivar,(/varatts(ivar)/))
      CASE ('oriheightv')
         CALL read_vars(oriheightv,filepars,ivar,(/varatts(ivar)/))
      CASE ('oriwidthv')
         CALL read_vars(oriwidthv,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarcoefv')
         CALL read_vars(wbarcoefv,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarcrestv')
         CALL read_vars(wbarcrestv,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarmodlv')
         CALL read_vars(wbarmodlv,filepars,ivar,(/varatts(ivar)/))

   END SELECT

ENDDO ivar_210

!
!3. Finalise
!-----------
!
!---close file
CALL close_filepars(filepars)
modfiles(io_weibar,1,1) = filepars

!---deallocate
DEALLOCATE(varatts)

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE read_weirs

!=================================================================

SUBROUTINE scalar_discharge(disscal,source,novars,itype)
!*****************************************************************
!
! *scalar_discharge* source term for a scalar discharge
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.11
!
! Description -
!
! Reference -
!
! Calling program - salinity_equation, temperature_equation
!
! Module calls -
!
!*****************************************************************
!
USE depths
USE grid
USE gridpars
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: itype, novars
REAL, INTENT(IN), DIMENSION(numdis,novars) :: disscal
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nz,novars) :: source


!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*disscal*   REAL    Data value of discharged scalar. Unit is [psi] for wet
!                    or [psi/s] for dry discharge
!*source*    REAL    Source term in the scalar transport equation       [psi/s]
!*novars*    INTEGER Number of disxcharge variables
!*itype*     REAL    Type of discharge
!                     = 1 => wet discharge
!                     = 2 => dry discharge
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, ivar, j, k, n, nloc, npcc


procname(pglev+1) = 'scalar_discharge'
CALL log_timer_in(npcc)

!
!1. Wet discharge
!----------------
!

IF (itype.EQ.1) THEN

   ivar_110: DO ivar=1,novars
   nloc_110: DO nloc=1,numdisloc
      n = indexdisloc(nloc); k = kdis(n)
      IF (disflag(n)) THEN
         i = idisloc(nloc); j = jdisloc(nloc); k = kdis(n)
         IF (k.EQ.0) THEN
            source(i,j,:,ivar) = source(i,j,:,ivar) + &
                               & (disspeed(n)/deptotatc(i,j))*disscal(n,ivar)
         ELSE
            source(i,j,k,ivar) = source(i,j,k,ivar) + &
                               & (disspeed(n)/delzatc(i,j,k))*disscal(n,ivar)
         ENDIF
      ENDIF
   ENDDO nloc_110
   ENDDO ivar_110

!
!2. Dry discharge
!----------------
!

ELSEIF (itype.EQ.2) THEN

   ivar_210: DO ivar=1,novars
   nloc_210: DO nloc=1,numdisloc
      n = indexdisloc(nloc); k = kdis(n)
      IF (disflag(n)) THEN
         i = idisloc(nloc); j = jdisloc(nloc); k = kdis(n)
         IF (k.EQ.0) THEN
            source(i,j,:,ivar) = source(i,j,:,ivar) + &
                               & disscal(n,ivar)/(deptotatc(i,j)*garea(i,j))
         ELSE
            source(i,j,k,ivar) = source(i,j,k,ivar) + &
                               & disscal(n,ivar)/(delzatc(i,j,k)*garea(i,j))
         ENDIF
      ENDIF
   ENDDO nloc_210
   ENDDO ivar_210

ENDIF

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE scalar_discharge

!=================================================================

SUBROUTINE surface_discharge
!*****************************************************************
!
! *surface_discharge* add volume discharge to the surface elevation 
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling program - surface_elevation
!
! Module calls -
!
!*****************************************************************
!
USE depths
USE iopars
USE structures
USE timepars
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
INTEGER :: i, j, n, nloc, npcc


procname(pglev+1) = 'surface_discharge'
CALL log_timer_in(npcc)

nloc_310: DO nloc=1,numdisloc
   n = indexdisloc(nloc)
   IF (disflag(n)) THEN
      i = idisloc(nloc); j = jdisloc(nloc)
      zeta(i,j) = zeta(i,j) + delt2d*disspeed(n)
   ENDIF
ENDDO nloc_310

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE surface_discharge

!=================================================================

SUBROUTINE thin_dams(cnode)
!*****************************************************************
!
! *thin_dams* update pointer arrays at thin dams
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling process - mask_function, pointer_arrays
!
! Module calls - 
!
!*****************************************************************
!
USE grid
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=1), INTENT(IN) :: cnode

!
! Name     Type         Purpose
!--------------------------------------------------------------------
!*cnode*   CHARACTER    Velocity node ('U' or 'V')
!
!--------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: i, iiloc, j, jjloc, npcc


procname(pglev+1) = 'thin_dams'
CALL log_timer_in(npcc)

SELECT CASE (cnode)

!
!1. U-nodes
!----------
!

CASE('U')

   iiloc_110: DO iiloc=1,numthinuloc
      i = ithinuloc(iiloc)
      j = jthinuloc(iiloc)
      nodeatu(i,j,:) = 1
   ENDDO iiloc_110

!
!2. V-nodes
!----------
!

CASE('V')

   jjloc_210: DO jjloc=1,numthinvloc
      i = ithinvloc(jjloc)
      j = jthinvloc(jjloc)
      nodeatv(i,j,:) = 1
   ENDDO jjloc_210

END SELECT

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE thin_dams

!========================================================================

SUBROUTINE update_dischr_data(iddesc,ifil,disdata,disvals,nodat,novars,&
                            & nosecsdat,update)
!************************************************************************
!
! *update_dischr_data* Update discharge locations or data
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12
!
! Reference -
!
! Calling program - discharges, salinity_equation, temperature_equation
!
! External calls - define_dischr_data
!
! Module calls - diff_dates, error_file, lim_dims, loop_index
!
!*************************************************************************
!
USE iopars
USE syspars
USE timepars
USE error_routines, ONLY: error_file
USE time_routines, ONLY: diff_dates, log_timer_in, log_timer_out
USE utility_routines, ONLY: lim_dims, loop_index

IMPLICIT NONE

!
!*Arguments
!
LOGICAL, INTENT(OUT) :: update
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
INTEGER (KIND=kndilong), INTENT(INOUT), DIMENSION(2) :: nosecsdat
REAL, INTENT(INOUT), DIMENSION(nodat,novars,2) :: disdata
REAL, INTENT(INOUT), DIMENSION(nodat,novars) :: disvals

!
! Name     Type     Purpose
!------------------------------------------------------------------------------
!*disdata*  REAL    Input data at old and new time level
!*disvals*  REAL    Input data (eventually) inperpolated in time
!*iddesc*   INTEGER Data file id
!*ifil*     INTEGER No. of data file
!*nodat*    INTEGER Number of discharge locations
!*novars*   INTEGER Number of input variables
!*nosecsdat LONGINT Number of seconds since start of simulation and old/new
!                   data time
!*update*   LOGICAL Returns .TRUE. if update is made
!
!------------------------------------------------------------------------------
!
!*Local variables
!
LOGICAL :: no_time_series
CHARACTER (LEN=lentime) :: ciodatetime
INTEGER :: ivar, l, npcc
INTEGER, DIMENSION(3) :: tlims
INTEGER, DIMENSION(MaxIOFiles) :: noread
REAL (KIND=kndrlong) :: fill_value, fill_value_eps, ratio1, ratio2
REAL, DIMENSION(novars) :: datcos, datsin

IF (.NOT.modfiles(iddesc,ifil,1)%defined) THEN
   update = .FALSE.
   RETURN
ENDIF

!
!1. Return if no update
!----------------------
!

tlims = modfiles(iddesc,ifil,1)%tlims
no_time_series = lim_dims(ABS(tlims)).EQ.1
IF ((nt.GT.tlims(1).AND.no_time_series).OR.(.NOT.loop_index(tlims,nt)).OR.&
  & modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
   update = .FALSE.
   RETURN
ENDIF
fill_value = MERGE(double_fill,modfiles(iddesc,ifil,1)%fill_value,&
                 & modfiles(iddesc,ifil,1)%fill_value.EQ.0.0)
fill_value_eps = 0.00001*ABS(fill_value)

update = .TRUE.

procname(pglev+1) = 'update_dischr_data'
SELECT CASE (iddesc)
   CASE (io_disloc,io_disvol,io_discur)
      CALL log_timer_in()
   CASE DEFAULT
      CALL log_timer_in(npcc)
END SELECT

!
!2.1 Time-independent data
!-------------------------
!

IF (no_time_series) THEN
   IF (nt.EQ.tlims(1)) THEN
      ciodatetime = cdatetime_undef
      DO WHILE (ciodatetime.NE.CDateTime)
         CALL define_dischr_data(iddesc,ifil,ciodatetime,disdata(:,:,1),&
                               & nodat,novars)
         IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
            CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
         ENDIF
         CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(2))
         IF (nosecsdat(2).GT.nosecsrun) GOTO 1002
         noread(ifil) = 1
      ENDDO
   ENDIF

!
!2.2 Time series data
!--------------------
!

ELSE

!  ---first series at initial time
   IF (nt.EQ.tlims(1)) THEN
      CALL define_dischr_data(iddesc,ifil,ciodatetime,disdata(:,:,2),&
                            & nodat,novars)
      IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
         CALL error_file(ierrno_fend,filepars= modfiles(iddesc,ifil,1))
      ENDIF
      CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(2))
      IF (nosecsdat(2).GT.nosecsrun) GOTO 1002
      noread(ifil) = 1
   ENDIF

!  ---update when needed
   IF (nt.LT.nstep) THEN
      DO WHILE (nosecsdat(2).LE.nosecsrun)
!        ---store old data
         disdata(:,:,1) = disdata(:,:,2)
         nosecsdat(1) = nosecsdat(2)
!        ---read new data
         CALL define_dischr_data(iddesc,ifil,ciodatetime,disdata(:,:,2),&
                               & nodat,novars)
         noread(ifil) = noread(ifil) + 1
         IF (modfiles(iddesc,ifil,1)%iostat.EQ.3) THEN
            IF (noread(ifil).EQ.2) THEN
               CALL error_file(ierrno_fend,filepars=modfiles(iddesc,ifil,1))
            ELSE
               update = .FALSE.
               GOTO 1001
            ENDIF
         ENDIF
         CALL diff_dates(CStartDateTime,ciodatetime,0,nosecsdat(2))
      ENDDO
   ENDIF

ENDIF

!
!3. Evaluate data at current time
!---------------------------------
!
!---without time interpolation
IF (no_time_series.OR.tlims(3).LT.0) THEN
   SELECT CASE (iddesc)
      CASE (io_disloc,io_disvol,io_discur)
         disvals = disdata(:,:,1)
      CASE DEFAULT
         disvals = MERGE(disdata(:,:,1),float_fill,&
                       & ABS(disdata(:,:,1)-fill_value).GT.fill_value_eps)
    END SELECT

!---with time interpolation
ELSE
   ratio2 = MERGE(0.0,(nosecsrun-nosecsdat(1))&
               & /REAL(nosecsdat(2)-nosecsdat(1)),tlims(3).LT.0)
   ratio1 = 1.0 - ratio2
   SELECT CASE (iddesc)
      CASE (io_disloc,io_disvol)
         disvals = ratio1*disdata(:,:,1) + ratio2*disdata(:,:,2)
      CASE (io_discur)
         ivar_310: DO ivar=1,novars
            IF (ivar.EQ.1) THEN
               disvals(:,1) = ratio1*disdata(:,1,1) + ratio2*disdata(:,1,2)
            ELSEIF (ivar.EQ.2) THEN
               datsin = ratio1*SIN(disdata(:,2,1))+ratio2*SIN(disdata(:,2,2))
               datcos = ratio1*COS(disdata(:,2,1))+ratio2*COS(disdata(:,2,2))
               disvals(:,2) = ATAN2(datsin,datcos)
            ENDIF
         ENDDO ivar_310
      CASE DEFAULT
         ivar_320: DO ivar=1,novars
         l_320: DO l= 1,nodat
            IF (ABS(disdata(l,ivar,1)-fill_value).LE.fill_value_eps) THEN
               disvals(l,ivar) = float_fill
            ELSEIF (ABS(disdata(l,ivar,2)-fill_value).LE.fill_value_eps) THEN
               disvals(l,ivar) = disdata(l,ivar,1)
            ELSE
               disvals(l,ivar) = ratio1*disdata(l,ivar,1) + &
                               & ratio2*disdata(l,ivar,2)
            ENDIF
         ENDDO l_320
         ENDDO ivar_320
   END SELECT
ENDIF

1001 CONTINUE

SELECT CASE (iddesc)
   CASE (io_disloc,io_disvol,io_discur)
      CALL log_timer_out()
   CASE DEFAULT
      CALL log_timer_out(npcc,itm_structs)
END SELECT


RETURN

!---exit code in case of error
1002 nerrs = 1
IF (errchk.AND.nerrs.LE.maxerrors) THEN
   WRITE (ioerr,'(A)') 'Read invalid date/time: '//ciodatetime
   WRITE (ioerr,'(A)') 'Should be no later than: '//CDateTime
ENDIF
CALL error_file(ierrno_read,filepars=modfiles(iddesc,ifil,1))

END SUBROUTINE update_dischr_data

!=================================================================

SUBROUTINE weirs_bstress
!*****************************************************************
!
! *weirs_bstress* bottom stress at weirs and barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling process - bottom_stress
!
! Module calls - Uvar_at_V, Vvar_at_U
!
!*****************************************************************
!
USE currents
USE depths
USE fluxes
USE grid
USE gridpars
USE iopars
USE physpars
USE structures
USE switches
USE array_interp, ONLY: Uvar_at_V, Vvar_at_U
USE time_routines, ONLY: log_timer_in, log_timer_out

!
!*Local variables
!
INTEGER :: i, iiloc, j, jjloc, k, npcc
REAL :: botrat, curatu, curatv, uvelatv, vvelatu

procname(pglev+1) = 'weirs_sink'
CALL log_timer_in(npcc)

!
!1. U-nodes
!----------
!

iiloc_110: DO iiloc=1,numwbaruloc
   i = iwbaruloc(iiloc); j = jwbaruloc(iiloc)

   IF (node2du(i,j).EQ.2) THEN
      
!
!1.1 2-D case
!------------
!
      IF (iopt_bstres_nodim.EQ.2) THEN
         vvelatu = Vvar_at_U(vmvel(i-1:i,j:j+1),i-1,j,nz,1,3)
         curatu = SQRT(umvel(i,j)**2+vvelatu**2)
         IF (iopt_bstres_drag.GT.2) THEN
            botrat =  MAX(deptotatu(i,j)/(enap*zroughatu(i,j)),zbtoz0lim)
            bdragcoefatu(i,j) =  (ckar/LOG(botrat))**2
         ENDIF
         bfricatu(i,j) = bdragcoefatu(i,j)*curatu
         ubstresatu(i,j) = bfricatu(i,j)*umvel(i,j)

!
!1.2 3-D case
!------------
!

      ELSEIF (iopt_bstres_nodim.EQ.3) THEN
         k = 1
         DO WHILE (k.LT.nz.AND.nodeatu(i,j,k).EQ.1)
            k = k + 1
         ENDDO
         vvelatu = Vvar_at_U(vvel(i-1:i,j:j+1,k),i-1,j,k,1,3)
         curatu = SQRT(uvel(i,j,k)**2+vvelatu**2)
         IF (iopt_bstres_drag.GT.2) THEN
            botrat =  MAX(0.5*delzatu(i,j,k)/zroughatu(i,j),zbtoz0lim)
            bdragcoefatu(i,j) =  (ckar/LOG(botrat))**2
         ENDIF
         bfricatu(i,j) = bdragcoefatu(i,j)*curatu
         ubstresatu(i,j) = bfricatu(i,j)*uvel(i,j,k)
      ENDIF

   ENDIF

ENDDO iiloc_110

!
!2. V-nodes
!----------
!

jjloc_210: DO jjloc=1,numwbarvloc
   i = iwbarvloc(jjloc); j = jwbarvloc(jjloc)

   IF (node2dv(i,j).EQ.2) THEN
      
!
!2.1 2-D case
!------------
!
      IF (iopt_bstres_nodim.EQ.2) THEN
         uvelatv = Uvar_at_V(umvel(i:i:1,j-1:j),i,j-1,nz,1,3)
         curatv = SQRT(uvelatv**2+vmvel(i,j)**2)
         IF (iopt_bstres_drag.GT.2) THEN
            botrat =  MAX(deptotatv(i,j)/(enap*zroughatv(i,j)),zbtoz0lim)
            bdragcoefatv(i,j) =  (ckar/LOG(botrat))**2
         ENDIF
         bfricatv(i,j) = bdragcoefatv(i,j)*curatv
         vbstresatv(i,j) = bfricatv(i,j)*vmvel(i,j)

!
!2.2 3-D case
!------------
!

      ELSEIF (iopt_bstres_nodim.EQ.3) THEN
         k = 1
         DO WHILE (k.LT.nz.AND.nodeatv(i,j,k).EQ.1)
            k = k + 1
         ENDDO
         uvelatv = Uvar_at_V(uvel(i:i+1,j-1:j,k),i,j-1,k,1,3)
         curatv = SQRT(uvelatv**2+vvel(i,j,k)**2)
         IF (iopt_bstres_drag.GT.2) THEN
            botrat =  MAX(0.5*delzatv(i,j,k)/zroughatv(i,j),zbtoz0lim)
            bdragcoefatv(i,j) =  (ckar/LOG(botrat))**2
         ENDIF
         bfricatv(i,j) = bdragcoefatv(i,j)*curatv
         vbstresatv(i,j) = bfricatv(i,j)*vvel(i,j,k)

      ENDIF

   ENDIF

ENDDO jjloc_210

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE weirs_bstress

!=================================================================

SUBROUTINE weirs_depth    
!*****************************************************************
!
! *weirs_depth* total water depth at weirs and barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling process - water_depths
!
! Module calls -
!
!*****************************************************************
!
USE currents
USE depths
USE grid
USE iopars
USE physpars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!

INTEGER :: i, ii, iiloc, j, jj, jjloc, npcc
REAL :: wdep


procname(pglev+1) = 'weirs_depth'
CALL log_timer_in(npcc)

!
!1. U-nodes
!----------
!

iiloc_110: DO iiloc=1,numwbaruloc
   i = iwbaruloc(iiloc); j = jwbaruloc(iiloc)
   ii = indexwbaru(iiloc)
   IF (node2du(i,j).EQ.2) THEN
!     ---upwind level
      IF (udvel(i,j).NE.0.0) THEN
         wdep = depmeanatu(i,j) + MERGE(zeta(i-1,j),zeta(i,j),udvel(i,j).GT.0.0)
      ELSE
         wdep = depmeanatu(i,j) + MAX(zeta(i-1,j),zeta(i,j))
      ENDIF
!     ---part above crest
      deptotatu(i,j) = MAX(wdep-wbarcrestu(ii),0.0)
!     ---add orifice section
      IF (oriwidthu(ii).NE.0.0) THEN
         wdep = MIN(wdep,oriwidthu(ii)+oriheightu(ii))
         deptotatu(i,j) = deptotatu(i,j) + wdep - oriheightu(ii)
      ENDIF
      deptotatu(i,j) = MAX(deptotatu(i,j),dmin_fld)
   ELSE
      deptotatu(i,j) = dmin_fld
   ENDIF
   
ENDDO iiloc_110

!
!2. V-nodes
!----------
!

jjloc_210: DO jjloc=1,numwbarvloc
   i = iwbarvloc(jjloc); j = jwbarvloc(jjloc)
   jj = indexwbarv(jjloc)
   IF (node2dv(i,j).EQ.2) THEN
!     ---upwind level
      IF (vdvel(i,j).NE.0.0) THEN
         wdep = depmeanatv(i,j) + MERGE(zeta(i,j-1),zeta(i,j),vdvel(i,j).GT.0.0)
      ELSE
         wdep = depmeanatv(i,j) + MAX(zeta(i,j-1),zeta(i,j))
      ENDIF
!     ---part above crest
      deptotatv(i,j) = MAX(wdep-wbarcrestv(jj),0.0)
!     ---add orifice section 
      IF (oriwidthv(jj).NE.0.0) THEN
         wdep = MIN(wdep,oriwidthv(jj)+oriheightv(jj))
         deptotatv(i,j) = deptotatv(i,j) + wdep - oriheightv(jj)
      ENDIF
      deptotatv(i,j) = MAX(deptotatv(i,j),dmin_fld)
   ELSE
      deptotatv(i,j) = dmin_fld
   ENDIF

ENDDO jjloc_210

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE weirs_depth

!=================================================================

SUBROUTINE weirs_loss
!*****************************************************************
!
! *weirs_loss* energy losses at weirs and barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling process - current_pred, current_2d
!
! Module calls -
!
!*****************************************************************
!
USE currents
USE depths
USE grid
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local Arrays
!
INTEGER :: i, ii, iiloc, j, jj, jjloc, npcc
REAL :: twothird = 2.0/3.0
REAL, DIMENSION(numwbaru) :: wbarelossu_new
REAL, DIMENSION(numwbarv) :: wbarelossv_new


procname(pglev+1) = 'weirs_loss'
CALL log_timer_in(npcc)

!
!1. Weirs at U-nodes
!-------------------
!

iiloc_110: DO iiloc=1,numwbaruloc
   i = iwbaruloc(iiloc); j = jwbaruloc(iiloc)
   ii = indexwbaru(iiloc)
   wbarelossu_new(ii) = 0.0
   IF (node2du(i,j).EQ.2) THEN

!
!1.1 First direction
!-------------------
!

      IF (udvel(i-1,j).GT.0.0) THEN

         IF (deptotatc(i-1,j).GT.wbarcrestu(ii)) THEN

!
!1.1.1 Free flow
!---------------
!
            IF (deptotatc(i,j)/deptotatc(i-1,j).LE.wbarmodlu(ii)) THEN
               wbarelossu_new(ii) = gaccatu(i,j)*ABS(deptotatc(i-1,j)-&
                   & (ABS(udvel(i-1,j)/wbarcoefu(ii))**twothird-wbarcrestu(ii)))

!
!1.1.2 Submerged flow
!--------------------
!
            ELSE
               wbarelossu_new(ii) = gaccatu(i,j)*(1.0-wbarmodlu(ii))*&
             & (udvel(i-1,j)/wbarcoefu(ii)/(deptotatc(i-1,j)-wbarcrestu(ii)))**2

               IF (oriwidthu(ii).GT.0.0) THEN
                  wbarelossu_new(ii) = wbarelossu_new(ii) + 0.5*(udvel(i-1,j)/&
                                & (oricoefu(ii)*oriwidthu(ii)))**2
               ENDIF
            ENDIF

!
!1.1.3 Orifice flow
!------------------
!
         ELSEIF (oriwidthu(ii).GT.0.0.AND.zeta(i-1,j).GE.zeta(i,j)) THEN
            wbarelossu_new(ii) = 0.5*(udvel(i-1,j)&
                          & /(oricoefu(ii)*oriwidthu(ii)))**2
         ENDIF

      ENDIF

!
!1.2 Second direction
!--------------------
!

      IF (udvel(i-1,j).LT.0.0) THEN

         IF (deptotatc(i,j).GT.wbarcrestu(ii)) THEN

!
!1.2.1 Free flow
!---------------
!
            IF ((deptotatc(i-1,j)/deptotatc(i,j)).LE.wbarmodlu(ii)) THEN
               wbarelossu_new(ii) = gaccatu(i,j)*ABS(deptotatc(i,j)-&
                     & (ABS(udvel(i,j)/wbarcoefu(ii))**twothird-wbarcrestu(ii)))

!
!1.2.2 Submerged flow
!--------------------
!
            ELSE
               wbarelossu_new(ii) = gaccatu(i,j)*(1.0-wbarmodlu(ii))*&
                 & (udvel(i,j)/wbarcoefu(ii)/(deptotatc(i,j)-wbarcrestu(ii)))**2

               IF (oriwidthu(ii).GT.0.0) THEN
                  wbarelossu_new(ii) = wbarelossu_new(ii) + 0.5*(udvel(i,j)/&
                                & (oricoefu(ii)*oriwidthu(ii)))**2
               ENDIF
            ENDIF

!
!1.2.3 Orifice flow
!------------------
!
         ELSEIF (oriwidthu(ii).GT.0.0.AND.zeta(i-1,j).LT.zeta(i,j)) THEN
            wbarelossu_new(ii) = 0.5*(udvel(i,j)/&
                               & (oricoefu(ii)*oriwidthu(ii)))**2
         ENDIF

      ENDIF

!
!1.2.4 Time relaxation
!---------------------
!

      wbarelossu(ii) = (1.0-wbarrlxu)*wbarelossu(ii) + &
                      & wbarrlxu*wbarelossu_new(ii)

   ENDIF

ENDDO iiloc_110

!
!2. Weirs at V-nodes
!-------------------
!

jjloc_210: DO jjloc=1,numwbarvloc
   i = iwbarvloc(jjloc); j = jwbarvloc(jjloc)
   jj = indexwbarv(jjloc)
   wbarelossv_new(jj) = 0.0
   IF (node2dv(i,j).EQ.2) THEN

!
!2.1 First direction
!-------------------
!

      IF (vdvel(i,j-1).GE.0.0) THEN

         IF (deptotatc(i,j-1).GT.wbarcrestv(jj)) THEN

!
!2.1.1 Free flow
!---------------
!

            IF ((deptotatc(i,j)/deptotatc(i,j-1)).LE.wbarmodlv(jj)) THEN
               wbarelossv_new(jj) = gaccatv(i,j)*ABS(deptotatc(i,j-1)-&
                   & (ABS(vdvel(i,j-1)/wbarcoefv(jj))**twothird-wbarcrestv(jj)))

!
!2.1.2 Submerged flow
!--------------------
!
            ELSE
               wbarelossv_new(jj) = gaccatv(i,j)*(1.0-wbarmodlv(jj))*&
             & (vdvel(i,j-1)/wbarcoefv(jj)/(deptotatc(i,j-1)-wbarcrestv(jj)))**2

               IF (oriwidthv(jj).GT.0.0) THEN
                  wbarelossv_new(jj) = wbarelossv_new(jj) + 0.5*(vdvel(i,j-1)/&
                                & (oricoefv(jj)*oriwidthv(jj)))**2
               ENDIF
            ENDIF

!
!2.1.3 Orifice flow
!------------------
! 
         ELSEIF (oriwidthv(jj).GT.0.0.AND.zeta(i,j-1).GE.zeta(i,j)) THEN
            wbarelossv_new(jj) = 0.5*(vdvel(i,j-1)/&
                               & (oricoefv(jj)*oriwidthv(jj)))**2
         ENDIF
  
      ENDIF

!
!2.2 Second direction
!--------------------
!

      IF (vdvel(i,j-1).LT.0.0) THEN

         IF (deptotatc(i,j).GT.wbarcrestv(jj)) THEN
!
!2.2.1 Free flow
!---------------
!
            IF ((deptotatc(i,j-1)/deptotatc(i,j)).LE.wbarmodlv(jj)) THEN
               wbarelossv_new(jj) = gaccatv(i,j)*ABS(deptotatc(i,j)-&
                     & (ABS(vdvel(i,j)/wbarcoefv(jj))**twothird-wbarcrestv(jj)))

!
!2.2.2 Submerged flow
!--------------------
!
            ELSE
               wbarelossv_new(jj) = gaccatv(i,j)*(1.0-wbarmodlv(jj))*&
                 & (vdvel(i,j)/wbarcoefv(jj)/(deptotatc(i,j)-wbarcrestv(jj)))**2

               IF (oriwidthv(jj).GT.0.0) THEN
                  wbarelossv_new(jj) = wbarelossv_new(jj) + 0.5*(vdvel(i,j)/&
                                & (oricoefv(jj)*oriwidthv(jj)))**2
               ENDIF
            ENDIF

!
!2.2.3 Orifice flow
!------------------
! 
         ELSEIF (oriwidthv(jj).GT.0.0.AND.zeta(i,j-1).LT.zeta(i,j)) THEN
            wbarelossv_new(jj) = 0.5*(vdvel(i,j)&
                              & /(oricoefv(jj)*oriwidthv(jj)))**2
         ENDIF

      ENDIF

!
!2.2. Time relaxation
!--------------------
!

      wbarelossv(jj) = (1.0-wbarrlxv)*wbarelossv(jj) + &
                     & wbarrlxv*wbarelossv_new(jj)

   ENDIF

ENDDO jjloc_210

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE weirs_loss

!=================================================================

SUBROUTINE weirs_mask    
!*****************************************************************
!
! *weirs_mask* update velocity pointers at weirs/barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling process - coherens_main
!
! External calls - update_pointer_arrays
!
! Module calls - 
!
!*****************************************************************
!
USE currents
USE depths
USE grid
USE gridpars
USE iopars
USE structures
USE switches
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Local variables
!
LOGICAL :: flowdir
INTEGER :: i, ic, idir, ii, iiloc, j, jc, jdir, jj, jjloc, klo, kup, kw, npcc
REAL :: sigcrest, siglo, sigup


procname(pglev+1) = 'weirs_mask'
CALL log_timer_in(npcc)

!
!1. U-nodes
!----------
!

iiloc_110: DO iiloc=1,numwbaruloc
   i = iwbaruloc(iiloc); j = jwbaruloc(iiloc)
   ii = indexwbaru(iiloc)
   idir_111: DO idir=1,2
      flowdir = MERGE(udvel(i-1,j).GE.0.0,udvel(i+1,j).LT.0.0,idir.EQ.1)
      IF (flowdir) THEN
         ic = MERGE(i-1,i,idir.EQ.1)
         sigcrest = wbarcrestu(ii)/deptotatc(ic,j)

!
!1.1 2-D case
!------------
!
         IF (iopt_grid_nodim.EQ.2) THEN
            nodeatu(i,j,:) = MERGE(1,2,sigcrest.GE.1)

!
!1.2 3-D case
!------------
!

         ELSE

!
!1.2.1 Weirs
!-----------
!

            IF (sigcrest.GE.1.0) THEN
               kw = nz+1
            ELSE
               kw = 1
               DO WHILE (kw.LE.nz.AND.gscoordatw(ic,j,kw+1).LT.sigcrest)
                  kw = kw + 1
               END DO
               IF (gscoordatc(ic,j,kw).LE.sigcrest) kw = kw + 1
            ENDIF
            IF (kw.GT.1) nodeatu(i,j,1:kw-1) = 1
            IF (kw.LE.nz) nodeatu(i,j,kw:nz) = 2

!
!1.2.2 Orifices
!--------------
!

            IF (oriwidthu(ii).NE.0.0) THEN
!              ---lower
               siglo = oriheightu(ii)/deptotatc(ic,j)
               IF (siglo.GE.1.0) THEN
                  klo = nz+1
               ELSE
                  klo = 1
                  DO WHILE (klo.LE.nz.AND.gscoordatw(ic,j,klo+1).LT.siglo)
                     klo = klo + 1
                  END DO
                  IF (gscoordatc(ic,j,klo).LE.siglo) klo = klo + 1
               ENDIF
!              ---upper
               sigup = (oriheightu(ii)+oriwidthu(ii))/deptotatc(ic,j)
               IF (sigup.GE.1.0) THEN
                  kup = nz
               ELSE
                  kup = klo
                  DO WHILE (kup.LE.nz.AND.gscoordatw(ic,j,kup+1).LT.sigup)
                     kup = kup + 1
                  END DO
                  IF (gscoordatc(ic,j,kup).GT.sigup) kup = kup - 1
               ENDIF
               kup = MIN(kup,kw-1)
               IF (klo.GT.1) nodeatu(i,j,1:klo-1) = 1
               IF (klo.LE.kup) nodeatu(i,j,klo:kup) = 2
            ENDIF
            
         ENDIF
      
      ENDIF

   ENDDO idir_111
   
ENDDO iiloc_110

!
!2. V-nodes
!----------
!

jjloc_210: DO jjloc=1,numwbarvloc
   i = iwbarvloc(jjloc); j = jwbarvloc(jjloc)
   jj = indexwbarv(jjloc)
   jdir_211: DO jdir=1,2
      flowdir = MERGE(vdvel(i,j-1).GT.0.0,vdvel(i,j+1).LT.0.0,idir.EQ.1)
      IF (flowdir) THEN
         jc = MERGE(j-1,j,jdir.EQ.1)
         sigcrest = wbarcrestv(jj)/deptotatc(i,jc)

!
!1.1 2-D case
!------------
!
         IF (iopt_grid_nodim.EQ.2) THEN
            nodeatv(i,j,:) = MERGE(1,2,sigcrest.GE.1)

!
!1.2 3-D case
!------------
!

         ELSE

!
!1.2.1 Weirs
!-----------
!

            IF (sigcrest.GE.1.0) THEN
               kw = nz+1
            ELSE
               kw = 1
               DO WHILE (kw.LE.nz.AND.gscoordatw(i,jc,kw+1).LT.sigcrest)
                  kw = kw + 1
               END DO
               IF (gscoordatc(i,jc,kw).LE.sigcrest) kw = kw + 1
            ENDIF
            IF (kw.GT.1) nodeatv(i,j,1:kw-1) = 1
            IF (kw.LE.nz) nodeatv(i,j,kw:nz) = 2
            
!
!1.2.2 Orifices
!--------------
!

            IF (oriwidthv(jj).NE.0.0) THEN
!              ---lower
               siglo = oriheightv(jj)/deptotatc(i,jc)
               IF (siglo.GE.1.0) THEN
                  klo = nz+1
               ELSE
                  klo = 1
                  DO WHILE (klo.LE.nz.AND.gscoordatw(i,jc,klo+1).LT.siglo)
                     klo = klo + 1
                  END DO
                  IF (gscoordatc(i,jc,klo).LE.siglo) klo = klo + 1
               ENDIF
!              ---upper
               sigup = (oriheightv(jj)+oriwidthv(jj))/deptotatc(i,jc)
               IF (sigup.GE.1.0) THEN
                  kup = nz
               ELSE
                  kup = klo
                  DO WHILE (kup.LE.nz.AND.gscoordatw(i,jc,kup+1).LT.sigup)
                     kup = kup + 1
                  END DO
                  IF (gscoordatc(i,jc,kup).GT.sigup) kup = kup - 1
               ENDIF
               kup = MIN(kup,kw-1)
               IF (klo.GT.1) nodeatv(i,j,1:klo-1) = 1
               IF (klo.LE.kup) nodeatv(i,j,klo:kup) = 2
            ENDIF
            
         ENDIF

      ENDIF
      
   ENDDO jdir_211
   
ENDDO jjloc_210

!
!3. Update pointer arrays
!------------------------
!

CALL update_pointer_arrays

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE weirs_mask

!=================================================================

SUBROUTINE weirs_sink(wsinkatu,wsinkatv,nzdim)   
!*****************************************************************
!
! *weirs_sink*  sink terms in 2-D/3-D momentum equations by energy loss
!               through weirs/barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.6
!
! Description -
!
! Reference -
!
! Calling process - current_pred, current_2d
!
! Module calls -
!
!*****************************************************************
!
USE currents
USE grid
USE gridpars
USE iopars
USE structures
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: nzdim
REAL, INTENT(INOUT), DIMENSION(ncloc,nrloc,nzdim) :: wsinkatu, wsinkatv

!
! Name      Type    Purpose
!------------------------------------------------------------------------------
!*nzdim*    INTEGER Vertical dimension (1 for 2-D, nz for 3-D case)
!*wsinkatu* REAL    Sink term in the U-momentum equation                   [1/s]
!*wsinkatv* REAL    Sink term in the V-momentum equation                   [1/s]
!
!------------------------------------------------------------------------------
!
!*Local Arrays
!
INTEGER :: i, ii, iiloc, j, jj, jjloc, k, npcc
REAL :: denom


procname(pglev+1) = 'weirs_sink'
CALL log_timer_in(npcc)

!
!1. U-nodes
!----------
!

iiloc_110: DO iiloc=1,numwbaruloc
   i = iwbaruloc(iiloc); j = jwbaruloc(iiloc)
   ii = indexwbaru(iiloc)
   denom = delxatu(i,j)*MAX(0.001,ABS(umvel(i,j)))
   k_111: DO k=1,nzdim
      IF (nodeatu(i,j,k).EQ.2) THEN
         wsinkatu(i,j,k) = wsinkatu(i,j,k) + wbarelossu(ii)/denom
      ENDIF
   END DO k_111
ENDDO iiloc_110

!
!2. V-nodes
!----------
!

jjloc_210: DO jjloc=1,numwbarvloc
   i = iwbarvloc(jjloc); j = jwbarvloc(jjloc)
   jj = indexwbarv(jjloc)
   denom = delyatv(i,j)*MAX(0.001,ABS(vmvel(i,j)))
   k_211: DO k=1,nzdim
      IF (nodeatv(i,j,k).EQ.2) THEN
         wsinkatv(i,j,k) = wsinkatv(i,j,k) + wbarelossv(jj)/denom
      ENDIF
   END DO k_211
ENDDO jjloc_210

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE weirs_sink

!=================================================================

SUBROUTINE write_dischr_data(iddesc,ifil,ciodatetime,disdata,nodat,novars)
!*****************************************************************
!
! *write_dischr_data* Write discharge data in standard format
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12.1
!
! Description -
!
! Reference -
!
! Calling program - define_dischr_data
!
! Module calls - error_alloc_struc, open_filepars, set_modfiles_atts,
!                set_modvars_atts, varatts_init, write_atts_mod, write_time,
!                write_vars
!
!*****************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE syspars
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: open_filepars, write_atts_mod, write_time, write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(IN) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, nodat, novars
REAL, INTENT(IN), DIMENSION(nodat,novars) :: disdata

!
! Name         Type Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER Data file id
!*ifil*        INTEGER No. of data file
!*ciodatetime* CHAR    Date/time in data file
!*disdata*     REAL    Input data
!*nodat*       INTEGER Number of discharge locations
!*novars*      INTEGER Number of input variables
!
!------------------------------------------------------------------------------
!
!* Local variables
!
INTEGER, SAVE :: maxvars = 4, nodesc = 5
INTEGER :: idesc, ivar, numvars
INTEGER, DIMENSION(novars) :: vecids
TYPE (FileParams) :: filepars
TYPE (VariableAtts), SAVE, ALLOCATABLE, DIMENSION(:,:) :: varatts


IF (.NOT.master) RETURN 

procname(pglev+1) = 'write_dischr_data'
CALL log_timer_in()

SELECT CASE (iddesc)
   CASE (io_disloc); idesc = 1
   CASE (io_disvol); idesc = 2
   CASE (io_discur); idesc = 3
   CASE (io_dissal); idesc = 4
   CASE (io_distmp); idesc = 5
END SELECT

filepars = modfiles(iddesc,ifil,2)
numvars = novars + 1

!
!1. Write header on first call
!-----------------------------
!

IF (filepars%iostat.EQ.0) THEN
   
!  ---file attributes
   CALL set_modfiles_atts(iddesc,ifil,2,filepars)

!  ---open file
   CALL open_filepars(filepars)
   
!  ---allocate
   IF (outflagdis) THEN
      IF (ALLOCATED(varatts)) DEALLOCATE(varatts)
      ALLOCATE (varatts(maxvars,nodesc),STAT=errstat)
      CALL error_alloc_struc('varatts',2,(/maxvars,nodesc/),'VariableAtts')
      CALL varatts_init(varatts)
      outflagdis = .FALSE.
   ENDIF
   
!  ---variable attributes
   CALL set_modvars_atts(iddesc,ifil,2,filepars,numvars,&
                       & varatts(1:numvars,idesc))

!  ---write atrtributes
   CALL write_atts_mod(filepars,varatts(1:numvars,idesc))
   
ENDIF
 
!
!2. Write data of discharges
!---------------------------
!
!---date/time
CALL write_time(ciodatetime,filepars)

!---data
vecids = (/(ivar,ivar=2,numvars)/)
CALL write_vars(disdata,filepars,0,varatts(2:numvars,idesc),vecids=vecids)

modfiles(iddesc,ifil,2) = filepars

CALL log_timer_out()


RETURN

END SUBROUTINE write_dischr_data

!=================================================================

SUBROUTINE write_dry_cells
!*****************************************************************
!
! *write_dry_cells* Write external file for dry cells
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12.1
!
! Description - write external files for dry cells 
!
! Reference -
!
! Calling process - define_global_grid
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, varatts_init,
!                write_atts_mod, write_vars
!
!*****************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE structures
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: ivar, npcc, numvars
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN 

procname(pglev+1) = 'write_dry_cells'
CALL log_timer_in(npcc)

!
!1. Write File Header
!--------------------
!
!---file attributes
filepars = modfiles(io_drycel,1,2)
CALL set_modfiles_atts(io_drycel,1,2,filepars)

numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL varatts_init(varatts)
CALL set_modvars_atts(io_drycel,1,2,filepars,numvars,varatts)

!---write
CALL write_atts_mod(filepars,varatts)
 
!
!2. Write coordinates of dry cells
!---------------------------------
!

ivar_210: DO ivar=1,numvars
   SELECT CASE (TRIM(varatts(ivar)%f90_name))
      CASE ('idry')
         CALL write_vars(idry,filepars,ivar,(/varatts(ivar)/))
      CASE ('jdry')
         CALL write_vars(jdry,filepars,ivar,(/varatts(ivar)/))
   END SELECT
ENDDO ivar_210

!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_drycel,1,2) = filepars
DEALLOCATE(varatts)

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE write_dry_cells

!=================================================================

SUBROUTINE write_thin_dams
!*****************************************************************
!
! *write_thin_dams* Write external file for thin dams
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12.1
!
! Description - write external files for thin dams  
!
! Reference -
!
! Calling process - initialise_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, write_atts_mod, write_vars
!
!*****************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE structures
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out

IMPLICIT NONE

!
!* Local variables
!
INTEGER :: ivar, npcc, numvars
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN 

procname(pglev+1) = 'write_thin_dams'
CALL log_timer_in(npcc)

!
!1. Write File Header
!--------------------
!
!---file attibutes
filepars = modfiles(io_thndam,1,2)
CALL set_modfiles_atts(io_thndam,1,2,filepars)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL varatts_init(varatts)
CALL set_modvars_atts(io_thndam,1,2,filepars,numvars,varatts)

!---write
CALL write_atts_mod(filepars,varatts)
 
!
!2. Write coordinates of thin dams
!---------------------------------
!

ivar_210: DO ivar=1,numvars
   SELECT CASE (TRIM(varatts(ivar)%f90_name))
      CASE ('ithinu')
         CALL write_vars(ithinu,filepars,ivar,(/varatts(ivar)/))
      CASE ('jthinu')
         CALL write_vars(jthinu,filepars,ivar,(/varatts(ivar)/))
      CASE ('ithinv')
         CALL write_vars(ithinv,filepars,ivar,(/varatts(ivar)/))
      CASE ('jthinv')
         CALL write_vars(jthinv,filepars,ivar,(/varatts(ivar)/))
   END SELECT
ENDDO ivar_210
   
!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_thndam,1,2) = filepars
DEALLOCATE(varatts)

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE write_thin_dams

!=================================================================

SUBROUTINE write_weirs
!*****************************************************************
!
! *write_weirs* Write external file for weirs/barriers
!
! Author - ANTEA
!
! Version - @(COHERENS)Structures_Model.f90 V2.12.1
!
! Description - write external files for weirs/barriers 
!
! Reference -
!
! Calling process - initialise_model
!
! Module calls - close_filepars, error_alloc_struc, open_filepars,
!                set_modfiles_atts, set_modvars_atts, varatts_init, write_vars
!
!*****************************************************************
!
USE datatypes
USE iopars
USE paralpars
USE structures
USE datatypes_init, ONLY: varatts_init
USE error_routines, ONLY: error_alloc_struc
USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                        & write_vars
USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
USE time_routines, ONLY: log_timer_in, log_timer_out


IMPLICIT NONE

!
!* Local variables
!
INTEGER :: ivar, npcc, numvars
TYPE (FileParams) :: filepars
TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts


IF (.NOT.master) RETURN 

procname(pglev+1) = 'write_weirs'
CALL log_timer_in(npcc)

!
!1. Write File Header
!--------------------
!
!---file attibutes
filepars = modfiles(io_weibar,1,2)
CALL set_modfiles_atts(io_weibar,1,2,filepars)
numvars = filepars%novars

!---open file
CALL open_filepars(filepars)

!---variable attributes
ALLOCATE (varatts(numvars),STAT=errstat)
CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
CALL varatts_init(varatts)
CALL set_modvars_atts(io_weibar,1,2,filepars,numvars,varatts)

!---write
CALL write_atts_mod(filepars,varatts)
 
!
!2. Write data of weirs/barriers
!-------------------------------
!
ivar_210: DO ivar=1,numvars
   SELECT CASE (TRIM(varatts(ivar)%f90_name))
      CASE ('iwbaru')
         CALL write_vars(iwbaru,filepars,ivar,(/varatts(ivar)/))
      CASE ('jwbaru')
         CALL write_vars(jwbaru,filepars,ivar,(/varatts(ivar)/))
      CASE ('oricoefu')
         CALL write_vars(oricoefu,filepars,ivar,(/varatts(ivar)/))
      CASE ('oriheightu')
         CALL write_vars(oriheightu,filepars,ivar,(/varatts(ivar)/))
      CASE ('oriwidthu')
         CALL write_vars(oriwidthu,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarcoefu')
         CALL write_vars(wbarcoefu,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarcrestu')
         CALL write_vars(wbarcrestu,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarmodlu')
         CALL write_vars(wbarmodlu,filepars,ivar,(/varatts(ivar)/))
      CASE ('iwbarv')
         CALL write_vars(iwbarv,filepars,ivar,(/varatts(ivar)/))
      CASE ('jwbarv')
         CALL write_vars(jwbarv,filepars,ivar,(/varatts(ivar)/))
      CASE ('oricoefv')
         CALL write_vars(oricoefv,filepars,ivar,(/varatts(ivar)/))
      CASE ('oriheightv')
         CALL write_vars(oriheightv,filepars,ivar,(/varatts(ivar)/))
      CASE ('oriwidthv')
         CALL write_vars(oriwidthv,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarcoefv')
         CALL write_vars(wbarcoefv,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarcrestv')
         CALL write_vars(wbarcrestv,filepars,ivar,(/varatts(ivar)/))
      CASE ('wbarmodlv')
         CALL write_vars(wbarmodlv,filepars,ivar,(/varatts(ivar)/))

   END SELECT

ENDDO ivar_210
   
!
!3. Finalise
!-----------
!

CALL close_filepars(filepars)
modfiles(io_weibar,1,2) = filepars
DEALLOCATE(varatts)

CALL log_timer_out(npcc,itm_structs)


RETURN

END SUBROUTINE write_weirs

!! AC 16.11.2022

SUBROUTINE write_mpv
   !*****************************************************************
   !
   ! *write_mpv* Write external file for mpv
   !
   ! Author - A. Capet
   !
   ! Version - @(COHERENS)Structures_Model.f90 V2.12.1
   !
   ! Description - write external files for mpv 
   !
   ! Reference -
   !
   ! Calling process - define_global_grid
   !
   ! Module calls - close_filepars, error_alloc_struc, open_filepars,
   !                set_modfiles_atts, set_modvars_atts, varatts_init,
   !                write_atts_mod, write_vars
   !
   !*****************************************************************
   !
   USE datatypes
   USE iopars
   USE paralpars
   USE structures
   USE datatypes_init, ONLY: varatts_init
   USE error_routines, ONLY: error_alloc_struc
   USE inout_routines, ONLY: close_filepars, open_filepars, write_atts_mod, &
                           & write_vars
   USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
   USE time_routines, ONLY: log_timer_in, log_timer_out
   
   IMPLICIT NONE
   
   !
   !* Local variables
   !
   INTEGER :: ivar, npcc, numvars
   TYPE (FileParams) :: filepars
   TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts
   
   
   IF (.NOT.master) RETURN 
   
   procname(pglev+1) = 'write_mpv'
   CALL log_timer_in(npcc)
   
   !
   !1. Write File Header
   !--------------------
   !
   !---file attributes
   filepars = modfiles(io_mpvcov,1,2)
   CALL set_modfiles_atts(io_mpvcov,1,2,filepars)
   
   numvars = filepars%novars
   
   !---open file
   CALL open_filepars(filepars)
   
   !---variable attributes
   ALLOCATE (varatts(numvars),STAT=errstat)
   CALL error_alloc_struc('varatts',1,(/numvars/),'VariableAtts')
   CALL varatts_init(varatts)
   CALL set_modvars_atts(io_mpvcov,1,2,filepars,numvars,varatts)
   
   !---write
   CALL write_atts_mod(filepars,varatts)
    
   !
   !2. Write partial cover of MPV
   !---------------------------------
   !
   ! AC 16.11.2022 ('The function write_mpv in Structures_Model.f90 is not completed')
   ivar_210: DO ivar=1,numvars
      SELECT CASE (TRIM(varatts(ivar)%f90_name))
         CASE ('idry')
            CALL write_vars(idry,filepars,ivar,(/varatts(ivar)/))
         CASE ('jdry')
            CALL write_vars(jdry,filepars,ivar,(/varatts(ivar)/))
      END SELECT
   ENDDO ivar_210
   
   !
   !3. Finalise
   !-----------
   !
   
   CALL close_filepars(filepars)
   modfiles(io_mpvcov,1,2) = filepars
   DEALLOCATE(varatts)
   
   CALL log_timer_out(npcc,itm_structs)
   
   
   RETURN
   
   END SUBROUTINE write_mpv

SUBROUTINE read_mpv
   !************************************************************************
   !
   ! *read_mpv* Read model mpv arrays in standard format
   !
   ! Author - A. Capet
   !
   ! Version - @(COHERENS)Grid_Arrays.f90  V2.12.1
   !
   ! Description -
   !
   ! Reference -
   !
   ! Calling program - initialise_model
   !
   ! Module calls - check_value_varatts, close_filepars, error_abort,
   !                error_alloc_struc, error_value_arr_struc, open_filepars,
   !                read_glbatts_mod, read_varatts_mod, read_vars,
   !                set_modfiles_atts, set_modvars_atts, varatts_init
   !
   !************************************************************************
  

   USE datatypes
   USE structures
   USE depths
   USE grid
   USE gridpars
   USE iopars
   USE syspars
   USE timepars
   USE check_model, ONLY: check_value_varatts
   USE datatypes_init, ONLY: varatts_init
   USE error_routines, ONLY: error_abort, error_alloc_struc, error_value_arr_struc
   USE inout_routines, ONLY: close_filepars, open_filepars, read_glbatts_mod, &
                           & read_varatts_mod, read_vars
   USE modvars_routines, ONLY: set_modfiles_atts, set_modvars_atts
   USE time_routines, ONLY: log_timer_in, log_timer_out

   USE inout_paral, ONLY: read_distribute_mod
   
   IMPLICIT NONE

   !*Local variables
   
   LOGICAL :: header
   INTEGER, SAVE :: novars, numvars 
   INTEGER :: ivar, nocoords1, nocoords2, novars1, novars2, numvars1, numvars2
   TYPE (FileParams) :: filepars, filepars1, filepars2
   TYPE (VariableAtts), ALLOCATABLE, DIMENSION(:) :: varatts, varatts1, varatts2

   INTEGER, DIMENSION(2) :: lbounds, ubounds
   INTEGER, DIMENSION(4) :: nhdims

   procname(pglev+1) = 'read_mpv'
   CALL log_timer_in()
 
   !1. Initialise
   !1.1 Data attributes    
   !-------------------
   !  
   header = modfiles(io_mpvcov,1,1)%header

   IF (header) THEN

      !  ---open data file
      filepars1 = modfiles(io_mpvcov,1,1)
      CALL open_filepars(filepars1)
      
      !  ---file attributes
      CALL read_glbatts_mod(filepars1)

      !  ---variable attributes
      nocoords1 = filepars1%nocoords
      novars1 = filepars1%novars

      numvars1 = nocoords1 + novars1

      !! AC 09.12.2022
      !! nocoords is read from the netcdf file, which includes coordinates variables.
      !! Yet for some reason, the check protocols assumes that nocoord > 0 automatially involves a time dimensions
      !! For now, I'll attempt forcing nocoords to zero
      nocoords1 = 0
      
      ALLOCATE (varatts1(numvars1),STAT=errstat)
      CALL error_alloc_struc('varatts1',1,(/numvars1/),'VariableAtts')
      CALL varatts_init(varatts1)
      CALL read_varatts_mod(filepars1,varatts1)
     
   ENDIF
     
!
!!1.2 Model attributes
!---------------------
!
   !---file attributes
   filepars2 = modfiles(io_mpvcov,1,1)
   CALL set_modfiles_atts(io_mpvcov,1,1,filepars2)

   !
   !---data attributes   
   nocoords2 = filepars2%nocoords

   novars2 = filepars2%novars

   numvars2 = nocoords2 + novars2
   
   ALLOCATE (varatts2(numvars2),STAT=errstat)
   CALL error_alloc_struc('varatts2',1,(/numvars2/),'VariableAtts')
   CALL varatts_init(varatts2)
   CALL set_modvars_atts(io_mpvcov,1,1,filepars2,numvars2,varatts2)

   !1.3 Check attributes
   
   IF (header) THEN
      CALL error_value_arr_struc(nocoords1,'modfiles','nocoords',nocoords2,3,&
                              & (/io_mpvcov,1,1/))
      CALL check_value_varatts(io_mpvcov,1,varatts1,varatts2,filepars1,filepars2)
      CALL error_abort('read_mpv',ierrno_read)
   ENDIF

   !1.4 Variables names used for idenfication

   IF (header) THEN
      novars = novars1 ; numvars= numvars1
      ALLOCATE (varatts(novars),STAT=errstat)
      CALL error_alloc_struc('varatts',1,(/novars/),'VariableAtts')
      CALL varatts_init(varatts)
      varatts = varatts1
      DEALLOCATE (varatts1)
      filepars = filepars1
   ELSE
      novars = novars2
      ALLOCATE (varatts(novars),STAT=errstat)
      CALL error_alloc_struc('varatts',1,(/novars/),'VariableAtts')
      CALL varatts_init(varatts1)
      varatts = varatts2
      filepars = filepars2
      CALL open_filepars(filepars)
   ENDIF

   !1.5 Deallocate
   !--------------
   !
   DEALLOCATE (varatts2)
!
!2. Read data
!------------
!
   
   ivar_210: DO ivar=1,numvars
      SELECT CASE (TRIM(varatts(ivar)%f90_name))                                                                                                                                                           
      CASE ('mpvcov')
         nhdims = 0
         lbounds = 1
         ubounds = (/nc,nr/)
         CALL read_distribute_mod(mpvcov,filepars,ivar,&
              & lbounds,ubounds,nhdims,'C  ',.TRUE.,&
              & (/varatts(ivar)/))

      END SELECT

   ENDDO ivar_210
   
   !3. Finalise
   !---close file
   CALL close_filepars(filepars)
   modfiles(io_mpvcov,1,1) = filepars

   !---deallocate    
   DEALLOCATE (varatts)

   CALL log_timer_out()

   RETURN
   
   END SUBROUTINE read_mpv
