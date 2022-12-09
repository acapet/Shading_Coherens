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
! *Usrdef_Surface_Data* User-defined surface data setup
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! $Date: 2020-02-04 08:57:32 +0000 (Die, 04 Feb 2020) $
!
! $Revision: 87 $
!
! Description - North Sea model (NoS)
!
! Reference -
!
! Subroutines - usrdef_surface_absgrd, usrdef_surface_relgrd,
!               usrdef_surface_data
!
!************************************************************************
!
!Scale and offset function

FUNCTION scale_and_offset(file_id, var_id)
!look in the netcdf file if there is a scale factor(1) and a offset(2)
USE netcdf

IMPLICIT NONE

INTEGER, INTENT(IN) :: file_id
INTEGER, INTENT(IN) :: var_id
REAL, DIMENSION(2) :: scale_and_offset
INTEGER :: status
status = nf90_get_att(file_id,var_id, "scale_factor", scale_and_offset(1))
if (status /= nf90_noerr) scale_and_offset(1) = 1
status = nf90_get_att(file_id,var_id, "add_offset", scale_and_offset(2))
if (status /= nf90_noerr) scale_and_offset(2) = 0

END FUNCTION

!========================================================================

SUBROUTINE usrdef_surface_absgrd(iddesc,ifil,n1dat,n2dat,xcoord,ycoord)
!************************************************************************
!
! *usrdef_surface_absgrd* Define coordinate arrays of surface grid(s)
!                         with respect to model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
!************************************************************************
!
USE datatypes
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, n1dat, n2dat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat) :: xcoord
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat) :: ycoord

!
! Name       Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*    INTEGER Grid file id
!*ifil*      INTEGER No. of grid file
!*n1dat*     INTEGER X-dimension of data grid
!*n2dat*     INTEGER Y-dimension of data grid
!*xcoord*    REAL    X-coordinates of data grid
!*ycoord*    REAL    Y-coordinates of data grid
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_surface_absgrd

!========================================================================

SUBROUTINE usrdef_surface_relgrd(iddesc,ifil,surfgridglb,nx,ny)
!************************************************************************
!
! *usrdef_surface_relgrd* Define relative coordinate array of surface grid(s)
!                         with respect to model grid
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description -
!
! Reference -
!
! Calling program - define_surface_input_grid, define_surface_output_grid
!
!************************************************************************
!
USE datatypes
USE gridpars

IMPLICIT NONE

!
!*Arguments
!
INTEGER, INTENT(IN) :: iddesc, ifil, nx, ny
TYPE (HRelativeCoords), INTENT(INOUT), DIMENSION(nx,ny) :: surfgridglb

!
! Name         Type    Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER Grid file id
!*ifil*        INTEGER No. of grid file
!*surfgridglb* DERIVED Relative coordinates of model grid to data grid or data
!                      grid to model grid
!*nx*          INTEGER X-dimension of data grid
!*ny*          INTEGER Y-dimension of data grid 
!
!------------------------------------------------------------------------------
!


RETURN

END SUBROUTINE usrdef_surface_relgrd

!========================================================================

SUBROUTINE usrdef_surface_data(iddesc,ifil,ciodatetime,surdata,n1dat,n2dat,&
                             & novars)
!************************************************************************
!
! *usrdef_surface_data* Define surface input data
!
! Author - Patrick Luyten
!
! Version - @(COHERENS)Usrdef_Surface_Data.f90  V2.0
!
! Description - North Sea model (NoS)
!
! Reference -
!
! Calling program - define_surface_data
!
!************************************************************************
!

USE iopars
USE paralpars
USE switches
USE error_routines, ONLY: error_alloc
USE inout_routines, ONLY: open_filepars
USE timepars
USE time_routines, ONLY: add_secs_to_date_int, convert_date, log_timer_in, &
                       & log_timer_out
USE syspars
USE physpars
!USE cf90_routines
USE netcdf

IMPLICIT NONE

!
!*Arguments
!
CHARACTER (LEN=lentime), INTENT(INOUT) :: ciodatetime
INTEGER, INTENT(IN) :: iddesc, ifil, novars, n1dat, n2dat
REAL, INTENT(INOUT), DIMENSION(n1dat,n2dat,novars) :: surdata

!
! Name         Type     Purpose
!------------------------------------------------------------------------------
!*iddesc*      INTEGER  Data file id
!*ifil*        INTEGER  No. of data file
!*ciodatetime* CHAR     Date/time in data file
!*surdata*     REAL     Data array
!*n1dat*       INTEGER  X-dimension of data array
!*n2dat*       INTEGER  Y-dimension of data array
!*novars*      INTEGER  Number of data parameters
!
!------------------------------------------------------------------------------
!
!*Local variables
!
INTEGER :: status
INTEGER, SAVE :: iunit,iunit2,chunksize,xrecs,yrecs
INTEGER, PARAMETER :: nodat = 9
INTEGER :: beginperiod, endperiod, i, irange, ivar, j
INTEGER, DIMENSION(7) :: iodatetime1, iodatetime2
REAL :: period, vpres
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: precip, qspec,dewt
! AC 14102022 : considering tcc instead of lcc, mcc, hcc '
!REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: clouddat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: clouddat
! ! ! ! ! ! 
INTEGER(KIND = 8), SAVE, ALLOCATABLE, DIMENSION(:,:) :: help_dat
REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: meteodat
INTEGER(KIND=8), SAVE, ALLOCATABLE, DIMENSION(:,:) :: udat
REAL, PARAMETER :: r1 = 0.62197, r2 = 1.0-r1
!DVDE
INTEGER, DIMENSION(9) :: flag
INTEGER :: totflag
CHARACTER(len=1) :: cvar
CHARACTER(len=5) :: typefct
CHARACTER(len=7) :: crange
CHARACTER(len=30) :: varname
!DVDE
CHARACTER(len=80) :: line

!parameters related to reading netcdf files:
TYPE(FileParams), SAVE :: filepars
INTEGER,SAVE :: ndims, nvars, numvars, nglobatts,norecords
INTEGER,SAVE :: dim_timeID, maxrecs,iyear,imont,iday,ihour,nextmonth,nextyear
INTEGER,SAVE ::time_plus,time_to_add
! AC 14102022 : considering tcc instead of lcc, mcc, hcc '
!INTEGER,SAVE :: record, stepID, tpID,lccID, mccID, hccID, uID, vID, tID,mslID,timeID,dID
!REAL,DIMENSION(2),SAVE :: tp_help,lcc_help, mcc_help, hcc_help, u10_help, v10_help, t2m_help,msl_help,d2m_help
INTEGER,SAVE :: record, stepID, tpID,tccID, uID, vID, tID,mslID,timeID,dID
REAL,DIMENSION(2),SAVE :: tp_help,tcc_help, u10_help, v10_help, t2m_help,msl_help,d2m_help
! ! ! !

INTEGER, DIMENSION(7), SAVE :: itime_ref, intdate
CHARACTER (LEN=1) :: pod
CHARACTER (LEN=4) :: cnextyear
CHARACTER (LEN=2) :: cnextmonth
CHARACTER (LEN=10) :: nameout
REAL, SAVE :: dt

INTERFACE
   FUNCTION scale_and_offset(file_id, var_id)
   !look in the netcdf file if there is a scale factor(1) and a offset(2)
      USE netcdf

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: file_id
      INTEGER, INTENT(IN) :: var_id
      REAL, DIMENSION(2) :: scale_and_offset
      INTEGER :: status

   END FUNCTION
END INTERFACE


procname(pglev+1) = 'usrdef_surface_data'
CALL log_timer_in()
!
!1. Iniitalise on first call
!-------------------------------
!

IF (modfiles(iddesc,ifil,1)%iostat.EQ.0) THEN

!  ---open data file   
   modfiles(iddesc,ifil,1)%iostat = 1
   status = nf90_open(modfiles(iddesc,ifil,1)%filename, nf90_NoWrite, iunit)

   !  ---allocate model data array
   IF (ALLOCATED(meteodat)) DEALLOCATE (meteodat)
   ALLOCATE (meteodat(n1dat,n2dat,nodat),STAT=errstat)

   CALL error_alloc('meteodat',3,(/n1dat,n2dat,nodat/),kndrtype)
   
   IF (ALLOCATED(help_dat)) DEALLOCATE (help_dat)
   ALLOCATE (help_dat(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('help_dat',2,(/n1dat,n2dat/),kndrtype)
   help_dat = 0.0

!  ---allocate readable data 
   IF (ALLOCATED(precip)) DEALLOCATE (precip)
   ALLOCATE (precip(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('precip',2,(/n1dat,n2dat/),kndrtype)

   IF (ALLOCATED(dewt)) DEALLOCATE (dewt)
   ALLOCATE (dewt(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('dewt',2,(/n1dat,n2dat/),kndrtype)

   IF (ALLOCATED(clouddat)) DEALLOCATE (clouddat)
   ! AC 14102022 : considering tcc instead of lcc, mcc, hcc '
   ! ALLOCATE (clouddat(n1dat,n2dat,3),STAT=errstat)
   ! CALL error_alloc('clouddat',3,(/n1dat,n2dat,3/),kndrtype)
   ALLOCATE (clouddat(n1dat,n2dat),STAT=errstat)
   CALL error_alloc('clouddat',2,(/n1dat,n2dat/),kndrtype)
   ! ! ! ! ! ! ! ! ! ! ! ! !
   
   meteodat(:,:,1) = 0.0
   meteodat(:,:,2) = 0.0
   meteodat(:,:,3) = 101300.0
   meteodat(:,:,4) = 285.0
   dewt = 0.005
   ! AC 14102022 : considering tcc instead of lcc, mcc, hcc '
   !clouddat(:,:,:) = 0.5
   clouddat(:,:) = 0.5
   ! ! ! ! ! ! ! ! ! ! ! ! ! !
   ndims=0; nvars=0; numvars=0; nglobatts=0
   dim_timeID=0; maxrecs=0;iyear=0;imont=0;iday=0

   IF (novars.GE.7) THEN 
      WRITE(*,*) 'ERROR - Use of precipitation not yet included'
      STOP
   ENDIF 
   ! --- determine attributes of the variables
   ! --- determine ID of the different variables in the netcdf file
   status = nf90_inq_varid(iunit, "time", timeID)
   status = nf90_inq_varid(iunit, "tp", tpID)
   tp_help = scale_and_offset(iunit,tpID)
   ! AC 14102022 : considering tcc instead of lcc, mcc, hcc '
   !status = nf90_inq_varid(iunit, "lcc", lccID)
   !lcc_help = scale_and_offset(iunit,lccID)
   !status = nf90_inq_varid(iunit, "mcc", mccID)
   !mcc_help = scale_and_offset(iunit,mccID)
   !status = nf90_inq_varid(iunit, "hcc", hccID)
   !hcc_help = scale_and_offset(iunit,hccID)
   status = nf90_inq_varid(iunit, "tcc", tccID)
   tcc_help = scale_and_offset(iunit,tccID)
   ! ! ! ! ! ! ! ! !
   status = nf90_inq_varid(iunit, "u10", uID)
   u10_help = scale_and_offset(iunit,uID)
   status = nf90_inq_varid(iunit, "v10", vID)
   v10_help = scale_and_offset(iunit,vID)
   status = nf90_inq_varid(iunit, "d2m", dID)
   d2m_help = scale_and_offset(iunit,dID)
   status = nf90_inq_varid(iunit, "t2m", tID)
   t2m_help = scale_and_offset(iunit,tID)
   status = nf90_inq_varid(iunit, "msl", mslID)
   msl_help = scale_and_offset(iunit,mslID)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! AC 12/10/2022 - Only one file for now !
   !time
   !   READ (runtitle(1:4),'(I4)') iyear
   !   READ (runtitle(5:6),'(I2)') imont
   iyear  = 2007
   imont  = 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   iday = 1
   itime_ref = (/iyear,imont,iday,0,0,0,0/)
   dt = 1.
   record = 0
   time_to_add = 0
   time_plus = 0
   status = nf90_inquire_dimension(iunit, timeID, len = maxrecs)
   print *, 'maxrecs ',maxrecs
   GOTO 999
ENDIF
!
!
IF (modfiles(io_metsur,ifil,1)%iostat.EQ.1) THEN
   modfiles(iddesc,ifil,1)%iostat=1
   !
   !2. Read meteo data
   !------------------
   !
   ! ---determine time
   record = record +1
   time_plus = 3600*(record-1)
   IF (record .EQ. (maxrecs+1)) THEN
      IF (imont .EQ. 12) THEN
         imont = 1
         iyear = iyear +1
      ELSE
         imont = imont +1
      ENDIF
      itime_ref = (/iyear,imont,iday,0,0,0,0/)
      time_plus = 0
      record = 1
      status = nf90_close(iunit)
      !!! AC : The file name here shouldn't be given explicitely, I guess? 
      ! status = nf90_open('METEO/era5_'//cnextyear//cnextmonth//'.nc', nf90_NoWrite, iunit)
      status = nf90_open(modfiles(io_metsur,1,1)%filename, nf90_NoWrite, iunit)
      !! 
   ENDIF
   
   CALL add_secs_to_date_int(itime_ref, intdate, time_plus, dt)
   ciodatetime = convert_date(intdate)
   
   ! ---read values data
   !  ---atmospheric pressure
   help_dat = 0.0

   status = nf90_get_var(iunit,mslID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   help_dat(:,:)=help_dat(:,n2dat:1:-1)
   meteodat(:,:,3) = help_dat(:,:)*msl_help(1)+msl_help(2)

   !  ---wind velocities
   help_dat = 0.0
   status = nf90_get_var(iunit,uID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   help_dat(:,:)   = help_dat(:,n2dat:1:-1)
   meteodat(:,:,1) = help_dat(:,:)*u10_help(1)+u10_help(2)
   
   help_dat =0.0
   status = nf90_get_var(iunit,vID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   help_dat(:,:)=help_dat(:,n2dat:1:-1)
   meteodat(:,:,2) = help_dat(:,:)*v10_help(1)+v10_help(2)

   !  ---precipitation
   help_dat = 0.0
   status = nf90_get_var(iunit,tpID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   help_dat(:,:)=help_dat(:,n2dat:1:-1)
   precip = help_dat(:,:)*tp_help(1)+tp_help(2)

   !   ---cloud cover
   ! AC 14102022 : considering tcc instead of lcc, mcc, hcc '
   !help_dat = 0.0
   !status = nf90_get_var(iunit,lccID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   !help_dat(:,:)=help_dat(:,n2dat:1:-1)
   !clouddat(:,:,1) = help_dat(:,:)*lcc_help(1)+lcc_help(2)
   
   !help_dat = 0.0
   !status = nf90_get_var(iunit,mccID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   !help_dat(:,:)=help_dat(:,n2dat:1:-1)
   !clouddat(:,:,2) = help_dat(:,:)*mcc_help(1)+mcc_help(2)
   
   !help_dat = 0.0
   !status = nf90_get_var(iunit,hccID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   !help_dat(:,:)=help_dat(:,n2dat:1:-1)
   !clouddat(:,:,3) = help_dat(:,:)*hcc_help(1)+hcc_help(2)

   help_dat = 0.0
   status = nf90_get_var(iunit,tccID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   help_dat(:,:)=help_dat(:,n2dat:1:-1)
   clouddat(:,:) = help_dat(:,:)*tcc_help(1)+tcc_help(2)
   
   !   --- air temperature
   help_dat = 0.0
   status = nf90_get_var(iunit,tID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   help_dat(:,:)=help_dat(:,n2dat:1:-1)
   meteodat(:,:,4) = help_dat(:,:)*t2m_help(1)+t2m_help(2)

   !   --- dewpoint temperature
   help_dat = 0.0
   status = nf90_get_var(iunit,dID,help_dat(:,:),start=(/1,1,record/),count=(/n1dat,n2dat,1/))
   help_dat(:,:)=help_dat(:,n2dat:1:-1)
   dewt(:,:) = help_dat(:,:)*d2m_help(1)+d2m_help(2)

    !3. Convert
    !----------
    !
    
    IF (iopt_temp.EQ.2) THEN
    
    !  ---temperature
       meteodat(:,:,4) = meteodat(:,:,4) - 273.15
    
    !  ---dewpoint to relative humidity
    !Alduchov, O. A., and R. E. Eskridge, 1996: Improved Magnus' form approximation of saturation vapor pressure. J. Appl. Meteor., 35, 601–609.
    !August, E. F., 1828: Ueber die Berechnung der Expansivkraft des Wasserdunstes. Ann. Phys. Chem., 13, 122–137.
    !Magnus, G., 1844: Versuche über die Spannkräfte des Wasserdampfs. Ann. Phys. Chem., 61, 225–247. 
       i_310: DO i=1,n1dat
       j_310: DO j=1,n2dat
          !convert to celcius first
          dewt(i,j) = dewt(i,j)- 273.15
          meteodat(i,j,5) = (exp((17.625 * dewt(i,j))/(243.04 + dewt(i,j)))/exp((17.625 * meteodat(i,j,4))/ &
               & (243.04 + meteodat(i,j,4))))
       ENDDO j_310
       ENDDO i_310
    !  ---cloud cover
    
       i_320: DO i=1,n1dat
          j_320: DO j=1,n2dat
!             AC 09.12.2022 : This was commented out to consider total cloud cover instead of low medium high components
!          IF (clouddat(i,j,2).NE.0.0) THEN
!             meteodat(i,j,6) = MAXVAL(clouddat(i,j,:))
!          ELSE
!             meteodat(i,j,6) = 1.0 - (1.0-clouddat(i,j,1))*(1.0-clouddat(i,j,3))
!          ENDIF
          meteodat(i,j,6) = clouddat(i,j)
          meteodat(i,j,6) = MAX(0.0,MIN(meteodat(i,j,6),1.0))
       ENDDO j_320
       ENDDO i_320
    ENDIF
    
    !
    !4. Store
    !--------
    surdata = meteodat(:,:,1:novars)
ENDIF 
999 CALL log_timer_out()
RETURN

END SUBROUTINE usrdef_surface_data

