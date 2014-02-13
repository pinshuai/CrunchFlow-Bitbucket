!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:44
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE read_diffusion(nout,nx,ny,nz)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE transport
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1

nxyz = nx*ny*nz

!  Default values here, overridden below

idiffus = 0
dzero = 0.0d0
dcoeff = 0.0d0
activation = 5.0d0
formation = 1.0d0
uli = 1.0d0

REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check "fix_diffusion" string which fixes diffusion coefficient (ignoring temperature)
  
  IF (ssch == 'fix_diffusion') THEN    ! Overrides other specifications of diffusion coefficient
    idiffus = 1
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      dcoeff = 0.0
      dzero = 0.0
      RETURN
    END IF
    
!  Looking for diffusion coefficient
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        dcoeff = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "fix_diffusion" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "fix_diffusion" '
      WRITE(*,*) ' Assuming diffusion coefficient = 0'
      WRITE(*,*)
      dcoeff = 0.0
    END IF
    
  ELSE
    GO TO 10
  END IF
  GO TO 10
  
  
ELSE         ! No string found
  GO TO 10
END IF

IF (idiffus == 1) GO TO 3000

1000 REWIND  nout

20  READ(nout,'(a)',END=2000) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check "calculate_diffusion" string which gives diffusion coefficient at 25C
  
  IF (ssch == 'calculate_diffusion') THEN
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      dzero = 0.0
      RETURN
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        dzero = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "calculate_diffusion" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "calculate_diffusion" '
      WRITE(*,*) ' Assuming diffusion coefficient = 0 '
      WRITE(*,*)
      dzero = 0
    END IF
    
  ELSE
    GO TO 20
  END IF
  GO TO 20
  
  
ELSE         ! No string found
  GO TO 20
END IF

2000 REWIND nout

30  READ(nout,'(a)',END=3000) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check for "diffusion_activation" string which gives activation energy (kcal) for diffusion coefficient
  
  IF (ssch == 'diffusion_activation') THEN
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      activation = 0.0
      RETURN
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        activation = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "diffusion_activation" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "diffusion_activation" '
      WRITE(*,*) ' Assuming canonical diffusion coefficient = 5 kcal '
      WRITE(*,*)
      activation = 5.0
    END IF
    
  ELSE
    GO TO 30
  END IF
  GO TO 30
  
ELSE         ! No string found
  GO TO 30
END IF

3000 REWIND nout

40  READ(nout,'(a)',END=4000) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check for "formation_factor"
  
  IF (ssch == 'formation_factor') THEN
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      formation = 1.0
      RETURN
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        formation = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "formation_factor" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "formation_factor" '
      WRITE(*,*) ' Assuming value = 1  '
      WRITE(*,*)
      formation = 1.0
    END IF
    
  ELSE
    GO TO 40
  END IF
  GO TO 40
  
ELSE         ! No string found
  GO TO 40
END IF

4000 REWIND nout

50  READ(nout,'(a)',END=5000) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check for "cementation_exponent"
  
  IF (ssch == 'cementation_exponent') THEN
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      uli = 1.0d0
      RETURN
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        uli = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "cementation_exponent" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "cementation_exponent" '
      WRITE(*,*) ' Assuming value = 1  '
      WRITE(*,*)
      uli = 1.0
    END IF
    
  ELSE
    GO TO 50
  END IF
  GO TO 50
  
ELSE         ! No string found
  GO TO 50
END IF


5000 RETURN
END SUBROUTINE read_diffusion
