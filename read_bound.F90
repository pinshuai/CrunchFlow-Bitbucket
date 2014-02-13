!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:25
 
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

SUBROUTINE read_bound(nout,nchem,nx,ny,nz,ncomp,nspec,ngas,nkin,  &
    nexchange,nexch_sec,nsurf,nsurf_sec)
USE crunchtype
USE params
USE concentration
USE mineral
USE strings
USE medium
USE temperature

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nchem
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nspec
INTEGER(I4B), INTENT(IN)                                    :: ngas
INTEGER(I4B), INTENT(IN)                                    :: nkin
INTEGER(I4B), INTENT(IN)                                    :: nexchange
INTEGER(I4B), INTENT(IN)                                    :: nexch_sec
INTEGER(I4B), INTENT(IN)                                    :: nsurf
INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: kk
INTEGER(I4B)                                                :: k
INTEGER(I4B)                                                :: nex
INTEGER(I4B)                                                :: ns


jc = -1

REWIND nout

10  READ(nout,'(a)',END=800) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
  IF (ssch == 'x_begin') THEN
    WRITE(*,*)
    WRITE(*,*) ' X = 0 boundary condition found'
    WRITE(*,*)
    IF (nx == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' NX = 1, so boundary condition ignored'
      WRITE(*,*)
!             pause
      GO TO 10
    END IF
!   Look for condition label following the boundary
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
!     Check the trailing string to see if it matches one of the
!       geochemical conditions given
      DO nco = 1,nchem
        IF (ssch == condlabel(nco)) THEN
          DO ik = 1,ncomp+nspec
            spb(ik,1) = spcond10(ik,nco)
          END DO
          DO kk = 1,ngas
            spbgas(kk,1) = spcondgas10(kk,nco)
          END DO
          DO k = 1,nkin
            volb(k,1) = volin(k,nco)
          END DO
          DO nex = 1,nexchange+nexch_sec
            spexb(nex,1) = spcondex10(nex,nco)*AqueousToBulkCond(nco)
          END DO
          DO ns = 1,nsurf+nsurf_sec
            spsurfb(ns,1) = spcondsurf10(ns,nco)*AqueousToBulkCond(nco)
          END DO
          GO TO 200
        END IF
      END DO
      WRITE(*,*)
      WRITE(*,*) ' Condition label for BC at X_BEGIN not found'
      WRITE(*,*)
      READ(*,*)
      STOP
      200         CONTINUE
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Blank string following X_BEGIN'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
!   Look for the type of boundary condition (Dirichlet, flux, etc.)
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'dirichlet') THEN
        jc(1) = 1
      ELSE IF (ssch == 'first') THEN
        jc(1) = 1
      ELSE IF (ssch == 'flux') THEN
        jc(1) = 2
      ELSE IF (ssch == 'second') THEN
        jc(1) = 2
      ELSE IF (ssch == 'neumann') THEN
        jc(1) = 2
      ELSE IF (ssch == 'third') THEN
        jc(1) = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont recognize this type of boundary condition'
        WRITE(*,*) ' At boundary X_BEGIN'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      jc(1) = 2     ! Assume no flux boundary
    END IF
    GO TO 10
  ELSE IF (ssch == 'x_end') THEN
    WRITE(*,*)
    WRITE(*,*) ' X = NX boundary condition found'
    WRITE(*,*)
    IF (nx == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' NX = 1, so boundary condition ignored'
      WRITE(*,*)
      GO TO 10
    END IF
!   Look for condition label following the boundary
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
!     Check the trailing string to see if it matches one of the
!       geochemical conditions given
      DO nco = 1,nchem
        IF (ssch == condlabel(nco)) THEN
          DO ik = 1,ncomp+nspec
            spb(ik,2) = spcond10(ik,nco)
          END DO
          DO kk = 1,ngas
            spbgas(kk,2) = spcondgas10(kk,nco)
          END DO
          DO k = 1,nkin
            volb(k,2) = volin(k,nco)
          END DO
          DO nex = 1,nexchange+nexch_sec
            spexb(nex,2) = spcondex10(nex,nco)*AqueousToBulkCond(nco)
          END DO
          DO ns = 1,nsurf+nsurf_sec
            spsurfb(ns,2) = spcondsurf10(ns,nco)*AqueousToBulkCond(nco)
          END DO
          GO TO 300
        END IF
      END DO
      WRITE(*,*)
      WRITE(*,*) ' Condition label for BC at X_END not found'
      WRITE(*,*)
      READ(*,*)
      STOP
      300         CONTINUE
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Blank string following X_END'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
!   Look for the type of boundary condition (Dirichlet, flux, etc.)
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'dirichlet') THEN
        jc(2) = 1
      ELSE IF (ssch == 'first') THEN
        jc(2) = 1
      ELSE IF (ssch == 'flux') THEN
        jc(2) = 2
      ELSE IF (ssch == 'second') THEN
        jc(2) = 2
      ELSE IF (ssch == 'neumann') THEN
        jc(2) = 2
      ELSE IF (ssch == 'third') THEN
        jc(2) = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont recognize this type of boundary condition'
        WRITE(*,*) ' At boundary X_END'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      jc(2) = 2     ! Assume no flux boundary
    END IF
    GO TO 10
  ELSE IF (ssch == 'y_begin') THEN
    WRITE(*,*)
    WRITE(*,*) ' Y = 0 boundary condition found'
    WRITE(*,*)
    IF (ny == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' NY = 1, so boundary condition ignored'
      WRITE(*,*)
      GO TO 10
    END IF
!   Look for condition label following the boundary
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
!     Check the trailing string to see if it matches one of the
!       geochemical conditions given
      DO nco = 1,nchem
        IF (ssch == condlabel(nco)) THEN
          DO ik = 1,ncomp+nspec
            spb(ik,3) = spcond10(ik,nco)
          END DO
          DO kk = 1,ngas
            spbgas(kk,3) = spcondgas10(kk,nco)
          END DO
          DO k = 1,nkin
            volb(k,3) = volin(k,nco)
          END DO
          DO nex = 1,nexchange+nexch_sec
            spexb(nex,3) = spcondex10(nex,nco)*AqueousToBulkCond(nco)
          END DO
          DO ns = 1,nsurf+nsurf_sec
            spsurfb(ns,3) = spcondsurf10(ns,nco)*AqueousToBulkCond(nco)
          END DO
          GO TO 400
        END IF
      END DO
      WRITE(*,*)
      WRITE(*,*) ' Condition label for BC at Y_BEGIN not found'
      WRITE(*,*)
      READ(*,*)
      STOP
      400         CONTINUE
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Blank string following Y_BEGIN'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
!   Look for the type of boundary condition (Dirichlet, flux, etc.)
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'dirichlet') THEN
        jc(3) = 1
      ELSE IF (ssch == 'first') THEN
        jc(3) = 1
      ELSE IF (ssch == 'flux') THEN
        jc(3) = 2
      ELSE IF (ssch == 'second') THEN
        jc(3) = 2
      ELSE IF (ssch == 'neumann') THEN
        jc(3) = 2
      ELSE IF (ssch == 'third') THEN
        jc(3) = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont recognize this type of boundary condition'
        WRITE(*,*) ' At boundary Y_BEGIN'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      jc(3) = 2     ! Assume no flux boundary
    END IF
    GO TO 10
  ELSE IF (ssch == 'y_end') THEN
    WRITE(*,*)
    WRITE(*,*) ' Y = NY boundary condition found'
    WRITE(*,*)
    IF (ny == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' NY = 1, so boundary condition ignored'
      WRITE(*,*)
      GO TO 10
    END IF
!   Look for condition label following the boundary
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
!     Check the trailing string to see if it matches one of the
!       geochemical conditions given
      DO nco = 1,nchem
        IF (ssch == condlabel(nco)) THEN
          DO ik = 1,ncomp+nspec
            spb(ik,4) = spcond10(ik,nco)
          END DO
          DO kk = 1,ngas
            spbgas(kk,4) = spcondgas10(kk,nco)
          END DO
          DO k = 1,nkin
            volb(k,4) = volin(k,nco)
          END DO
          DO nex = 1,nexchange+nexch_sec
            spexb(nex,4) = spcondex10(nex,nco)*AqueousToBulkCond(nco)
          END DO
          DO ns = 1,nsurf+nsurf_sec
            spsurfb(ns,4) = spcondsurf10(ns,nco)*AqueousToBulkCond(nco)
          END DO
          GO TO 500
        END IF
      END DO
      WRITE(*,*)
      WRITE(*,*) ' Condition label for BC at Y_END not found'
      WRITE(*,*)
      READ(*,*)
      STOP
      500         CONTINUE
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Blank string following Y_END'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
!   Look for the type of boundary condition (Dirichlet, flux, etc.)
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'dirichlet') THEN
        jc(4) = 1
      ELSE IF (ssch == 'first') THEN
        jc(4) = 1
      ELSE IF (ssch == 'flux') THEN
        jc(4) = 2
      ELSE IF (ssch == 'second') THEN
        jc(4) = 2
      ELSE IF (ssch == 'neumann') THEN
        jc(4) = 2
      ELSE IF (ssch == 'third') THEN
        jc(4) = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont recognize this type of boundary condition'
        WRITE(*,*) ' At boundary Y_END'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      jc(4) = 2     ! Assume no flux boundary
    END IF
    GO TO 10
  ELSE IF (ssch == 'z_begin') THEN
    WRITE(*,*)
    WRITE(*,*) ' Z = 0 boundary condition found'
    WRITE(*,*)
    IF (nz == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' NZ = 1, so boundary condition ignored'
      WRITE(*,*)
      GO TO 10
    END IF
!   Look for condition label following the boundary
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
!     Check the trailing string to see if it matches one of the
!       geochemical conditions given
      DO nco = 1,nchem
        IF (ssch == condlabel(nco)) THEN
          DO ik = 1,ncomp+nspec
            spb(ik,5) = spcond10(ik,nco)
          END DO
          DO kk = 1,ngas
            spbgas(kk,5) = spcondgas10(kk,nco)
          END DO
          DO k = 1,nkin
            volb(k,5) = volin(k,nco)
          END DO
          DO nex = 1,nexchange+nexch_sec
            spexb(nex,5) = spcondex10(nex,nco)*AqueousToBulkCond(nco)
          END DO
          DO ns = 1,nsurf+nsurf_sec
            spsurfb(ns,5) = spcondsurf10(ns,nco)*AqueousToBulkCond(nco)
          END DO
          GO TO 600
        END IF
      END DO
      WRITE(*,*)
      WRITE(*,*) ' Condition label for BC at Z_BEGIN not found'
      WRITE(*,*)
      READ(*,*)
      STOP
      600         CONTINUE
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Blank string following Z_BEGIN'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
!   Look for the type of boundary condition (Dirichlet, flux, etc.)
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'dirichlet') THEN
        jc(5) = 1
      ELSE IF (ssch == 'first') THEN
        jc(5) = 1
      ELSE IF (ssch == 'flux') THEN
        jc(5) = 2
      ELSE IF (ssch == 'second') THEN
        jc(5) = 2
      ELSE IF (ssch == 'neumann') THEN
        jc(5) = 2
      ELSE IF (ssch == 'third') THEN
        jc(5) = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont recognize this type of boundary condition'
        WRITE(*,*) ' At boundary Z_BEGIN'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      jc(5) = 2     ! Assume no flux boundary
    END IF
    GO TO 10
  ELSE IF (ssch == 'z_end') THEN
    WRITE(*,*)
    WRITE(*,*) ' Z = NZ boundary condition found'
    WRITE(*,*)
    IF (nz == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' NZ = 1, so boundary condition ignored'
      WRITE(*,*)
      GO TO 10
    END IF
!   Look for condition label following the boundary
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
!     Check the trailing string to see if it matches one of the
!       geochemical conditions given
      DO nco = 1,nchem
        IF (ssch == condlabel(nco)) THEN
          DO ik = 1,ncomp+nspec
            spb(ik,6) = spcond10(ik,nco)
          END DO
          DO kk = 1,ngas
            spbgas(kk,6) = spcondgas10(kk,nco)
          END DO
          DO k = 1,nkin
            volb(k,6) = volin(k,nco)
          END DO
          DO nex = 1,nexchange+nexch_sec
            spexb(nex,6) = spcondex10(nex,nco)*AqueousToBulkCond(nco)
          END DO
          DO ns = 1,nsurf+nsurf_sec
            spsurfb(ns,6) = spcondsurf10(ns,nco)*AqueousToBulkCond(nco)
          END DO
          GO TO 700
        END IF
      END DO
      WRITE(*,*)
      WRITE(*,*) ' Condition label for BC at Z_END not found'
      WRITE(*,*)
      READ(*,*)
      STOP
      700         CONTINUE
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Blank string following Z_END'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
!   Look for the type of boundary condition (Dirichlet, flux, etc.)
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'dirichlet') THEN
        jc(6) = 1
      ELSE IF (ssch == 'first') THEN
        jc(6) = 1
      ELSE IF (ssch == 'flux') THEN
        jc(6) = 2
      ELSE IF (ssch == 'second') THEN
        jc(6) = 2
      ELSE IF (ssch == 'neumann') THEN
        jc(6) = 2
      ELSE IF (ssch == 'third') THEN
        jc(6) = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont recognize this type of boundary condition'
        WRITE(*,*) ' At boundary Z_END'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      jc(6) = 2     ! Assume no flux boundary
    END IF
    GO TO 10
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Cannot interpret string in boundary read'
    WRITE(*,*) ssch
    WRITE(*,*) ' Ignoring this line'
    WRITE(*,*)
    GO TO 10
  END IF
ELSE              !  Blank string, read next line
  GO TO 10
END IF

800 CONTINUE
IF (nx > 1 .AND. jc(1) == -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' No boundary conditions found for JX = 1 '
  WRITE(*,*) 
  STOP
END IF
IF (nx > 1 .AND. jc(2) == -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' No boundary conditions found for JX = NX '
  WRITE(*,*) 
  STOP
END IF
IF (ny > 1 .AND. jc(3) == -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' No boundary conditions found for JY = 1 '
  WRITE(*,*) 
  STOP
END IF
IF (ny > 1 .AND. jc(4) == -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' No boundary conditions found for JY = NY '
  WRITE(*,*) 
  STOP
END IF
IF (nz > 1 .AND. jc(5) == -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' No boundary conditions found for JZ = 1 '
  WRITE(*,*) 
  STOP
END IF
IF (nz > 1 .AND. jc(6) == -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' No boundary conditions found for JZ = NZ '
  WRITE(*,*) 
  STOP
END IF
       
RETURN
END SUBROUTINE read_bound
