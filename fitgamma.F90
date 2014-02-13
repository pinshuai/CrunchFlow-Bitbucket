!******************        GIMRT98      ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 09:56:41
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!********************** C.I. Steefel *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE fitgamma(nbasis,ntemp,adh,bvec,vec)
USE crunchtype
USE params

IMPLICIT NONE

interface
SUBROUTINE ludcmp90(a,indx,d,n)
USE crunchtype
REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
INTEGER(I4B), INTENT(IN)                                   :: n
REAL(DP), intent(out)                                      :: d
END SUBROUTINE ludcmp90
END interface

interface
SUBROUTINE lubksb90(a,indx,b,n)
USE crunchtype
REAL(DP), DIMENSION(:,:), intent(in)                       :: a
INTEGER(I4B), DIMENSION(:), intent(in)                     :: indx
REAL(DP), DIMENSION(:), intent(inout)                      :: b
INTEGER(I4B), INTENT(IN)                                   :: n
END SUBROUTINE lubksb90
END interface

!  External variables needed for dimensions

INTEGER(I4B)                                               :: nbasis
INTEGER(I4B)                                               :: ntemp

!  External variables and arrays

REAL(DP), DIMENSION(ntemp), INTENT(IN)                     :: adh
REAL(DP), DIMENSION(nbasis), INTENT(OUT)                   :: bvec
REAL(DP), DIMENSION(nbasis,ntemp), INTENT(IN)              :: vec

!  Internal variables and arrays

REAL(DP), DIMENSION(nbasis,nbasis)                         :: w
REAL(DP)                                                   :: dd

INTEGER(I4B), DIMENSION(nbasis)                            :: indx
INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: k
INTEGER(I4B)                                               :: j
CHARACTER (LEN=1)                                          :: trans
INTEGER(I4B)                                               :: info
INTEGER(I4B), PARAMETER                                    :: ione=1

trans = 'N'

DO j = 1, nbasis
  bvec(j) = 0.d0
  DO i = 1, ntemp
    bvec(j) = bvec(j) + adh(i)*vec(j,i)
  END DO
END DO

DO j = 1, nbasis
  DO k = j, nbasis
    w(j,k) = 0.d0
    DO i = 1, ntemp
      w(j,k) = w(j,k) + vec(j,i)*vec(k,i)
    END DO
    IF (j /= k) w(k,j) = w(j,k)
  END DO
END DO

!!CALL ludcmp90(w,indx,dd,nbasis)
!!CALL lubksb90(w,indx,bvec,nbasis)


!!  call ludcmp90(fj,indx,det,neqn)
!!  call lubksb90(fj,indx,beta,neqn)

CALL dgetrf(nbasis,nbasis,w,nbasis,indx,info)
CALL dgetrs(trans,nbasis,ione,w,nbasis,indx,bvec,nbasis,info)

!!  CALL dgetrf(neqn,neqn,fj,neqn,indx,info)
!!  CALL dgetrs(trans,neqn,ione,fj,neqn,indx,beta,neqn,info)

RETURN
END SUBROUTINE fitgamma
!**************************************************************
