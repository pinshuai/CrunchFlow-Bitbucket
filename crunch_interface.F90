MODULE crunch_interface

INTERFACE
  SUBROUTINE CourantStepAlt(nx,ny,nz,dtmaxcour)
  USE crunchtype
  IMPLICIT NONE    
  INTEGER(I4B), INTENT(IN)                                 :: nx
  INTEGER(I4B), INTENT(IN)                                 :: ny
  INTEGER(I4B), INTENT(IN)                                 :: nz
  REAL(DP), INTENT(OUT)                                    :: dtmaxcour
  END SUBROUTINE CourantStepAlt
END INTERFACE

INTERFACE
  SUBROUTINE ludcmp90(a,indx,d,n)
  USE crunchtype
  REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
  INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
  INTEGER(I4B), INTENT(IN)                                   :: n
  REAL(DP), intent(out)                                      :: d
  END SUBROUTINE ludcmp90
END INTERFACE

INTERFACE
  SUBROUTINE lubksb90(a,indx,b,n)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP),  DIMENSION(:,:), INTENT(IN)                          :: a
  REAL(DP),  DIMENSION(:), INTENT(IN OUT)                        :: b
  INTEGER(I4B),  DIMENSION(:),INTENT(IN)                         :: indx
  INTEGER(I4B), INTENT(IN)                                       :: n
  END SUBROUTINE lubksb90
END INTERFACE

INTERFACE
  SUBROUTINE tridag_ser(a,b,c,r,u)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(DP), DIMENSION(:), INTENT(OUT) :: u
  END SUBROUTINE tridag_ser
END INTERFACE

INTERFACE
  FUNCTION locate(xx,x)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: x
  INTEGER(I4B) :: locate
  END FUNCTION locate
END INTERFACE

!!INTERFACE
!!SUBROUTINE assemble(nx,ny,nz,ncomp,nspec,nkin,nrct,ngas,ikin,  &
!!    nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,nn,delt,   &
!!    user,amatpetsc,blockfl,blockl,blockm,blockr,blockfr)
!!  USE crunchtype
!!  IMPLICIT NONE
!!  INTEGER(I4B), INTENT(IN)                      :: nx
!!  INTEGER(I4B), INTENT(IN)                      :: ny
!!  INTEGER(I4B), INTENT(IN)                      :: nz
!!  INTEGER(I4B), INTENT(IN)                      :: ncomp
!!  INTEGER(I4B), INTENT(IN)                      :: nspec
!!  INTEGER(I4B), INTENT(IN)                      :: nkin
!!  INTEGER(I4B), INTENT(IN)                      :: nrct
!!  INTEGER(I4B), INTENT(IN)                      :: ngas
!!  INTEGER(I4B), INTENT(IN)                      :: ikin
!!  INTEGER(I4B), INTENT(IN)                      :: nexchange
!!  INTEGER(I4B), INTENT(IN)                      :: nexch_sec
!!  INTEGER(I4B), INTENT(IN)                      :: nsurf
!!  INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
!!  INTEGER(I4B), INTENT(IN)                      :: npot
!!  INTEGER(I4B), INTENT(IN)                      :: ndecay
!!  INTEGER(I4B), INTENT(IN)                      :: nn
!!  REAL(DP), INTENT(IN)                          :: delt

!!!********************** Add in PETSC declarations for f90 variables *****
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockfl
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockl
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockm
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockr
!!REAL(DP), DIMENSION(ncomp+nexchange+nsurf+npot,ncomp+nexchange+nsurf+npot) :: blockfr
!!!************************* End PETSc declarations ****************************
!!  END SUBROUTINE assemble
!!END INTERFACE

INTERFACE
  SUBROUTINE readModFlowParameters(nout,nchem,nparams,jxTemp,jyTemp,jzTemp,  &
              conditionNum,modflowstring,lenstring)
  USE crunchtype
  USE strings
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: nout
  INTEGER(I4B), INTENT(IN)                                    :: nchem
  INTEGER(I4B), INTENT(OUT)                                   :: nparams
  INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jxTemp
  INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jyTemp
  INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jzTemp
  INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: conditionNum
  CHARACTER (LEN=mls), INTENT(IN)                             :: modflowstring
  INTEGER(I4B), INTENT(IN)                                    :: lenstring
  END SUBROUTINE readModFlowParameters
END INTERFACE

INTERFACE
  FUNCTION imaxloc(arr)
  USE crunchtype
  REAL(DP), DIMENSION(:), INTENT(IN) :: arr
  INTEGER(I4B) :: imaxloc
  END FUNCTION imaxloc
END INTERFACE

INTERFACE
  SUBROUTINE SurfaceConcentration(jx,jy,jz,ncomp,nsurf,nsurf_sec)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: jx
  INTEGER(I4B), INTENT(IN)                                    :: jy
  INTEGER(I4B), INTENT(IN)                                    :: jz
  INTEGER(I4B), INTENT(IN)                                    :: ncomp
  INTEGER(I4B), INTENT(IN)                                    :: nsurf
  INTEGER(I4B), INTENT(IN)                                    :: nsurf_sec
  END SUBROUTINE SurfaceConcentration
END INTERFACE

INTERFACE
  SUBROUTINE ExchangeConcentration(jx,jy,jz,ncomp,nexchange,nexch_sec)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: jx
  INTEGER(I4B), INTENT(IN)                                    :: jy
  INTEGER(I4B), INTENT(IN)                                    :: jz
  INTEGER(I4B), INTENT(IN)                                    :: ncomp
  INTEGER(I4B), INTENT(IN)                                    :: nexchange
  INTEGER(I4B), INTENT(IN)                                    :: nexch_sec
  END SUBROUTINE ExchangeConcentration
END INTERFACE

END MODULE crunch_interface
