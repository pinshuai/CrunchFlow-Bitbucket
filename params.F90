MODULE params

  USE crunchtype

  INTEGER(I4B), PARAMETER :: mcomp=50
  INTEGER(I4B), PARAMETER :: mrct=500
  INTEGER(I4B), PARAMETER :: mspec=1000
  INTEGER(I4B), PARAMETER :: mgas=75
  INTEGER(I4B), PARAMETER :: msurf=25
  INTEGER(I4B), PARAMETER :: mtot=mcomp+mspec
  INTEGER(I4B), PARAMETER :: mexch=15
  INTEGER(I4B), PARAMETER :: meqn=mcomp+mexch+msurf
  INTEGER(I4B), PARAMETER :: mexch_sec=150
  INTEGER(I4B), PARAMETER :: msurf_sec=250
  INTEGER(I4B), PARAMETER :: msurftot=msurf+msurf_sec
  INTEGER(I4B), PARAMETER :: mx=1000
  INTEGER(I4B), PARAMETER :: mp=meqn*mx*meqn*3
  INTEGER(I4B), PARAMETER :: my=1000
  INTEGER(I4B), PARAMETER :: mz=1
  INTEGER(I4B), PARAMETER :: mwidth=2
  INTEGER(I4B), PARAMETER :: mxyz=mx*my*mz
  INTEGER(I4B), PARAMETER :: nc=50
  INTEGER(I4B), PARAMETER :: meqneq=nc+mexch+msurf
  INTEGER(I4B), PARAMETER :: nm=mrct
  INTEGER(I4B), PARAMETER :: mcmplx=mspec
  INTEGER(I4B), PARAMETER :: ng=mgas
  INTEGER(I4B), PARAMETER :: msec=mcmplx+nm+ng+msurf_sec
  INTEGER(I4B), PARAMETER :: ndim=nc+mcmplx+nm+ng+msurf_sec
  INTEGER(I4B), PARAMETER :: nleq=4  
  INTEGER(I4B), PARAMETER :: ntmp=8 
  INTEGER(I4B), PARAMETER :: mbasis=1
  INTEGER(I4B), PARAMETER :: mls=132 
  INTEGER(I4B), PARAMETER :: mchem=5000
  INTEGER(I4B), PARAMETER :: mhet=100
  INTEGER(I4B), PARAMETER :: mperm=100
  INTEGER(I4B), PARAMETER :: mzone=60
  INTEGER(I4B), PARAMETER :: mreact=10
  INTEGER(I4B), PARAMETER :: mpre=35
  INTEGER(I4B), PARAMETER :: maqkin=25

  REAL(DP), PARAMETER :: rgasKCAL=0.001987d0      !!  kcal/deg-mol
  REAL(DP), PARAMETER :: rgas=.00831470       !!  kJ/deg-mol
  REAL(DP), PARAMETER :: tk25=1.0_dp/298.15_dp 
  REAL(DP), PARAMETER :: clg=2.30258509299405
  REAL(DP), PARAMETER :: secyr=365.0d0*24.0d0*60.0d0*60.0d0
  REAL(DP), PARAMETER :: ed=5.0
  REAL(DP), PARAMETER :: d0=1.e-09 
  REAL(DP), PARAMETER :: t_half = -0.6931471805599453

END MODULE params

