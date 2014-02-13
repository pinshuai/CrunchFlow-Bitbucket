!************** (C) COPYRIGHT 1995 Carl I. Steefel *******************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:01:43
 
!                      All Rights Reserved

!  OSRT (Operator Splitting Reactive Transport) IS PROVIDED "AS IS"
!  AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED. THE USER ASSUMES ALL RISKS
!  OF USING OSRT. THERE IS NO CLAIM OF THE MERCHANTABILITY OR FITNESS
!  FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO USERS AT
!  ANY SITES OTHER THAN YOUR OWN.
!**********************************************************************

SUBROUTINE newfile(fn,suf,fnv,nint,i)
USE CrunchType
USE CrunchFunctions
USE params
USE strings

IMPLICIT NONE

CHARACTER (LEN=mls), INTENT(IN)                    :: fn
CHARACTER (LEN=mls), INTENT(IN)                    :: suf
CHARACTER (LEN=mls), INTENT(OUT)                   :: fnv
INTEGER(I4B), INTENT(IN)                           :: nint
INTEGER(I4B), INTENT(IN)                           :: i

INTEGER(I4B)                                       :: n
INTEGER(I4B)                                       :: in
INTEGER(I4B)                                       :: ls
INTEGER(I4B)                                       :: lfileno
INTEGER(I4B)                                       :: lsuffix
INTEGER(I4B)                                       :: lfnv
INTEGER(I4B)                                       :: ncheck

CHARACTER (LEN=mls)                                :: dumstring
CHARACTER (LEN=mls)                                :: dummy


CALL stringlen(fn,ls)

dumstring = IntegerToCharacter(nint)

CALL stringlen(dumstring,lfileno)

CALL stringlen(suf,lsuffix)

fnv = ' '

fnv(1:ls) = fn(1:ls)
lfnv = ls + lfileno
fnv(ls+1:ls+lfnv) = dumstring(1:lfileno)
fnv(lfnv+1:lfnv+lsuffix) = suf(1:lsuffix)


RETURN
END SUBROUTINE newfile
!**********************************************************

