module c_lbfgs

use iso_c_binding, only: c_int, c_double, c_bool

implicit none

contains

subroutine c_lbfgs_step(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG) bind(c)

    integer(c_int)  :: N,M,IPRINT(2),IFLAG
    real(c_double)  :: X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
    real(c_double)  :: F,EPS,XTOL
    integer(c_int)  :: DIAGCO
    external LB2
    call LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)

end subroutine

end module
