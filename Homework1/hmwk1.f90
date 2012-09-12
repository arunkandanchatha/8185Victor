MODULE nrtype
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    INTEGER, PARAMETER :: SP = KIND(1.0D0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
    INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
    INTEGER, PARAMETER :: LGT = KIND(.true.)
        INTEGER, PARAMETER :: wp=DP
    REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
    REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
    REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
    REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
    REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
    REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
    REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
    REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp
    TYPE sprs2_sp
        INTEGER(I4B) :: n,len
        REAL(SP), DIMENSION(:), POINTER :: val
        INTEGER(I4B), DIMENSION(:), POINTER :: irow
        INTEGER(I4B), DIMENSION(:), POINTER :: jcol
    END TYPE sprs2_sp
    TYPE sprs2_dp
        INTEGER(I4B) :: n,len
        REAL(DP), DIMENSION(:), POINTER :: val
        INTEGER(I4B), DIMENSION(:), POINTER :: irow
        INTEGER(I4B), DIMENSION(:), POINTER :: jcol
    END TYPE sprs2_dp
END MODULE nrtype

MODULE  hmwks
    use nrtype
    implicit none

    contains

    subroutine mystop()
        STOP 0
    end subroutine mystop

    subroutine H2()
    ! This subroutine performs the linear interpoliation of exp

        INTEGER(I4B) :: j, calcInt
        REAL(DP), DIMENSION (1:11) :: initialValues=(/(exp(j/10.),j=0,10,1)/)
        REAL(DP) :: x, y, calcVal

        print *, "HOMEWORK 2"
        write (*,*) "For what value of 'x' (between 0 and 1) do you want e^x?"
        flush(6)
        read(*,*) x

        !error checking
        IF (x<0 .OR. x>1) THEN
            print *, "Invalid x"
            call mystop()
        END IF

        !find proper interval using linear interpolation
        calcInt = INT(x*10)+1
        y = (x*10)-INT(x*10)

        calcVal = initialvalues(calcInt) + y*(initialValues(calcInt+1) - initialValues(calcInt))
        print *, "e^",x,"=",calcval
    end subroutine H2

    subroutine H3()
    ! This subroutine uses bisection to calculate x such that sin(2x)-2x=0

        INTEGER(I4B) :: counter=1, MAXCOUNT=10000
        REAL(DP) :: lower, upper, mid, sens
        REAL(DP) :: eps=100

        print *, "HOMEWORK 3: Solving sin(2*Pi*x)-2x=0"
        write(*,*) "Starting Lower Bound: "
        flush(6)
        read (*,*) lower
        write(*,*) "Starting Upper Bound: "
        flush(6)
        read (*,*) upper
        write(*,*) "Allowable error (in decimal form): "
        flush(6)
        read (*,*) sens

        DO WHILE ((ABS(eps) >= sens) .AND. (counter < MAXCOUNT))
            mid = (upper+lower)/2
            eps = sin(2.0*PI_D*mid)-2*mid

            IF (eps .GT. 0.0D0) THEN
                lower = mid
            ELSE
                upper = mid
            END IF

            counter = counter + 1
        ENDDO

        IF (counter >= MAXCOUNT) THEN
            print *, "Could not find solution in ",counter," iterations."
            call mystop()
        ENDIF

        print *, "Solution: ",mid
    end subroutine H3

end module hmwks

program main
    use hmwks
    implicit none
    call H2()
    call H3()
end program main

