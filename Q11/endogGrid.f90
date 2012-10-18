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

! This module specifies the model. Two samples are given here:
!     1) Endogenous growth model with no labour
!     2) Endogenous growth model with labour
!
! USAGE: This module is required for the endogenousGrid module defined below.
!        You must specify all subroutines defined below
module modelDefinition
    use nrtype
    implicit none

    !----------------------------------------------------------------
    ! 1. Calibration
    !----------------------------------------------------------------

    real (DP), parameter :: delta  = 0.0196                         ! Depreciation
    real (DP), parameter :: alpha  = 0.4                            ! Capital Share
    real (DP), parameter :: tau    = 2.d0                           ! risk aversion
    real (DP), parameter :: beta   = 0.9896                         ! Discount factor
    real (DP), parameter :: theta   = 0.357                         ! weight on consumption
    real (DP), parameter :: rho    = 0.95                           ! Autoregressive
    real (DP), parameter :: sigma  = 0.007                          ! Variance

    real(DP) :: k_ss, l_ss, c_ss, r_ss, y_ss

    private::alpha, beta, tau, delta, theta, rho, sigma
    private:: k_ss, l_ss, c_ss, r_ss, y_ss

contains

    subroutine sub_modelclean()
    end subroutine sub_modelclean

    subroutine sub_modelstop(calling)
        ! a personal stop subroutine. Makes it easier to edit behaviour of stop. All
        ! functions and subroutines should call this.
        !
        ! INPUTS: calling - a string indicating where this subroutine was called from

        CHARACTER (LEN=*), intent(in) :: calling
        print *, "STOP: ", calling
        call sub_modelclean()
        STOP 0
    end subroutine sub_modelstop

    function sub_modelOutput(inputVars)
        ! a function that returns the output given a set of required inputs. The typical
        ! minimal inputs are:
        !   1) z - a stochastic shock
        !   2) k - the k at which to evaluate the output
        !   3) l - the l at which to evaluate the output
        real(DP) :: sub_modelOutput
        real(DP), dimension(:), intent(IN) :: inputVars
        real(DP) :: z, k

        z=inputVars(1)
        k=inputVars(2)

        sub_modelOutput = exp(z)*k**alpha+(1-delta)*k
    end function sub_modelOutput

    function sub_modelf_kp(inputVars)
        ! a function that returns the marginal productivity of capital, given a set
        ! of required inputs. The typical minimal inputs are:
        !   1) z - a stochastic shock
        !   2) k - the k at which to evaluate the output
        !   3) l - the l at which to evaluate the output

        real(DP) :: sub_modelf_kp
        real(DP), dimension(:), intent(IN) :: inputVars
        real(DP) :: z, k

        z=inputVars(1)
        k=inputVars(2)

        sub_modelf_kp = -alpha*exp(z)*k**(alpha-1)-(1-delta)
    end function sub_modelf_kp

    function sub_modelf_lp(inputVars)
        ! a function that returns the marginal productivity of labour, given a set
        ! of required inputs. The typical minimal inputs are:
        !   1) z - a stochastic shock
        !   2) k - the k at which to evaluate the output
        !   3) l - the l at which to evaluate the output
        real(DP) :: sub_modelf_lp
        real(DP), dimension(:), intent(IN) :: inputVars
        real(DP) :: z, k, l

        z=inputVars(1)
        k=inputVars(2)
        l=inputVars(3)

        sub_modelf_lp = (1-alpha)*exp(z)*k**(alpha)*l**(-alpha)
    end function sub_modelf_lp

    subroutine steadystate(outputVars)
            ! This subroutine contains the calculations for steady state values.
            !
            ! OUTPUTS: outputVars - a vector of output values.
            !
            ! Note: To pass arrays, see the initValueFn subroutine below
        implicit none

        real(DP), dimension(:), INTENT(OUT) :: outputVars

#ifdef HAS_LABOUR
    real(DP) :: phi, omega, big_phi

    phi     = (1/alpha*(1/beta-1+delta))**(1/(1-alpha))
    omega   = (phi**(1-alpha)-delta)
    big_phi = theta/(1-theta)*(1-alpha)*phi**(-alpha)

    k_ss = big_phi/(omega+phi*big_phi)
    l_ss = phi*k_ss
    c_ss = omega*k_ss
    r_ss = alpha*(k_ss**(alpha-1))*(l_ss**(1-alpha))
    y_ss = k_ss**alpha*(l_ss**(1-alpha))

    outputVars(1) = y_ss
    outputVars(2) = c_ss
    outputVars(3) = k_ss
    outputVars(4) = l_ss
    outputVars(5) = r_ss

    print *, 'Steady State values'
    print *, 'Output: ', y_ss, 'Capital: ', k_ss, 'Labor: ', l_ss

#else
        k_ss = ((1/beta-1+delta)/alpha)**(1/(alpha-1))
        y_ss = k_ss**alpha
        c_ss = k_ss**alpha-delta*k_ss
        outputVars(1) = y_ss
        outputVars(2) = c_ss
        outputVars(3) = k_ss

        print *, 'Steady State values'
        print *, 'Output: ', y_ss, 'Capital: ', k_ss

#endif
    end subroutine steadystate

    subroutine initValueFn(inputVars, outputVars)
            ! This subroutine contains the calculations to initialize the initial guess of an increasing
            ! value function.
            !
            ! INPUTS: inputVars - a vector of input values. typicaly
            !                         1) step1:  the step size
            !                         2) m: the number of points on the k grid
            !                         3) n: the number of exogenous states
            ! OUTPUTS: outputVars - a vector (which is actually a linearized matrix) representing the value function
            !
            ! Note: To pass arrays, see the initValueFn subroutine below
        implicit none

        real(DP), dimension(:), INTENT(IN) :: inputVars
        real(DP), dimension(:), INTENT(OUT) :: outputVars

#ifdef HAS_LABOUR
    real(DP) :: step1, valueinitial, theta
    integer :: index_k, i1, m, n

    step1 = inputVars(1)
    m = inputVars(2)
    n = inputVars(3)

    valueinitial = (1/(1-beta))*((c_ss**theta*(1-l_ss)**theta)**(1-tau))/(1-tau)
    forall (index_k = 1:m, i1=1:n) outputVars((index_k-1)*n+i1) = index_k/step1 + valueinitial

#else
        real(DP) :: step1, valueinitial
        integer :: index_k, i1, m, n

        step1 = inputVars(1)
        m = inputVars(2)
        n = inputVars(3)

        valueinitial = (1/(1-beta))*c_ss**(1-tau)/(1-tau)

        forall (index_k = 1:m, i1=1:n) outputVars((index_k-1)*n+i1) = index_k/step1 + valueinitial
#endif
    end subroutine initValueFn

    SUBROUTINE cashInHand(outputVars, m, n, grid_y, grid_k, grid_z)
            !
            ! This subroutine defines tomorrow's "cash in hand" using model parameters.
            !
            ! INPUTS: m: the number of points on the k grid
            !         n: the number of exogenous states
            !         grid_k: the vector of size m of capital grid points
            !         grid_y: The vector of size n of output given shocks
            !         grid_z: The vector of discretized shocks (optional)
            ! OUTPUTS: outputVars - a m-by-n matrix representing the cash in hand function
        real(DP), dimension(:), INTENT(IN) :: grid_k
        real(DP), dimension(:), INTENT(IN) :: grid_y
        REAL(DP), dimension(:), INTENT(IN), OPTIONAL :: grid_z
        integer, INTENT(IN) :: m, n
        real(DP), dimension(:,:), INTENT(OUT) :: outputVars

        integer :: index_k, index_z

        IF (size(grid_k,dim=1)/=m) THEN
            PRINT '(a,i3)', 'cashInHand: grid_k must be a vector of size ',m
            call sub_modelstop('program terminated in cash in hand')
        END IF

        IF (size(grid_y,dim=1)/=n) THEN
            PRINT '(a,i3)', 'cashInHand: grid_y must be a vector of size ',n
            call sub_modelstop('program terminated in cash in hand')
        END IF

        IF (present(grid_z)) THEN
            IF (size(grid_z,dim=1)/=n) THEN
                PRINT '(a,i3)', 'cashInHand: grid_z must be a vector of size ',n
                call sub_modelstop('program terminated in cash in hand')
            END IF
        END IF

        IF ( (size(outputVars,dim=1)/=m) .and. (size(outputVars,dim=2)/=n) ) THEN
            PRINT '(a,i3)', 'cashInHand: outputVars must be a matrix of size ',m, '-by-',n
            call sub_modelstop('program terminated in cash in hand')
        END IF

#ifndef HAS_LABOUR
        forall (index_k = 1:m, index_z=1:n)
            outputVars(index_k,index_z)=&
                &exp(grid_y(index_z))*grid_k(index_k)**alpha+(1-delta)*grid_k(index_k)
        end forall
#else
            PRINT '(a,i3)', 'cashInHand: has_labour, but cash in hand not defined. ERROR!'
            call sub_modelstop('program terminated in cash in hand')
#endif
    END SUBROUTINE cashInHand

    FUNCTION getCStar(length_grid_k,n_tauchen,valuefn,stepVec) RESULT(y)
        ! Get optimal consumption given the value function
        ! Note: The last step is from the envelope condition
        INTEGER, INTENT(IN) :: length_grid_k
        INTEGER, INTENT(IN) :: n_tauchen
        REAL(DP), DIMENSION(length_grid_k,n_tauchen), INTENT(IN) :: valuefn
        REAL(DP), DIMENSION(length_grid_k), INTENT(IN) :: stepVec
        REAL(DP), DIMENSION(length_grid_k,n_tauchen) :: y
        REAL(DP), DIMENSION(length_grid_k,n_tauchen) :: D
        INTEGER :: index_k

        D(1,:)=(valuefn(2,:)-valuefn(1,:))/stepVec(1)
        D(length_grid_k,:)=(valuefn(length_grid_k,:)-valuefn(length_grid_k-1,:))/stepVec(length_grid_k-1)
        do index_k = 2,length_grid_k-1
            D(index_k,:)=(valuefn(index_k+1,:)-valuefn(index_k-1,:))/(stepVec(index_k)+stepVec(index_k-1))
        enddo
        y=D**(-1/tau)
    END FUNCTION getCStar

    FUNCTION getVStar(length_grid_k,n_tauchen,Cstar,valuefn) RESULT(y)
        ! Get optimal consumption given the value function
        ! Note: The last step is from the envelope condition
        INTEGER, INTENT(IN) :: length_grid_k
        INTEGER, INTENT(IN) :: n_tauchen
        REAL(DP), DIMENSION(length_grid_k,n_tauchen), INTENT(IN) :: valuefn, Cstar
        REAL(DP), DIMENSION(length_grid_k,n_tauchen) :: y

        y=1/(1-tau)*Cstar**(1-tau)+valuefn
    END FUNCTION


END MODULE modelDefinition

! This module contains everything we neeed to perform a value function interation using the endogenous
! grid method. The main program below gives an example of usage.
!
! USAGE: You must define all the routines in the modelDefinition module above.
MODULE  endogenousGrid
    use nrtype
    use modelDefinition
    implicit none


    !----------------------------------------------------------------
    ! 1. Calibration
    !----------------------------------------------------------------

    real (DP), parameter :: delta  = 0.0196                         ! Depreciation
    real (DP), parameter :: alpha  = 0.4                            ! Capital Share
    real (DP), parameter :: tau    = 2.d0                           ! risk aversion
    real (DP), parameter :: beta   = 0.9896                         ! Discount factor
    real (DP), parameter :: theta   = 0.357                         ! weight on consumption
    real (DP), parameter :: rho    = 0.95                           ! Autoregressive
    real (DP), parameter :: sigma  = 0.007                          ! Variance

    !----------------------------------------------------------------
    ! 2. Numerical parameters
    !----------------------------------------------------------------

    integer :: n_tauchen = 41                                   ! Number Tauchen points
    integer :: Tsimul    = 50000                                ! Number of simulations for the EEerror

    real (DP), parameter :: toler   = 1e-6                      ! Numerical tolerance
    integer, parameter :: MAX_ITER = 5000                       ! max number of iterations for convergence
    real (DP), allocatable, dimension (:) :: stepVec
    integer :: i1

    private :: sub_myinterp1, sub_kendogenousnewton
    private :: i1

contains

    subroutine sub_mystop(calling)
        ! a personal stop subroutine. Makes it easier to edit behaviour of stop. All
        ! functions and subroutines should call this.
        !
        ! INPUTS: calling - a string indicating where this subroutine was called from

        CHARACTER (LEN=*), intent(in) :: calling
        print *, "STOP: ", calling
        call clean()
        STOP 0
    end subroutine sub_mystop

    subroutine initialize(steps)
        integer, intent(in) :: steps
        allocate(stepVec(steps))
    end subroutine initialize

    subroutine clean()
        deallocate(stepVec)
    end subroutine clean


    subroutine normal_01_cdf ( x, cdf )

        !*****************************************************************************80
        !
        !! NORMAL_01_CDF evaluates the Normal 01 CDF.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    10 February 1999
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Reference:
        !
        !    AG Adams,
        !    Algorithm 39,
        !    Areas Under the Normal Curve,
        !    Computer Journal,
        !    Volume 12, pages 197-198, 1969.
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) X, the argument of the CDF.
        !
        !    Output, real ( kind = 8 ) CDF, the value of the CDF.
        !
        implicit none

        real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
        real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
        real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
        real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
        real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
        real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
        real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
        real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
        real ( kind = 8 ), parameter :: b1 = 3.8052D-08
        real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
        real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
        real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
        real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
        real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
        real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
        real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
        real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
        real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
        real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
        real ( kind = 8 ) cdf
        real ( kind = 8 ) q
        real ( kind = 8 ) x
        real ( kind = 8 ) y
        !
        !  |X| <= 1.28.
        !
        if ( abs ( x ) <= 1.28D+00 ) then

            y = 0.5D+00 * x * x

            q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
                + a6 / ( y + a7 ) ) ) )
        !
        !  1.28 < |X| <= 12.7
        !
        else if ( abs ( x ) <= 12.7D+00 ) then

            y = 0.5D+00 * x * x

            q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
                + b2 / ( abs ( x ) + b3 &
                + b4 / ( abs ( x ) - b5 &
                + b6 / ( abs ( x ) + b7 &
                - b8 / ( abs ( x ) + b9 &
                + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
        !
        !  12.7 < |X|
        !
        else

            q = 0.0D+00

        end if
        !
        !  Take account of negative X.
        !
        if ( x < 0.0D+00 ) then
            cdf = q
        else
            cdf = 1.0D+00 - q
        end if

        return
    end subroutine normal_01_cdf

    subroutine sub_tauchen(y,pi,picum,rho,sigma,cover_tauchen)
        ! Tauchen method to approximate univariate AR(1) process by Markov chain
        !      x(t) = rho x(t-1)+ sigma sqrt(1-rho^2) e(t),   e(t)~N(0,1)
        ! INPUTS: rho - serial correlation coefficient,
        !         sigma - coefficient of variation
        !         cover_tauchen - a spread variable for the sigma, not really sure about this
        ! OUTPUT: y is the income states
        !         pi is an n-by-n matrix of Markov transition probabilities
        !         picum is an n-by-n matrix representing cdf of pi


        implicit none

        real(DP), intent(in) :: rho, sigma, cover_tauchen
        real(DP), dimension (:,:), intent(out) :: pi , picum
        real(DP), dimension (:), intent(out) :: y

        integer :: yc, tc, tcc, n, i1

        real(DP) :: sigy, mu, upp, low, temp, temp2

        !
        ! basic error checking
        !
        n=size(y,dim=1)
        IF (size(y,dim=1)/=n) THEN
            PRINT '(a,i3)', 'tauchen: y must be a vector of size ',n
            call sub_mystop('program terminated by tauchen')
        END IF
        IF ((size(pi,dim=1)/=n) .or. (size(pi,dim=2)/=n)) THEN
            PRINT '(a,i3,a,i3)', 'tauchen: pi must be a square matrix of size ',n,' x ',n
            call sub_mystop('program terminated by tauchen')
        END IF
        IF ((size(picum,dim=1)/=n) .or. (size(picum,dim=2)/=n)) THEN
            PRINT '(a,i3,a,i3)', 'tauchen: picum must be a square matrix of size ',n,' x ',n
            call sub_mystop('program terminated by tauchen')
        END IF

        !
        ! Define the discrete states of the Markov chain
        !
        sigy = sigma/sqrt(1.0-rho**2)

        y(n) = cover_tauchen*sigy
        y(1) = -y(n)

        do yc = 2, n-1
            y(yc) = y(1)+(y(n)-y(1))*(yc-1.0)/(n-1.0)
        end do

        !
        ! Compute the transition matrix
        !
        do tc = 1,n
            mu = rho*y(tc)
            upp = (y(2)+y(1))/2.0
            upp = (upp-mu)/sigma
            call normal_01_cdf(upp, temp)
            pi(tc,1) = temp
            low = (y(n)+y(n-1))/2.0
            low = (low-mu)/sigma
            call  normal_01_cdf(low, temp)
            pi(tc,n) = 1.0-temp

            do tcc = 2, n-1
                low = (y(tcc-1)+y(tcc))/2.0
                upp = (y(tcc+1)+y(tcc))/2.0
                low = (low-mu)/sigma
                upp = (upp-mu)/sigma
                call normal_01_cdf(upp, temp)
                call normal_01_cdf(low, temp2)
                pi(tc,tcc) = temp-temp2

            end do
        end do

        !
        ! Compute the CDF of the transition matrix.
        !
        picum(:,1)=pi(:,1)
        do tc = 2,n
            picum(:,tc)=picum(:,tc-1)+pi(:,tc)
        enddo

        open(unit=98,file='Tauchen.txt')
        write(98,*) 'rho, sigma, n'
        write(98,*) rho, sigma, n
        write(98,*) 'Income States'
        write(98,*) y(1:n)
        write(98,*) 'Transition Matrix'

        do tc=1,n
            write(98,*) pi(tc,:)
        end do

        write(98,*) 'Transition Matrix CDF'

        do tc=1,n
            write(98,*) picum(tc,:)
        end do

        close(98)

        open(unit=16,   file='transition.txt', status = 'replace')
        write (16, '(41f20.10)') (pi(i1, :), i1 = 1,n)
        close(16)

    end subroutine sub_tauchen

    SUBROUTINE sub_grid_generation(x,xcentre,xbounds,s)
        ! Purpose: Generate grid x on [xcentre*(1-xbounds),xcentre*(1+xbounds)] using spacing parameter s set as follows:
        ! s=1       linear spacing
        ! s>1       left skewed grid spacing with power s
        ! 0<s<1     right skewed grid spacing with power s
        ! s<0       geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
        ! s=-1      logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
        ! s=0       logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
        !
        ! INPUTS: xcentre - The centre of the grid
        !         xbounds - how far away from the centre the bounds go.
        !         s - skewness of grid - see above
        ! OUTPUT: x - the generated grid

        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(OUT) :: x
        REAL(DP), INTENT(IN) :: xcentre,xbounds, s
        REAL(DP) :: c ! growth rate of grid subintervals for logarithmic spacing
        REAL(DP) :: xmax, xmin
        INTEGER :: n,i

        n=size(x)
        xmax=xcentre*(1+xbounds);
        xmin=xcentre*(1-xbounds);

        FORALL(i=1:n) x(i)=(i-1)/real(n-1,WP)
        IF (s>0.0_WP) THEN
            x=x**s*(xmax-xmin)+xmin
        ELSE
            IF (s==-1.0_WP) THEN
                c=xmax-xmin+1
            ELSE
                c=-s
            END IF
            x=((xmax-xmin)/(c-1))*(c**x)-((xmax-c*xmin)/(c-1))
        END IF
        FORALL(i=2:n) stepVec(i-1) = x(i)-x(i-1)

    END SUBROUTINE sub_grid_generation

    subroutine sub_myinterp1(x,f_x,xp,length_x,interp_value)
        ! Purpose:
        !
        ! INPUTS: x
        !         f_x
        !         xp
        ! OUTPUT: interp_value

        implicit none

        integer, intent(in) :: length_x

        real(DP), dimension(length_x), intent(in) :: x, f_x, xp
        real(DP), dimension(length_x), intent(out) :: interp_value

        integer, dimension(1) :: x_min
        integer :: index_k, x_index
        real (DP) :: t


        do index_k = 1,length_x

            x_min = minloc(abs(x-xp(index_k)))
            x_index = x_min(1)

            if (xp(index_k)<x(x_index)) x_index = x_index-1
            x_index = max(x_index,1)
            x_index = min(x_index,length_x-1)
            t = (f_x(x_index+1)-f_x(x_index))/(x(x_index+1)-x(x_index))
            interp_value(index_k) = t*(xp(index_k)-x(x_index))+f_x(x_index)

        enddo
    end subroutine sub_myinterp1

    subroutine sub_myinterp3(x,f_x,length_index,xp,interp_value)

        implicit none

        integer, intent(in) :: length_index

        real (dp), dimension(length_index), intent(in) :: x
        real (dp), dimension(length_index, n_tauchen), intent(in) :: f_x
        real (dp), intent(in) :: xp
        real (dp), dimension(n_tauchen), intent(out) :: interp_value

        real (dp), dimension(n_tauchen) :: t

        if (xp>x(2) .and. length_index == 3) then
            t = (f_x(3,:)-f_x(2,:))/(x(3)-x(2))
            interp_value = t*(xp-x(2))+f_x(2,:)
        else
            t = (f_x(2,:)-f_x(1,:))/(x(2)-x(1))
            interp_value = t*(xp-x(1))+f_x(1,:)
        endif

    end subroutine sub_myinterp3


    subroutine sub_kendogenousnewton(kend,Yend,z,kguess)

        implicit none

        real(DP), intent(in) :: Yend,z,kguess
        real(DP), intent(out) :: kend
        real(DP), dimension(2) :: inputVars
        integer :: go_on,iter

        real(DP) :: f_k, f_kp, kp

        kend=kguess

        inputVars(1)=z

        go_on = 1
        iter=0
        do while (go_on == 1)
            iter = iter+1

             ! difference between target production (Yend) and estimate at
             ! given k (kend)
            inputVars(2)=kend
            f_k  = Yend-sub_modelOutput(inputVars)
            f_kp = sub_modelf_kp(inputVars)
            kp = kend-(f_k/f_kp)

            if (abs(kp-kend)<toler) go_on = 0

            kend = kp
            if(iter>MAX_ITER) then
                PRINT '(a,i3,a,i3)', 'sub_kendogenousnewton: too many iterations'
                call sub_mystop('program terminated by sub_kendogenousnewton')
            endif

        enddo

    end subroutine sub_kendogenousnewton

    subroutine sub_kendogenousnewton2(kend,z,labor,Cstar,kprime)

        implicit none

        real(dp), intent(in) :: z,labor,Cstar,kprime
        real(dp), intent(out) :: kend

        integer :: go_on

        real(dp) :: f_k, f_kp, kp

        !Uses as initial guess the value of capital tomorrow.
        kend=kprime

        go_on = 1
        do while (go_on == 1)

            f_k  = kprime + Cstar-exp(z)*kend**alpha*labor**(1-alpha)-(1-delta)*kend
            f_kp = -alpha*exp(z)*kend**(alpha-1)*labor**(1-alpha)-(1-delta)
            kp = kend-(f_k/f_kp)

            if (abs(kp-kend)<toler) go_on = 0
            kend = kp
        enddo

    end subroutine sub_kendogenousnewton2

    subroutine sub_valueend1(valuefn,grid_k,length_grid_k,y,transition,l_ss)
                ! The actual value iteration function when we have labour
                ! INPUTS: valuefn - The initial guess of the value function, size m-by-n
                !         grid_k - grid of possible capital values (..,ks,..), size m
                !         length_grid_k - basically, m
                !         y - the income states - size n
                !         transition - the transition matrix - size n-by-n
                !         l_ss - the steady state value of labour
                ! OUTPUT: valuefn - the final value function - size m-by-n


        implicit none

        integer, intent(in) :: length_grid_k

        real(dp), dimension(length_grid_k,n_tauchen), intent(inout) :: valuefn
        real(dp), dimension(n_tauchen,n_tauchen), intent(in) :: transition
        real(dp), dimension(length_grid_k), intent(in) :: grid_k
        real(dp), dimension(n_tauchen), intent(in) :: y
        real(dp), intent(in) :: l_ss

        integer :: iter, index_k, index_z

        real(dp), dimension(length_grid_k,n_tauchen) :: cih, D, Cstar, Vend, Vend1, Yend, value1
        real(dp) :: diff, B1, B2, B3

        iter = 0
        diff = 1000.d0

        ! Need to define "cash in hand" tomorrow.
        ! BLEH! function specific. Still want to change this!
        do index_k = 1,length_grid_k
            do index_z = 1,n_tauchen
                cih(index_k,index_z)=exp(y(index_z))*grid_k(index_k)**alpha*l_ss**(1-alpha)+(1-delta)*grid_k(index_k)
            enddo
        enddo

        ! These parameters simplify the calculations later

        B1=(theta-1)*(1-tau)
        B2=1/(theta*(1-tau)-1)
        B3=theta*(1-tau)

        !Here we start the iterations.

        iter=0
        diff=1

        do while((diff>toler) .and. (iter<MAX_ITER))
            iter=iter+1
            !Compute the derivatives of Vtilda at the grid points only.
            !
            ! A smarter way to do this - use the Euler condition. Would be nice if
            ! we could just pass in a function pointer to the euler condition and calculate this.
            D(1,:)=(valuefn(2,:)-valuefn(1,:))/stepVec(1)
            D(length_grid_k,:)=(valuefn(length_grid_k,:)-valuefn(length_grid_k-1,:))/stepVec(length_grid_k-1)
            do index_k = 2,length_grid_k-1
                D(index_k,:)=(valuefn(index_k+1,:)-valuefn(index_k-1,:))/(stepVec(index_k)+stepVec(index_k-1))
            enddo

            Cstar=(1/theta*D*(1-l_ss)**B1)**B2
            Vend=1/(1-tau)*Cstar**B3*(1-l_ss)**(-B1)+valuefn

            do index_z = 1,n_tauchen
                Yend(:,index_z)=Cstar(:,index_z)+grid_k
                call sub_myinterp1(Yend(:,index_z),Vend(:,index_z), cih(:,index_z), length_grid_k, Vend1(:,index_z))
            enddo

            do index_k = 1,length_grid_k
                do index_z = 1,n_tauchen
                    value1(index_k,index_z)=beta*dot_product(transition(index_z,:),Vend1(index_k,:))
                enddo
            enddo

            diff = maxval(abs(valuefn-value1))

            valuefn=value1
            if (mod(iter,10)==0) then
                print *, 'sub_valueend1 Iteration: ', iter, 'Tolerance ', diff
                flush(6)
            endif
        enddo

        !The only variable that we need out of this program is Vend1 since it is already defined on grid_k.
        valuefn=Vend1*(1-beta)

    end subroutine sub_valueend1

    subroutine sub_valueend2(valuefn,grid_k,length_grid_k,y,transition,kend,g_l,diffstandard)

        implicit none

        integer, intent(in) :: length_grid_k

        real(dp), dimension(length_grid_k,n_tauchen), intent(inout) :: valuefn
        real(dp), dimension(length_grid_k,n_tauchen), intent(inout) :: kend, g_l
        real(dp), dimension(n_tauchen,n_tauchen), intent(in) :: transition
        real(dp), dimension(length_grid_k), intent(in) :: grid_k
        real(dp), dimension(n_tauchen), intent(in) :: y
        real(dp), intent(in) :: diffstandard

        integer :: iter, index_k, index_z, deriv

        real(dp), dimension(length_grid_k,n_tauchen) :: D, Cstar, Vend, Vend1, value1, L1
        real(dp), dimension(length_grid_k,4) :: X
        real(dp), dimension(4,4) :: XXinv
        real(dp), dimension(4) :: betahat
        real(dp) :: diff, B1, B2, B3

        iter = 0
        diff = 1000.d0

        ! Need to interpolate the policy function and the value function at the new Kend values

        do index_z = 1,n_tauchen
            call sub_myinterp1(grid_k, g_l(:,index_z), kend(:,index_z), length_grid_k, L1(:,index_z))
        enddo
        Vend1=valuefn

        ! These parameters simplify the calculations later

        B1=(theta-1)*(1-tau)
        B2=1/(theta*(1-tau)-1)
        B3=theta*(1-tau)

        !Options for the derivative
        deriv=1 !1 is the usual method, 2 is the polynomial.

        !Here we start the iterations.
        iter=0
        diff=1

        do while(diff>0.1d0*diffstandard)

            iter=iter+1
            do index_k = 1,length_grid_k
                do index_z = 1,n_tauchen
                    valuefn(index_k,index_z)=beta*dot_product(transition(index_z,:),Vend1(index_k,:))
                enddo
            enddo

            !Compute the derivatives of Vtilda at the grid points only.
            if (deriv==1) then
                !Compute the derivatives of Vtilda at the grid points only.
                !
                ! A smarter way to do this - use the Euler condition. Would be nice if
                ! we could just pass in a function pointer to the euler condition and calculate this.
                D(1,:)=(valuefn(2,:)-valuefn(1,:))/stepVec(1)
                D(length_grid_k,:)=(valuefn(length_grid_k,:)-valuefn(length_grid_k-1,:))/stepVec(length_grid_k-1)
                do index_k = 2,length_grid_k-1
                    D(index_k,:)=(valuefn(index_k+1,:)-valuefn(index_k-1,:))/(stepVec(index_k)+stepVec(index_k-1))
                enddo
            else
                !Derivative using the regression.
                do index_z =1,n_tauchen
                    betahat = matmul(XXinv,matmul(transpose(X),valuefn(:,index_z)))
                    D(:,index_z) = betahat(2) + 2.d0*betahat(3)*grid_k + 3.d0*betahat(4)*grid_k**2
                enddo
            endif

            Cstar=(1/theta*D*(1-L1)**B1)**B2
            Vend=1/(1-tau)*Cstar**B3*(1-L1)**(-B1)+valuefn
            do index_k = 1,length_grid_k
                do index_z = 1,n_tauchen
                    call sub_kendogenousnewton2(kend(index_k,index_z),y(index_z),L1(index_k,index_z),&
                        & Cstar(index_k,index_z),grid_k(index_k))
                enddo
            enddo

            do index_z = 1,n_tauchen
                call sub_myinterp1(kend(:,index_z), Vend(:,index_z), grid_k, length_grid_k, value1(:,index_z))
                call sub_myinterp1(grid_k, g_l(:,index_z), kend(:,index_z), length_grid_k, L1(:,index_z))
            enddo

            diff = maxval(abs(value1-Vend1))
            Vend1=value1
#if 0
            if (mod(iter,10)==0) then
                print *, 'sub_valueEnd2 - Iteration: ', iter, 'Tolerance ', diff
                flush(6)
            endif
#endif
        enddo

        !The only variable that we need out of this program is Vend1 since it is already defined on grid_k.
        valuefn=value1*(1-beta)

    end subroutine sub_valueend2

    subroutine sub_value(valuefn, g_k, g_c, grid_k, y, transition)
            ! The actual value iteration function
            ! INPUTS: valuefn - The initial guess of the value function, size m-by-n
            !         grid_k - grid of possible capital values (..,ks,..), size m
            !         y - the income states - size n
            !         transition - the transition matrix - size n-by-n
            ! OUTPUT: g_c - interpolated value of consumption - size m-by-n
            !         g_k - interpolated value of capital - size m-by-n
            !         valuefn - the final value function - size m-by-n

        implicit none

        real(DP), dimension(:,:), intent(inout) :: valuefn
        real(DP), dimension(:,:), intent(out) :: g_k, g_c
        real(DP), dimension(:,:), intent(in) :: transition
        real(DP), dimension(:), intent(in) :: grid_k
        real(DP), dimension(:), intent(in) :: y

        integer :: iter, index_k, index_z
        integer :: length_grid_k,m,n

        real(DP), allocatable, dimension(:,:) :: cih, D, Cstar, Vend, Vend1, Yend, value1, kend
        real(DP) :: diff

        length_grid_k=size(grid_k)
        m=length_grid_k
        n=n_tauchen
        !
        !do size checks
        !
        IF ((size(valuefn,dim=1)/=m) .or. (size(valuefn,dim=2)/=n)) THEN
            PRINT '(a,i3,a,i3)', 'sub_value: valuefn must be a matrix of size ',m,' x ',n
            call sub_mystop('program terminated by sub_value')
        END IF
        IF (size(y,dim=1)/=n) THEN
            PRINT '(a,i3)', 'sub_value: y must be a vector of size ',n
            call sub_mystop('program terminated by sub_value')
        END IF
        IF (size(transition,dim=1)/=n .or. size(transition,dim=2)/=n) THEN
            PRINT '(a,i3,a,i3)', 'sub_value: transition must be a square matrix of size ',n,' x ',n
            call sub_mystop('program terminated by sub_value')
        END IF
        IF ((size(g_k,dim=1)/=m) .or. (size(g_k,dim=2)/=n)) THEN
            PRINT '(a,i3,a,i3)', 'sub_value: g_k must be a matrix of size ',m,' x ',n
            call sub_mystop('program terminated by sub_value')
        END IF
        IF ((size(g_c,dim=1)/=m) .or. (size(g_c,dim=2)/=n)) THEN
            PRINT '(a,i3,a,i3)', 'sub_value: g_c must be a matrix of size ',m,' x ',n
            call sub_mystop('program terminated by sub_value')
        END IF

        !
        !allocate arrays
        !
        allocate(cih(m,n))
        allocate(D(m,n))
        allocate(Cstar(m,n))
        allocate(Vend(m,n))
        allocate(Vend1(m,n))
        allocate(Yend(m,n))
        allocate(value1(m,n))
        allocate(kend(m,n))

        iter = 0
        diff = 1000.d0

        ! define the cash in hand array
        call cashInHand(cih, length_grid_k, n_tauchen, y, grid_k)

        !Here we start the iterations.
        print *,"starting iterations"
        flush(6)
        do while((diff>toler) .and. (iter<MAX_ITER))

            iter=iter+1

            !Get optimal consumption
            !
            Cstar = getCStar(length_grid_k,n_tauchen,valuefn,stepVec)

            !get optimal value function
            Vend=getVStar(length_grid_k,n_tauchen,Cstar,valuefn)

            do index_z = 1,n_tauchen
                Yend(:,index_z)=Cstar(:,index_z)+grid_k
                call sub_myinterp1(Yend(:,index_z),Vend(:,index_z), cih(:,index_z), length_grid_k, Vend1(:,index_z))
            enddo

            do index_k = 1,length_grid_k
                do index_z = 1,n_tauchen
                    value1(index_k,index_z)=beta*dot_product(transition(index_z,:),Vend1(index_k,:))
                enddo
            enddo

            diff = maxval(abs(value1-valuefn))
            valuefn=value1
            if (mod(iter,50)==0) then
                print *, 'subvalue: Iteration: ', iter, 'Tolerance ', diff
                flush(6)
            end if
        enddo

        ! At this point the program recovers the endogenous grid for capital.
        do index_k = 1,length_grid_k
            do index_z = 1,n_tauchen
                call sub_kendogenousnewton(kend(index_k,index_z),Yend(index_k,index_z),y(index_z),grid_k(index_k))
            enddo
        enddo

        ! Here we interpolate capital on the grid so that we can compare the results to the standard algorithm.
        do index_z = 1,n_tauchen
            call sub_myinterp1(kend(:,index_z),Vend(:,index_z), grid_k, length_grid_k, valuefn(:,index_z))
            call sub_myinterp1(kend(:,index_z),grid_k, grid_k, length_grid_k, g_k(:,index_z))
            call sub_myinterp1(kend(:,index_z),Cstar(:,index_z), grid_k, length_grid_k, g_c(:,index_z))
        enddo

        !
        !deallocate arrays
        !
        deallocate(cih)
        deallocate(D)
        deallocate(Cstar)
        deallocate(Vend)
        deallocate(Vend1)
        deallocate(Yend)
        deallocate(value1)
        deallocate(kend)

    end subroutine sub_value

    subroutine sub_newton(kp,k,z,l)

        implicit none

        real(dp), intent(in) :: kp, k, z
        real(dp), intent(out) :: l

        integer :: go_on, iter

        real(dp) :: f_l, f_lp, lp

        l = 0.3d0
        go_on = 1

        iter = 1
        do while ((go_on == 1) .and. (iter<MAX_ITER))
            iter=iter+1

            f_l  = kp-exp(z)*(k**alpha)*(l**(1-alpha))-(1-delta)*k+(theta/(1-theta))*(1-alpha)*exp(z)*(k**alpha)*(l**(-alpha))*(1-l)
            f_lp = -(1-alpha)*exp(z)*(k**alpha)*(l**(-alpha))+(-alpha)*(theta/(1-theta))*(1-alpha)*exp(z)*(k**alpha)&
                &*(l**(-1-alpha))*(1-l)-(theta/(1-theta))*(1-alpha)*exp(z)*(k**alpha)*(l**(-alpha))
            lp = l-(f_l/f_lp)

            if (abs(lp-l)<toler .or. lp>0.95 .or. lp<0.01) go_on = 0
            l = lp

            if (mod(iter,10)==0) then
                print *, 'subnewton: Iteration: ', iter, 'Tolerance ', abs(lp-l)
                flush(6)
            end if

        enddo

        l = min(0.95,l)
        l = max(0.01,l)

    end subroutine sub_newton

    subroutine sub_refinement(index_kp,index_z,k,z,grid_k,length_grid_k,transition,value,kp,mylabor,value_max_so_far)

        implicit none

        integer, intent(in) :: length_grid_k, index_kp, index_z

        real (dp), dimension(length_grid_k,n_tauchen), intent(in) :: value
        real (dp), dimension(n_tauchen,n_tauchen), intent(in) :: transition
        real (dp), dimension(length_grid_k), intent(in) :: grid_k
        real (dp), intent(in) :: k, z
        real (dp), intent(out) :: kp, mylabor, value_max_so_far

        integer :: length_index, index_kp_low, index_kp_high, index_kp_now, go_on, iter

        real (dp), dimension(n_tauchen) :: myvalue
        real (dp) :: c, kp_low, kp_high, kp_guess, kp_guess_new, value_comp, mylabor_guess

        index_kp_now = max(index_kp-1,1)
        kp = grid_k(index_kp_now)

        call sub_newton(kp,k,z,mylabor)

        c = exp(z)*(k**alpha)*(mylabor**(1-alpha))+(1-delta)*k-kp
        value_max_so_far = (1-beta)*(((c**theta)*((1-mylabor)**(1-theta)))**(1-tau))/(1-tau)+&
            &beta*dot_product(transition(index_z,:),value(index_kp_now,:))

        index_kp_low  = max(1,index_kp-2)
        index_kp_high = min(length_grid_k,index_kp)
        length_index = index_kp_high-index_kp_low+1

        kp_low  = grid_k(index_kp_low)
        kp_high = grid_k(index_kp_high)
        kp_guess = kp_low+0.6*(kp_high-kp_low)

        go_on = 1
        iter = 0
        do while ((go_on == 1) .and. (iter<MAX_ITER))
            iter = iter+1

            call sub_myinterp3(grid_k(index_kp_low:index_kp_high),value(index_kp_low:index_kp_high,:),length_index,kp_guess,myvalue)

            call sub_newton(kp_guess,k,z,mylabor_guess)

            c = exp(z)*(k**alpha)*(mylabor_guess**(1-alpha))+(1-delta)*k-kp_guess

            if (c<0.0) then
                c = toler/10
            endif

            value_comp = (1-beta)*(((c**theta)*((1-mylabor_guess)**(1-theta)))**(1-tau))/(1-tau)+&
                &beta*dot_product(transition(index_z,:),myvalue)

            if (value_comp>value_max_so_far) then
                if (kp_guess<kp) then
                    kp_high = kp
                else
                    kp_low = kp
                endif
                kp_guess_new = kp_low+0.6*(kp_high-kp_low)
                kp = kp_guess
                value_max_so_far = value_comp
            else
                if (kp_guess<kp) then
                    kp_low = kp_guess
                else
                    kp_high = kp_guess
                endif
                kp_guess_new = kp_low+0.6*(kp_high-kp_low)
            endif

            if (abs(kp_high-kp_low)<toler) go_on = 0
            kp_guess = kp_guess_new

#if 0
            if (mod(iter,10)==0) then
                print *, 'subrefinement: Iteration: ', iter, 'Tolerance ', kp_high-kp_low
                flush(6)
            end if
#endif
        enddo

        kp = kp_guess
        mylabor = mylabor_guess

    end subroutine sub_refinement

    subroutine sub_myinterp(x,f_x,xp,length_x,interp_value)

        implicit none

        integer, intent(in) :: length_x

        real(dp), dimension(length_x), intent(in) :: x, f_x
        real(dp), intent(in) :: xp
        real(dp), intent(out) :: interp_value

        integer, dimension(1) :: x_min
        integer :: x_index

        real (dp) :: t

        x_min   = minloc(abs(x-xp))
        x_index = x_min(1)
        if (xp<x(x_index)) x_index = x_index-1
        x_index = max(x_index,1)
        x_index = min(x_index,length_x-1)

        t = (f_x(x_index+1)-f_x(x_index))/(x(x_index+1)-x(x_index))

        interp_value = t*(xp-x(x_index))+f_x(x_index)

    end subroutine sub_myinterp

    subroutine sub_golden(kend, K, grid_k, g_k1, length_grid_k)

        ! Golden Search Algorithm.
        implicit none

        integer, intent(in) :: length_grid_k

        real (dp), dimension(length_grid_k), intent(in) :: grid_k, g_k1
        real (dp), intent(in) :: K
        real (dp), intent(out) :: kend

        integer :: go_on, iter

        real (dp) :: kend1, A, B, C, D, fB, fC

        A=K-1
        D=K+1
        B=0.6*A+0.4*D
        C=0.4*A+0.6*D

        call sub_myinterp(grid_k,g_k1,B,length_grid_k,kend1)
        fB=-abs(kend1-K)
        call sub_myinterp(grid_k,g_k1,C,length_grid_k,kend1)
        fC=-abs(kend1-K)

        go_on=1

        iter = 0
        do while((go_on==1) .and. (iter<MAX_ITER))
            iter = iter+1

            if (fB>fC) then
                D=C
                C=B
                fC=fB
                B=0.6*C+0.4*A
                call sub_myinterp(grid_k,g_k1,B,length_grid_k,kend1)
                fB=-abs(kend1-K)
            else
                A=B
                B=C
                fB=fC
                C=0.6*B+0.4*D
                call sub_myinterp(grid_k,g_k1,C,length_grid_k,kend1)
                fC=-abs(kend1-K)
            endif

            if (abs(D-A)<toler) go_on=0

#if 0
            if (mod(iter,10)==0) then
                print *, 'subgolden: Iteration: ', iter, 'Tolerance ', abs(D-A)
                flush(6)
            end if
#endif

        enddo

        kend=B

    end subroutine sub_golden

    subroutine sub_valuestandard(valuefn,g_k,g_c,g_l,diff,grid_k,length_grid_k,y,transition)

        implicit none

        integer, intent(in) :: length_grid_k
        real(dp), dimension(length_grid_k,n_tauchen), intent(inout) :: valuefn,g_k,g_c,g_l
        real(dp), dimension(n_tauchen,n_tauchen), intent(in) :: transition
        real(dp), dimension(length_grid_k), intent(in) :: grid_k
        real(dp), dimension(n_tauchen), intent(in) :: y
        real(dp), intent(out) :: diff

        integer :: index_k, index_kp, index_z, previous_kp, index_maxk_so_far, flag

        real(dp), dimension(length_grid_k,n_tauchen) :: value_prov, g_lprov, g_kprov
        real(dp) :: diff2, k, kp, z, mylabor, c, value_max_so_far, value_comp, value_refinement, mylaborold


        do index_k = 1, length_grid_k               ! Capital grid
            k = grid_k(index_k)
            do index_z = 1, n_tauchen
                if (index_z == 1) then
                    previous_kp = 1
                else
                    previous_kp = index_maxk_so_far
                endif

                z = y(index_z)
                value_max_so_far = -100000000.d0

                flag = 1
                do index_kp = previous_kp, length_grid_k

                    kp = grid_k(index_kp)

                    call sub_newton(kp,k,z,mylabor)

                    mylaborold = mylabor

                    c = exp(z)*(k**alpha)*(mylabor**(1-alpha))+(1-delta)*k-kp
                    if (c>0) then

                        value_comp = (1-beta)*(((c**theta)*((1-mylabor)**(1-theta)))**(1-tau))/(1-tau)&
                            &+beta*dot_product(transition(index_z,:),valuefn(index_kp,:))

                        if (value_comp>value_max_so_far) then

                            value_max_so_far = value_comp
                            g_lprov(index_k,index_z) = mylabor
                            g_kprov(index_k,index_z) = kp
                            index_maxk_so_far = index_kp
                        else

                            if (flag == 1 .and. length_grid_k>20) then

                                call sub_refinement(index_kp,index_z,k,z,grid_k,length_grid_k,transition,&
                                    &valuefn,kp,mylabor,value_refinement)

                                flag = flag+1

                                if (value_max_so_far<value_refinement) then
                                    value_max_so_far = value_refinement
                                    g_lprov(index_k,index_z) = mylabor
                                    g_kprov(index_k,index_z) = kp
                                endif
                            endif
                        endif
                    else
                        exit
                    endif
                enddo
                value_prov(index_k,index_z) = value_max_so_far
            enddo
        enddo

        diff  = maxval(abs(value_prov - valuefn))
        diff2 = maxval(abs(g_lprov - g_l))
        valuefn = value_prov

        g_l = g_lprov
        g_k = g_kprov
        do index_k = 1, length_grid_k
            do index_z = 1, n_tauchen
                g_c(index_k,index_z)=exp(y(index_z))*grid_k(index_k)**alpha*g_l(index_k,index_z)**(1-alpha)&
                    &+(1-delta)*grid_k(index_k)-g_k(index_k,index_z)
            enddo
        enddo

    end subroutine sub_valuestandard

END MODULE endogenousGrid

program main
    use nrtype
    use endogenousGrid
    use modelDefinition

    implicit none

    real(DP), allocatable, dimension(:,:)  :: transition, transitioncdf, valuefn, g_k, g_c
    real(DP), allocatable, dimension(:)  :: y, grid_k, inputVars, outputVars
    integer :: index_k, length_grid_k, i1
    real(DP) :: k_ss, c_ss, y_ss
    real(DP) :: cover, cover_tauchen, step1

    ! variables used for timing
    real (4) :: elapt, ta(2)
    real(DP) :: t0,t1,t2,delta_t
    ! end timing variables

#ifdef HAS_LABOUR
    real(DP), allocatable, dimension(:,:)  :: g_l, kend
    real(DP) :: l_ss, r_ss
#endif

    call CPU_TIME(t0)
    call CPU_TIME(t1)

    print *, 'Beginning of the Program Value Function Iteration'
    print *, ' '
    flush(6)

    !----------------------------------------------------------------
    ! 1. Computation of the Steady State
    !----------------------------------------------------------------

#ifndef HAS_LABOUR
    allocate(outputVars(3))

    call steadyState(outputVars)

    y_ss = outputVars(1)
    c_ss = outputVars(2)
    k_ss = outputVars(3)

    deallocate(outputVars)

#else
    allocate(outputVars(5))

    call steadyState(outputVars)

    y_ss = outputVars(1)
    c_ss = outputVars(2)
    k_ss = outputVars(3)
    l_ss = outputVars(4)
    r_ss = outputVars(5)

    deallocate(outputVars)
#endif

    !----------------------------------------------------------------
    ! 2. Tauchen Points
    !----------------------------------------------------------------
    allocate(transition(n_tauchen,n_tauchen))
    allocate(transitioncdf(n_tauchen,n_tauchen))
    allocate(y(n_tauchen))

    cover_tauchen = 3.d0
    call sub_tauchen(y,transition,transitioncdf,rho,sigma,cover_tauchen)

    print *, ' '
    print *, 'For the stochastic shock, Tauchen`s approximation method with', n_tauchen , 'points will be used.'
    print *, ' '
    flush(6)

    !----------------------------------------------------------------
    ! 3. Endogenous Grid algorithm
    !----------------------------------------------------------------

    print *, 'Beginning the iteration of the value function with an initial linear guess '
    print *, ' '
    flush(6)

    cover = 0.25d0         ! Coverage of the grid (+- of steady state k)
    length_grid_k = 1000
    step1 = length_grid_k

    call initialize(length_grid_k)

    allocate(grid_k(length_grid_k))
    allocate(valuefn(length_grid_k,n_tauchen))
    allocate(g_k(length_grid_k,n_tauchen))
    allocate(g_c(length_grid_k,n_tauchen))

#ifdef HAS_LABOUR
    allocate(g_l(length_grid_k,n_tauchen))
    allocate(kend(length_grid_k,n_tauchen))
#endif

    call sub_grid_generation(grid_k, k_ss, cover, 2.D0)

    allocate(inputVars(3))
    allocate(outputVars(length_grid_k*n_tauchen))

    inputVars(1)=step1
    inputVars(2)=length_grid_k
    inputVars(3)=n_tauchen

    call initValueFn(inputVars, outputVars)

    forall(index_k=1:length_grid_k, i1=1:n_tauchen) valuefn(index_k,i1)=outputVars((index_k-1)*n_tauchen+i1)

    deallocate(inputVars)
    deallocate(outputVars)

#ifndef HAS_LABOUR
    !Here we call the subroutine to carry out the endogenous grid
    call sub_value(valuefn,g_k,g_c,grid_k,y,transition)
#else

    !Endogenous grid algortihm holding labor constant at the s.s.
    call sub_valueend1(valuefn,grid_k,length_grid_k,y,transition,l_ss)

    !First standard iteration to recover policy functions.
    call sub_valuestandard(valuefn,g_k,g_c,g_l,diff,grid_k,length_grid_k,y,transition)

    iter = 0
    do while ((diff>toler) .and. (iter<MAX_ITER))

        ! The first step is take the input from the standard algorithm and compute the values of Kend such that
        ! g_k(end)=k_grid
        do index_k = 1,length_grid_k
            do index_z = 1,n_tauchen
                call sub_golden(kend(index_k,index_z), grid_k(index_k), grid_k, g_k(:,index_z), length_grid_k)
            enddo
        enddo

        !Now do the endogenous grid for a few more extra steps
        valuefn=valuefn/(1-beta)
        call sub_valueend2(valuefn,grid_k,length_grid_k,y,transition,kend,g_l,diff)

        !One standard iteration to recover policy functions.
        call sub_valuestandard(valuefn,g_k,g_c,g_l,diff,grid_k,length_grid_k,y,transition)

        iter = iter+1
        if (mod(iter,1)==0) then
            print *, 'main Iteration: ', iter, 'Tolerance ', diff
            flush(6)
        end if

    enddo

#endif

    call CPU_TIME(t2)
    delta_t=t2+t0-2.0d0*t1
    print *, ' '
    print *, 'Elapsed time in the main program =', delta_t, 'sec'
    print *, ' '

    !----------------------------------------------------------------
    ! 4. Output
    !----------------------------------------------------------------

    open(unit=11,   file='value.txt', status = 'replace')
    write (11, '(5f20.10)') (valuefn(i1, :), i1 = 1,length_grid_k)
    close(11)

    open(unit=12,   file='g_c.txt', status = 'replace')
    write (12, '(5f20.10)') (g_c(i1, :), i1 = 1,length_grid_k)
    close(12)

    open(unit=13,   file='g_k.txt', status = 'replace')
    write (13, '(5f20.10)') (g_k(i1, :), i1 = 1,length_grid_k)
    close(13)

    !----------------------------------------------------------------
    ! 5. cleanup
    !----------------------------------------------------------------
    call clean()

    elapt = etime(ta)
    write (*,*), ' '
    write (*,'(a17, x, f8.2, x, a20)'), ' Program has used', elapt, 'seconds of CPU time,'
    write (*,'(a5, x, f8.2, x, a24, x, f5.2, x, a23)'), ' with', ta(1), 'seconds of user time and', ta(2), 'seconds of system time.'

    print *, ' '
    print *, 'End of the program'
    print *, ' '

end program main
