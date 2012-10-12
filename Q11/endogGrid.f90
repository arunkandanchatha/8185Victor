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

! This module contains everything we neeed to perform a value function interation using the endogenous
! grid method. The main program below gives an example of usage.
!
! USAGE:
MODULE  endogenousGrid
    use nrtype
    implicit none


    !----------------------------------------------------------------
    ! 1. Calibration
    !----------------------------------------------------------------

    real (DP), parameter :: delta  = 0.0196                         ! Depreciation
    real (DP), parameter :: alpha  = 0.4                            ! Capital Share
    real (DP), parameter :: tau    = 2.d0                           ! risk aversion
    real (DP), parameter :: beta   = 0.9896                         ! Discount factor
    real (DP), parameter :: rho    = 0.95                           ! Autoregressive
    real (DP), parameter :: sigma  = 0.007                          ! Variance

    !----------------------------------------------------------------
    ! 2. Numerical parameters
    !----------------------------------------------------------------

    integer :: n_tauchen = 41                                   ! Number Tauchen points
    integer :: Tsimul    = 50000                                ! Number of simulations for the EEerror

    real (DP), parameter :: toler   = 1e-6                      ! Numerical tolerance
    real (DP), allocatable, dimension (:) :: stepVec
    integer :: i1

    private :: sub_myinterp1, sub_kendogenousnewton
    private :: i1, toler

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
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(OUT) :: x
        REAL(DP), INTENT(IN) :: xcentre,xbounds, s
        REAL(DP) :: c ! growth rate of grid subintervals for logarithmic spacing
        REAL(DP) :: xmax, xmin
        INTEGER :: n,i

        n=size(x)
        xmax=xcentre*(1+xbounds);
        xmin=xcentre*(1-xbounds);

        print *, "New grid"
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

    subroutine sub_kendogenousnewton(kend,Yend,z,kguess)

        implicit none

        real(DP), intent(in) :: Yend,z,kguess
        real(DP), intent(out) :: kend

        integer :: go_on

        real(DP) :: f_k, f_kp, kp

        kend=kguess

        go_on = 1
        do while (go_on == 1)

            !ugh. Again, specific to the function. We should change this as well
            f_k  = Yend-exp(z)*kend**alpha-(1-delta)*kend
            f_kp = -alpha*exp(z)*kend**(alpha-1)-(1-delta)
            kp = kend-(f_k/f_kp)

            if (abs(kp-kend)<0.0000001) go_on = 0
            kend = kp
        enddo

    end subroutine sub_kendogenousnewton

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

        !
        ! Need to define "cash in hand" tomorrow.
        !
        forall (index_k = 1:length_grid_k, index_z=1:n_tauchen)
            cih(index_k,index_z)=exp(y(index_z))*grid_k(index_k)**alpha+(1-delta)*grid_k(index_k)
        end forall

            !Here we start the iterations.

        print *,"starting iterations"
        flush(6)
        do while((diff>toler) .and. (iter<5000))

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
            Cstar=D**(-1/tau)

            !again, this is for a very specific value function. Any way  we could
            ! pass in a function pointer?
            Vend=1/(1-tau)*Cstar**(1-tau)+valuefn

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
                print *, 'Iteration: ', iter, 'Tolerance ', diff
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
END MODULE endogenousGrid

program main
    use nrtype
    use endogenousGrid
    implicit none

    real(DP), allocatable, dimension(:,:)  :: transition, transitioncdf, valuefn, g_k, g_c
    real(DP), allocatable, dimension(:)  :: y, grid_k

    integer :: index_k, length_grid_k, i1

    real(DP) :: k_ss, c_ss, y_ss
    real(DP) :: cover, cover_tauchen, valueinitial, step1

    real (4) :: elapt, ta(2) !This is used in order to time the main block of the program.
    real(dp) :: t0,t1,t2,delta_t

    call CPU_TIME(t0)
    call CPU_TIME(t1)

    print *, 'Beginning of the Program Value Function Iteration'
    print *, ' '

    !----------------------------------------------------------------
    ! 1. Computation of the Steady State
    !----------------------------------------------------------------

    k_ss = ((1/beta-1+delta)/alpha)**(1/(alpha-1))
    y_ss = k_ss**alpha
    c_ss = k_ss**alpha-delta*k_ss

    print *, 'Steady State values'
    print *, 'Output: ', y_ss, 'Capital: ', k_ss

    allocate(transition(n_tauchen,n_tauchen))
    allocate(transitioncdf(n_tauchen,n_tauchen))
    allocate(y(n_tauchen))

    !----------------------------------------------------------------
    ! 2. Tauchen Points
    !----------------------------------------------------------------

    cover_tauchen = 3.d0
    call sub_tauchen(y,transition,transitioncdf,rho,sigma,cover_tauchen)

    print *, ' '
    print *, 'For the stochastic shock, Tauchen`s approximation method with', n_tauchen , 'points will be used.'
    print *, ' '
    flush(6)

    !----------------------------------------------------------------
    ! 3. Endogenous Grid algorithm
    !----------------------------------------------------------------

    cover = 0.25d0         ! Coverage of the grid (+- of steady state k)
    length_grid_k = 1000
    step1 = length_grid_k


    call initialize(length_grid_k)

    allocate(grid_k(length_grid_k))
    allocate(valuefn(length_grid_k,n_tauchen))
    allocate(g_k(length_grid_k,n_tauchen))
    allocate(g_c(length_grid_k,n_tauchen))

    call sub_grid_generation(grid_k, k_ss, cover, 2.D0)

    !
    ! This is the initial increasing guess for the value function.
    !
    valueinitial = (1/(1-beta))*c_ss**(1-tau)/(1-tau)

    do index_k = 1,length_grid_k
        valuefn(index_k,:) = index_k/step1 + valueinitial
    enddo

    !Here we call the subroutine to carry out the endogenous grid
    call sub_value(valuefn,g_k,g_c,grid_k,y,transition)

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
