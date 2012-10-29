#define NEOCLASSICAL
!#define AIYAGARI
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

!
! This module puts a wrapper around a "root" function. This way, we can
! use it while still having it depend on a set of values we DONT want
! the solver to change.
! Note: I could do this with global variables, but they are always dangerous
!       This is sort-of like creating an object (for those who are familiar
!       with OO programming
!
! USAGE: The following main function provides an example of how to use this module
!
!        PROGRAM main
!
!            USE mywrapper
!
!            PROCEDURE(template_ValueFunction), POINTER :: func
!            PROCEDURE(template_RootFunction), POINTER :: rootFunc
!            REAL(8), DIMENSION(3) :: myKs, k_prime, myValue
!
!            myKs=(/ 1,2,3 /)
!            k_prime=(/ 1,2,3 /)
!            myValue=(/ 1,2,3 /)
!
!            func => valueFunction
!            rootFunc => rootFunction
!
!            call  init(myKs, myValue, func, rootFunc)
!            CALL execute_solve(k_prime)
!            call destroy()
!        CONTAINS
!
!            SUBROUTINE rootFunction(funcPointer,input,output)
!                implicit none
!
!                PROCEDURE(valueFunctionWrapper), POINTER, INTENT(in) :: funcPointer
!               REAL(8), DIMENSION(:), INTENT(inout) :: input
!                REAL(8), DIMENSION(:), INTENT(in) :: output
!
!                input = input+1
!            END SUBROUTINE rootFunction
!
!
!            FUNCTION valueFunction(n,k,k_prime) RESULT(y)
!                INTEGER, INTENT(in) :: n
!                REAL(8), DIMENSION(n), INTENT(in) :: k, k_prime
!                REAL, DIMENSION(n) :: y
!
!                y=k+k_prime
!            END FUNCTION valueFunction
!        END PROGRAM main
MODULE myRootWrapper
    use nrtype
    implicit none

    !depending on whether your root function solves for arrays or elements
    REAL(DP), allocatable ,DIMENSION(:) :: currentK, valueFun
    PROCEDURE(template_ValueFunction), POINTER :: fp_value
    PROCEDURE(template_RootFunction), POINTER :: fp_root
    INTEGER::sizeCheck !variable to make sure we don't check size every time we solve value function
    INTEGER::arraySize !variable to store size at initialization

    private :: currentK, valueFun
    private :: fp_value, fp_root
    private :: sizeCheck, arraySize

    ABSTRACT INTERFACE
        !
        ! An interface to your value function. It should have the same
        ! parameter declarations as the actual function
        ! to prevent potential memory leaks, need to pass in the size of the array
        FUNCTION template_ValueFunction(n,k,k_prime) RESULT(y)
            INTEGER, INTENT(in) :: n
            REAL(8), DIMENSION(n), INTENT(in) :: k, k_prime
            REAL, DIMENSION(n) :: y
        END FUNCTION template_ValueFunction
        !
        ! An interface to your solver (i.e. root) function. It should have the same
        ! parameter declarations as the actual function
        ! Given that I don't know how the one you are using works, here I assume
        ! that it takes a function pointer, the input variables, and the output variables.
        ! it modifies the input variables so that the result of the function pointer is
        ! the output variables
        SUBROUTINE template_RootFunction(funcPointer,input,output)
            PROCEDURE(valueFunctionWrapper), POINTER, INTENT(in) :: funcPointer
            REAL(8), DIMENSION(:), INTENT(inout) :: input
            REAL(8), DIMENSION(:), INTENT(in) :: output
        END SUBROUTINE template_RootFunction

    END INTERFACE

CONTAINS

    SUBROUTINE init(myKs, myValue, fp_valueFun, fp_rootFun)
        !initialization of variables
        REAL(DP), DIMENSION(:), INTENT(in) :: myKs, myValue
        PROCEDURE(template_ValueFunction), POINTER, INTENT(in) :: fp_valueFun
        PROCEDURE(template_RootFunction), POINTER, INTENT(in) :: fp_rootFun

         !error check
        arraySize = size(myKs,dim=1)
        if(arraySize /= size(myValue,dim=1)) then
            stop 'value and capital array do not match in size'
        end if

        allocate(currentK(arraySize))
        allocate(valueFun(arraySize))

        currentK = myKs
        valueFun = myValue

        fp_value => fp_valueFun
        fp_root =>  fp_rootFun
    END SUBROUTINE init

    SUBROUTINE destroy()
        !deallocate variables
        sizeCheck=0
        deallocate(currentK)
        deallocate(valueFun)
    END SUBROUTINE destroy

    FUNCTION valueFunctionWrapper(k_prime) RESULT(v)
        REAL(DP), DIMENSION(arraySize), INTENT(in) :: k_prime
        REAL(DP), DIMENSION(arraySize) :: v
        !error check
        if(sizeCheck /= 1) THEN
            sizeCheck = 1
            if(arraySize /= size(k_prime,dim=1)) then
                stop 'valueFunctionWrapper: value and next period capital array do not match in size'
            end if
        end if

        v = fp_value(arraySize,currentK,k_prime)
    END FUNCTION valueFunctionWrapper

    SUBROUTINE execute_solve(k_prime)
        REAL(DP), DIMENSION(arraySize), INTENT(inout) :: k_prime

        PROCEDURE(valueFunctionWrapper), POINTER :: fp_vfw

        fp_vfw => valueFunctionWrapper

        call fp_root(fp_vfw, k_prime, valueFun)
    END SUBROUTINE execute_solve

END MODULE myRootWrapper



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

    real (DP), allocatable, dimension (:) :: shocks                 ! the actual shock values
    real(DP), allocatable, dimension(:,:) :: transitionMatrix       ! The markov matrix representing probability of moving
                                                                    ! from one state to another
    real (DP), allocatable, dimension (:) :: grid_k                 ! the capital values on which we want to obtain the estimate
    real (DP), allocatable, dimension (:) :: stepVec                ! distance between subsequent capital values
    integer :: shockCount, capitalCount                             ! the size of the capital array and shock array
#ifdef AIYAGARI
    real(DP) :: w, r
    private:: w,r
#endif

    private::alpha, beta, tau, delta, theta, rho, sigma, transitionMatrix, stepVec, grid_k
    private:: k_ss, l_ss, c_ss, r_ss, y_ss

contains

    subroutine sub_modelclean()
        deallocate(stepVec)
        deallocate(transitionMatrix)
        deallocate(shocks)
        deallocate(grid_k)
    end subroutine sub_modelclean

#ifdef NEOCLASSICAL
    subroutine sub_modelInitialize(steps,states)
#else
#ifdef AIYAGARI
    subroutine sub_modelInitialize(steps,states,wage)
        integer, intent(in) :: wage
#else
    subroutine sub_modelInitialize(steps,states)
#endif
#endif
        integer, intent(in) :: states
        integer, intent(in) :: steps

        shockCount = states
        capitalCount = steps
        allocate(stepVec(steps))
        allocate(grid_k(steps))
        allocate(shocks(states))
        allocate(transitionMatrix(states,states))
#ifdef AIYAGARI
        r=(1-.001)/beta-1
        w=wage
#endif
    end subroutine sub_modelInitialize

    subroutine sub_modelSetStates(capital_grid,shocks_vector,transitionFunction)
        REAL(DP), dimension(:), intent(in) :: capital_grid
        REAL(DP), dimension(:), intent(in) :: shocks_vector
        REAL(DP), dimension(:,:), intent(in) :: transitionFunction
        integer :: n,m,index_k

        n=size(grid_k,dim=1)
        IF (size(capital_grid,dim=1)/=n) THEN
            PRINT '(a,i3,a,i3)', 'sub_modelSetStates: capital_grid must be a matrix of size ',n
            call sub_modelstop('program terminated by sub_modelSetStates')
        END IF

        m=size(transitionMatrix,dim=1)
        IF (size(shocks_vector,dim=1)/=m) THEN
            PRINT '(a,i3,a,i3)', 'sub_modelSetStates: shocks_vector must be a matrix of size ',m
            call sub_modelstop('program terminated by sub_modelSetStates')
        END IF
        IF ((size(transitionFunction,dim=1)/=m) .or. (size(transitionFunction,dim=1)/=m)) THEN
            PRINT '(a,i3,a,i3)', 'sub_modelSetStates: transitionFunction must be a matrix of size ',n,'-by-',m
            call sub_modelstop('program terminated by sub_modelSetStates')
        END IF

        grid_k(:)=capital_grid(:)
        do index_k = 2,n
            stepVec(index_k-1)=grid_k(index_k)-grid_k(index_k-1)
        end do

        shocks(:)=shocks_vector(:)

        do index_k = 1,m
            transitionMatrix(index_k,:)=transitionFunction(index_k,:)
        end do
    end subroutine sub_modelSetStates


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
        !   2) k - the capital/asset level at which to evaluate the output
        !   3) l - the l at which to evaluate the output
        real(DP) :: sub_modelOutput
        real(DP), dimension(:), intent(IN) :: inputVars

#ifdef AIYAGARI
        ! in AIYAGARI, we don't have a production function. Rather, we have a first order condition that
        !               c+a'= (1+r)a+sw
        ! where r is the market clearing interest rate, a is current asset levels, s is the idiosyncratic shock, and w
        ! is the wage.
        !
        !The RHS of this equation can be considered the output function.
        real(DP) :: z, a

        z=inputVars(1)
        a=inputVars(2)

        sub_modelOutput = (1+r)*a+z*w
#endif

#ifdef NEOCLASSICAL
        real(DP) :: z, k

        z=inputVars(1)
        k=inputVars(2)

        sub_modelOutput = exp(z)*k**alpha+(1-delta)*k
#endif
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

#ifdef NEOCLASSICAL
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

    subroutine initFn(inputVars, outputVars)
            ! This subroutine contains the calculations to initialize the initial guess of an increasing
            ! function.
            !
            ! INPUTS: inputVars - a vector of input values. typicaly
            !                         1) step1:  the step size
            !                         2) m: the number of points on the k grid
            ! OUTPUTS: outputVars - a vector (which is actually a linearized matrix) representing the initialized function
            !
            ! Note: To pass arrays, see the initValueFn subroutine below
        implicit none

        real(DP), dimension(:), INTENT(IN) :: inputVars
        real(DP), dimension(:), INTENT(OUT) :: outputVars

#ifdef NEOCLASSICAL
        real(DP) :: step1, valueinitial
        integer :: index_k, i1, m, n

        step1 = inputVars(1)
        m = inputVars(2)
        n = inputVars(3)

        valueinitial = (1/(1-beta))*c_ss**(1-tau)/(1-tau)

        forall (index_k = 1:m, i1=1:n) outputVars((index_k-1)*n+i1) = index_k/step1 + valueinitial
#endif

#ifdef AIYAGARI
        step1 = inputVars(1)
        m = inputVars(2)
        n = inputVars(3)

        valueinitial = (1/(1-beta))*c_ss**(1-tau)/(1-tau)

        forall (index_k = 1:m, i1=1:n) outputVars((index_k-1)*n+i1) = 0.

#endif
    end subroutine initFn

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

        real(DP) :: step1, valueinitial
        integer :: index_k, i1, m, n

        step1 = inputVars(1)
        m = inputVars(2)
        n = inputVars(3)

        valueinitial = (1/(1-beta))*c_ss**(1-tau)/(1-tau)

        forall (index_k = 1:m, i1=1:n) outputVars((index_k-1)*n+i1) = index_k/step1 + valueinitial
    end subroutine initValueFn


    SUBROUTINE getNextGuess(fnToSolve, myGuess, onCapital, targetCapital)
            ! This subroutine gives the next guess of the function we are endogenizing
            !
            ! INPUTS: fnToSolve - the current values of the array
            ! OUTPUTS:
            !           myGuess - the guess of the next value
            !           onCapital - the capital levels for which the guess is evaluating
            !           targetCapital - the capital levels for which we want to evaluate according to model initialization
        real(DP), dimension(:,:), INTENT(IN) :: fnToSolve
        real(DP), dimension(:,:), INTENT(OUT) :: myGuess
        real(DP), dimension(:,:), INTENT(OUT) :: onCapital
        real(DP), dimension(:,:), INTENT(OUT) :: targetCapital

        integer :: m,n, index_z

        real(DP), allocatable, dimension(:,:) :: cih, Cstar

        m=size(grid_k)
        n=size(shocks)

        allocate(cih(m,n))
        allocate(Cstar(m,n))

        ! define the cash in hand array
        call cashInHand(cih, m,n)

        !Get optimal consumption
        !
        Cstar = getCStar(m,n,fnToSolve)

        !get optimal value function
        myGuess=getVStar(m,n,Cstar,fnToSolve)

        do index_z = 1,n
            onCapital(:,index_z)=Cstar(:,index_z)+grid_k
            targetCapital(:,index_z)=cih(:,index_z)
        enddo

        deallocate(cih)
        deallocate(Cstar)

    end SUBROUTINE getNextGuess

    SUBROUTINE cashInHand(outputVars, m, n)
            !
            ! This subroutine defines tomorrow's "cash in hand" using model parameters.
            !
            ! INPUTS: m: the number of points on the k grid
            !         n: the number of exogenous states
            ! OUTPUTS: outputVars - a m-by-n matrix representing the cash in hand function
        integer, INTENT(IN) :: m, n
        real(DP), dimension(:,:), INTENT(OUT) :: outputVars
        real(DP), allocatable, dimension(:) :: newInput
        integer :: index_k, index_z

        IF (size(grid_k,dim=1)/=m) THEN
            PRINT '(a,i3)', 'cashInHand: grid_k must be a vector of size ',m
            call sub_modelstop('program terminated in cash in hand')
        END IF

        IF (size(shocks,dim=1)/=n) THEN
            PRINT '(a,i3)', 'cashInHand: shocks must be a vector of size ',n
            call sub_modelstop('program terminated in cash in hand')
        END IF
        IF ( (size(outputVars,dim=1)/=m) .and. (size(outputVars,dim=2)/=n) ) THEN
            PRINT '(a,i3)', 'cashInHand: outputVars must be a matrix of size ',m, '-by-',n
            call sub_modelstop('program terminated in cash in hand')
        END IF

        allocate(newInput(2))
        do index_k = 1,m
            newInput(2)=grid_k(index_k)
            do index_z = 1,n
                newInput(1)=shocks(index_z)
                outputVars(index_k,index_z)=sub_modelOutput(newInput)
            end do
        end do
        deallocate(newInput)

    END SUBROUTINE cashInHand

    FUNCTION getCStar(length_grid_k,n_tauchen,fnToSolve) RESULT(y)
        ! Get optimal consumption given the value function
        ! Note: The last step is from the envelope condition
        INTEGER, INTENT(IN) :: length_grid_k
        INTEGER, INTENT(IN) :: n_tauchen
        REAL(DP), DIMENSION(length_grid_k,n_tauchen), INTENT(IN) :: fnToSolve
        REAL(DP), DIMENSION(length_grid_k,n_tauchen) :: y
        REAL(DP), DIMENSION(length_grid_k,n_tauchen) :: D
        INTEGER :: index_k

#ifdef NEOCLASSICAL
        !use slope between points as an estimate of the derivative at that point
        D(1,:)=(fnToSolve(2,:)-fnToSolve(1,:))/stepVec(1)
        D(length_grid_k,:)=(fnToSolve(length_grid_k,:)-fnToSolve(length_grid_k-1,:))/stepVec(length_grid_k-1)
        do index_k = 2,length_grid_k-1
            D(index_k,:)=(fnToSolve(index_k+1,:)-fnToSolve(index_k-1,:))/(stepVec(index_k)+stepVec(index_k-1))
        enddo
        y=D**(-1/tau)
#endif
#ifdef AIYAGARI
    !In Aiyagari, we can directly use the euler conditition.
        REAL(DP) :: temp
        do index_k=1,length_grid_k
            D(index_k,:)=1/(1+r)*((1+r)*beta*())+
#endif
    END FUNCTION getCStar

    FUNCTION getVStar(length_grid_k,n_tauchen,Cstar,valuefn) RESULT(y)
        INTEGER, INTENT(IN) :: length_grid_k
        INTEGER, INTENT(IN) :: n_tauchen
        REAL(DP), DIMENSION(length_grid_k,n_tauchen), INTENT(IN) :: valuefn, Cstar
        REAL(DP), DIMENSION(length_grid_k,n_tauchen) :: y

        y=1/(1-tau)*Cstar**(1-tau)+valuefn
    END FUNCTION

    FUNCTION sub_evaluate(valuefn) RESULT(y)
        REAL(DP), DIMENSION(capitalCount,shockCount), INTENT(IN) :: valuefn
        REAL(DP), DIMENSION(capitalCount,shockCount) :: y

        integer :: index_k, index_z

            do index_k = 1,capitalCount
                do index_z = 1,shockCount
                    y(index_k,index_z)=beta*dot_product(transitionMatrix(index_z,:),valuefn(index_k,:))
                enddo
            enddo
    end FUNCTION sub_evaluate


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
    real (DP), allocatable, dimension (:) :: stepVec            ! possible exogenous shocks
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

    SUBROUTINE sub_grid_generation(x,xcentre,xbounds,s, gentype)
        ! Purpose: Generate grid x on [xcentre*(1-xbounds),xcentre*(1+xbounds)] using spacing parameter s set as follows:
        ! s=1       linear spacing
        ! s>1       left skewed grid spacing with power s
        ! 0<s<1     right skewed grid spacing with power s
        ! s<0       geometric spacing with distances changing by a factor -s^(1/(n-1)), (>1 grow, <1 shrink)
        ! s=-1      logarithmic spacing with distances changing by a factor (xmax-xmin+1)^(1/(n-1))
        ! s=0       logarithmic spacing with distances changing by a factor (xmax/xmin)^(1/(n-1)), only if xmax,xmin>0
        !
        ! INPUTS: Inputs depend on gentype. If gentype is missing:
        !               xcentre - The centre of the grid
        !               xbounds - how far away from the centre the bounds go.
        !               s - skewness of grid - see above
        !         If gentype is provided and equals 1
        !               xcentre - The min point of the grid
        !               xbounds - The max point of the grid
        !               s - skewness of grid - see above
        ! OUTPUT: x - the generated grid

        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(OUT) :: x
        REAL(DP), INTENT(IN) :: xcentre,xbounds, s
        INTEGER, OPTIONAL, INTENT(IN) :: gentype
        REAL(DP) :: c ! growth rate of grid subintervals for logarithmic spacing
        REAL(DP) :: xmax, xmin
        INTEGER :: n,i

        if(present(gentype)) then
            if(gentype==1) then
                xmin=xcentre
                xmax=xbounds
            else
                PRINT *, 'grid_generation: unsuported gentype: ',gentype
                call sub_mystop('program terminated by grid_generation')
            endif
        else
            xmax=xcentre*(1+xbounds);
            xmin=xcentre*(1-xbounds);
        endif
        n=size(x)

    print *, 'Beginning forall, n=',n
    print *, ' '
    flush(6)

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
        ! INPUTS: x: Indexes for which f_x has values
        !         f_x: values for each x value in "x"
        !         xp: x values at which we want to find f_x (which we do through linear interpolation)
        !         length_x: the number of elements in x
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

    !a preferred algorithm
    ! given initial guess of value function
    ! get new guess (inputs: initial guess, states; outputs: values, index for those values)
    ! interpolate between new guess and the proper initial values
    ! get new value function estimate
    subroutine sub_endogenize(fnToSolve, grid_k, y, transition)
            ! The generalized endogenous grid method
            ! INPUTS: fnToSolve - The initial guess of the value function, size m-by-n
            !         grid_k - grid of possible capital values (..,ks,..), size m
            !         y - the possible shocks - size n
            !         transition - the transition matrix - size n-by-n
            ! OUTPUT: valuefn - the final value function - size m-by-n

        implicit none

        real(DP), dimension(:,:), intent(inout) :: fnToSolve
        real(DP), dimension(:), intent(in) :: grid_k
        real(DP), dimension(:), intent(in) :: y
        real(DP), dimension(:,:), intent(in) :: transition

        integer :: iter, index_k, index_z
        integer :: length_grid_k,m,n

        real(DP), allocatable, dimension(:,:) :: oncapital,newguess, targetcapital, Vend1, value1, kend
        real(DP) :: diff

        length_grid_k=size(grid_k)
        m=length_grid_k
        n=n_tauchen

        !
        !allocate arrays
        !
        allocate(oncapital(m,n))
        allocate(newguess(m,n))
        allocate(targetcapital(m,n))
        allocate(Vend1(m,n))
        allocate(value1(m,n))
        allocate(kend(m,n))

        iter = 0
        diff = 1000.d0

        call sub_modelSetStates(grid_k,y,transition)

        !Here we start the iterations.
        print *,"starting iterations"
        flush(6)
        do while((diff>toler) .and. (iter<MAX_ITER))

            iter=iter+1

            !cstar is the next guess but on different capital
            !Yend is the capital levels for which it is defined
            call getNextGuess(fnToSolve,newguess,oncapital,targetcapital)

            do index_z = 1,n
                call sub_myinterp1(oncapital(:,index_z),newguess(:,index_z), targetcapital, length_grid_k, Vend1(:,index_z))
            enddo

            value1=sub_evaluate(Vend1)
            diff = maxval(abs(value1-fnToSolve))
            fnToSolve=value1
            if (mod(iter,50)==0) then
                print *, 'subvalue: Iteration: ', iter, 'Tolerance ', diff
                flush(6)
            end if
        enddo

        ! At this point the program recovers the endogenous grid for capital.
        do index_k = 1,length_grid_k
            do index_z = 1,n_tauchen
                call sub_kendogenousnewton(kend(index_k,index_z),oncapital(index_k,index_z),y(index_z),grid_k(index_k))
            enddo
        enddo

        ! Here we interpolate capital on the grid so that we can compare the results to the standard algorithm.
!        do index_z = 1,n_tauchen
!            call sub_myinterp1(kend(:,index_z),Vend(:,index_z), grid_k, length_grid_k, fnToSolve(:,index_z))
!            call sub_myinterp1(kend(:,index_z),grid_k, grid_k, length_grid_k, g_k(:,index_z))
!            call sub_myinterp1(kend(:,index_z),Cstar(:,index_z), grid_k, length_grid_k, g_c(:,index_z))
!        enddo

        !
        !deallocate arrays
        !
        deallocate(oncapital)
        deallocate(newguess)
        deallocate(targetcapital)
        deallocate(Vend1)
        deallocate(value1)
        deallocate(kend)
    end subroutine sub_endogenize

    subroutine sub_value(fnToSolve, g_k, g_c, grid_k, y, transition)
            ! The actual value iteration function
            ! INPUTS: fnToSolve - The initial guess of the value function, size m-by-n
            !         grid_k - grid of possible capital values (..,ks,..), size m
            !         y - the possible shocks - size n
            !         transition - the transition matrix - size n-by-n
            ! OUTPUT: g_c - interpolated value of consumption - size m-by-n
            !         g_k - interpolated value of capital - size m-by-n
            !         valuefn - the final value function - size m-by-n

        implicit none

        real(DP), dimension(:,:), intent(inout) :: fnToSolve
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
        IF ((size(fnToSolve,dim=1)/=m) .or. (size(fnToSolve,dim=2)/=n)) THEN
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
        call sub_modelSetStates(grid_k,y,transition)
        call cashInHand(cih, length_grid_k, n_tauchen)

        !Here we start the iterations.
        print *,"starting iterations"
        flush(6)
        do while((diff>toler) .and. (iter<MAX_ITER))

            iter=iter+1

            !Get optimal consumption
            !
            Cstar = getCStar(length_grid_k,n_tauchen,fnToSolve)

            !get optimal value function
            Vend=getVStar(length_grid_k,n_tauchen,Cstar,fnToSolve)

            do index_z = 1,n_tauchen
                Yend(:,index_z)=Cstar(:,index_z)+grid_k
                call sub_myinterp1(Yend(:,index_z),Vend(:,index_z), cih(:,index_z), length_grid_k, Vend1(:,index_z))
            enddo

            ! This is sort of model specific
            do index_k = 1,length_grid_k
                do index_z = 1,n_tauchen
                    value1(index_k,index_z)=beta*dot_product(transition(index_z,:),Vend1(index_k,:))
                enddo
            enddo

            diff = maxval(abs(value1-fnToSolve))
            fnToSolve=value1
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
            call sub_myinterp1(kend(:,index_z),Vend(:,index_z), grid_k, length_grid_k, fnToSolve(:,index_z))
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

    call CPU_TIME(t0)
    call CPU_TIME(t1)

    print *, 'Beginning of the Program Value Function Iteration'
    print *, ' '
    flush(6)

    !----------------------------------------------------------------
    ! 1. Tauchen Points
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

    cover = 0.25d0         ! Coverage of the grid (+- of steady state k)
    length_grid_k = 1000
    step1 = length_grid_k

    !----------------------------------------------------------------
    ! 1.5 Initialize modules
    !----------------------------------------------------------------

    !    call initialize(length_grid_k)
    call sub_modelInitialize(length_grid_k,n_tauchen)

    !----------------------------------------------------------------
    ! 2. Computation of the Steady State
    !----------------------------------------------------------------

    allocate(outputVars(3))

    call steadyState(outputVars)

    y_ss = outputVars(1)
    c_ss = outputVars(2)
    k_ss = outputVars(3)

    deallocate(outputVars)

    !----------------------------------------------------------------
    ! 3. Endogenous Grid algorithm
    !----------------------------------------------------------------

    print *, 'Beginning the iteration of the value function with an initial linear guess '
    print *, ' '
    flush(6)

    allocate(grid_k(length_grid_k))
    allocate(valuefn(length_grid_k,n_tauchen))
    allocate(g_k(length_grid_k,n_tauchen))
    allocate(g_c(length_grid_k,n_tauchen))

    print *, 'Beginning grid gen'
    print *, ' '
    flush(6)
    call sub_grid_generation(grid_k, k_ss, cover, 2.D0)

    allocate(inputVars(3))
    allocate(outputVars(length_grid_k*n_tauchen))

    inputVars(1)=step1
    inputVars(2)=length_grid_k
    inputVars(3)=n_tauchen

    print *, 'Beginning init'
    print *, ' '
    flush(6)
    call initFn(inputVars, outputVars)
    call initValueFn(inputVars, outputVars)

    forall(index_k=1:length_grid_k, i1=1:n_tauchen) valuefn(index_k,i1)=outputVars((index_k-1)*n_tauchen+i1)

    deallocate(inputVars)
    deallocate(outputVars)

    !Here we call the subroutine to carry out the endogenous grid
    print *, 'Beginning endogenize'
    print *, ' '
    flush(6)
    call sub_endogenize(valuefn,grid_k,y,transition)
!    call sub_value(valuefn,g_k,g_c,grid_k,y,transition)

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
