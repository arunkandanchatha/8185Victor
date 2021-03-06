!#define NEOCLASSICAL
#define AIYAGARI
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
!    real (DP), parameter :: beta   = 0.9896                         ! Discount factor
    real (DP), parameter :: beta   = 0.96                         ! Discount factor
    real (DP), parameter :: theta   = 0.357                         ! weight on consumption
    real (DP), parameter :: rho    = 0.95                           ! Autoregressive
    real (DP), parameter :: sigma  = 0.007                          ! Variance
#ifdef AIYAGARI
    real(DP), parameter :: r = .035
#endif

    real(DP) :: k_ss, l_ss, c_ss, r_ss, y_ss

    real (DP), allocatable, dimension (:) :: shocks                 ! the actual shock values
    real(DP), allocatable, dimension(:,:) :: transitionMatrix       ! The markov matrix representing probability of moving
                                                                    ! from one state to another
    real (DP), allocatable, dimension (:) :: grid_k                 ! the capital values on which we want to obtain the estimate
    real (DP), allocatable, dimension (:) :: stepVec                ! distance between subsequent capital values
    integer :: shockCount, capitalCount                             ! the size of the capital array and shock array
#ifdef AIYAGARI
    real(DP) :: w
    private:: w,r
#endif

    private::alpha, beta, tau, delta, theta, rho, sigma, transitionMatrix, stepVec, grid_k
    private:: k_ss, l_ss, c_ss, r_ss, y_ss

#ifdef NEOCLASSICAL
    private::cashInHand,getCStar,getVStar
#endif

contains

    subroutine sub_modelclean()
        deallocate(stepVec)
        deallocate(transitionMatrix)
        deallocate(shocks)
        deallocate(grid_k)
    end subroutine sub_modelclean

#ifdef AIYAGARI
    subroutine sub_modelInitialize(steps,states,wage)
        integer, intent(in) :: wage
#else
        subroutine sub_modelInitialize(steps,states)
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
        w=wage
#endif
    end subroutine sub_modelInitialize

    subroutine sub_modelSetStates(capital_grid,shocks_vector,transitionFunction)
        REAL(DP), dimension(:), intent(in) :: capital_grid
        REAL(DP), dimension(:), intent(in) :: shocks_vector
        REAL(DP), dimension(:,:), intent(in) :: transitionFunction
        integer :: n,m,index_k

        n=capitalCount
        IF (size(capital_grid,dim=1)/=n) THEN
            PRINT '(a,i3,a,i3)', 'sub_modelSetStates: capital_grid must be a matrix of size ',n
            call sub_modelstop('program terminated by sub_modelSetStates')
        END IF

        m=shockCount
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
        flush(6)
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
#else
#ifdef NEOCLASSICAL
        real(DP) :: z, k

        z=inputVars(1)
        k=inputVars(2)

        sub_modelOutput = exp(z)*k**alpha+(1-delta)*k
#else
        call sub_modelstop("Undefined model called in sub_modelOutput")
#endif
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

#ifdef AIYAGARI
        real(DP) :: z
        z=inputVars(1)
        sub_modelf_kp = -1.0
#else
#ifdef NEOCLASSICAL
        real(DP) :: z, k
        z=inputVars(1)
        k=inputVars(2)

        sub_modelf_kp = -alpha*exp(z)*k**(alpha-1)-(1-delta)
#else
        call sub_modelstop("Undefined model called in sub_modelf_kp")
#endif
#endif
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

#ifndef NEOCLASSICAL
        call sub_modelstop("Undefined model called in sub_modelf_lp")
#endif
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
#else
        call sub_modelstop("Undefined model calling steadystate")
        outputVars(1) = 0
#endif
    end subroutine steadystate

    subroutine initGuess(valueFn)
            ! This subroutine contains the calculations to initialize the initial guess of an increasing
            ! function. For example, if neoclassical, then the value function. If Aiyagari, then the asset function
            !
            ! OUTPUTS: outputVars - a matrix representing the initialized function
            !
        implicit none

        real(DP), dimension(capitalCount,shockCount), INTENT(OUT) :: valueFn

#ifdef NEOCLASSICAL
        real(DP) :: valueinitial
        integer :: index_k, i1

        step1 = capitalCount
        valueinitial = (1/(1-beta))*c_ss**(1-tau)/(1-tau)

        forall (index_k = 1:capitalCount, i1=1:shockCount) valueFn(index_k,i1) = index_k/step1 + valueinitial
#else
#ifdef AIYAGARI
        integer :: index_k, i1

        forall (index_k = 1:capitalCount, i1=1:shockCount) valueFn(index_k,i1) = 0.
#else
        call sub_modelstop("Undefined model calling initGuess")
#endif
#endif
    end subroutine initGuess

    SUBROUTINE getNextGuess(fnToSolve, myGuess, onCapital, targetCapital)
            ! This subroutine gives the next guess of the function we are endogenizing
            !
            ! INPUTS: fnToSolve - the current values of the array
            ! OUTPUTS:
            !           myGuess - the guess of the next value
            !           onCapital - the capital levels for which the guess is evaluating
            !           targetCapital - the capital levels for which we want to evaluate according to model initialization
        real(DP), dimension(capitalCount,shockCount), INTENT(IN) :: fnToSolve
        real(DP), dimension(capitalCount,shockCount), INTENT(OUT) :: myGuess, onCapital,targetCapital

#ifdef NEOCLASSICAL
        integer :: index_z
        real(DP), dimension(capitalCount,shockCount) :: cih, Cstar

        ! define the cash in hand array
        call cashInHand(cih)

        !Get optimal consumption
        !
        Cstar = getCStar(fnToSolve)

        !get optimal value function
        myGuess=getVStar(Cstar,fnToSolve)

        do index_z = 1,shockCount
            onCapital(:,index_z)=Cstar(:,index_z)+grid_k
            targetCapital(:,index_z)=cih(:,index_z)
        enddo
#else
#ifdef AIYAGARI
        ! We actually need to perform a little trick. The fnToSolve is a'', not a' (which is
        ! stored in grid_k. We will use a' to find a, then return the following values:
        !   myGuess: The theoretical a', which we will keep as the original a', i.e. grid_k
        !   onCapital: This is the "a" given a'.
        !   targetCapital: a' is the next policy function given a, and we want to then find the
        !                  value of the policy function on a'. So this would be grid_k as well
        !
        !Note: sub_evaluate() will need to know this to properly return the next iteration of the
        !      endogenous function
        integer :: i,j
        real(DP), dimension(capitalCount,shockCount) :: eul_temp

        do i=1,capitalCount
            do j=1,shockCount
                eul_temp(i,j) = ((1.0_DP+r)*grid_k(i)+shocks(j)*w-fnToSolve(i,j))**(-tau)
            end do
        end do

        eul_temp = transpose((1.0_DP+r)*beta*matmul(transitionMatrix,transpose(eul_temp)))
        do i=1,capitalCount
            do j=1,shockCount
                oncapital(i,j)=(grid_k(i)-shocks(j)*w+eul_temp(i,j)**(-1/tau))/(1.0_DP+r)
                myGuess(i,j)=grid_k(i)
                targetCapital(i,j)=grid_k(i)
            end do
        end do
#else
        call sub_modelstop("Undefined model calling getNextGuess")
#endif
#endif
    end SUBROUTINE getNextGuess

    FUNCTION sub_evaluate(valuefn) RESULT(y)
            ! This subroutine assumes that you have already called getNextGuess. It then takes
            ! the interpolated value returned by that function, and finds the next
            ! guess guess of the function we are endogenizing
            !
            ! INPUTS: valuefn - the current values of the array
            ! OUTPUTS: y - the guess of the next value
        REAL(DP), DIMENSION(capitalCount,shockCount), INTENT(IN) :: valuefn
        REAL(DP), DIMENSION(capitalCount,shockCount) :: y

#ifdef NEOCLASSICAL
        integer :: index_k, index_z
        do index_k = 1,capitalCount
            do index_z = 1,shockCount
                y(index_k,index_z)=beta*dot_product(transitionMatrix(index_z,:),valuefn(index_k,:))
            enddo
        enddo
#else
#ifdef AIYAGARI
            !In Aiyagari, valuefn is a''
        integer :: i,j
        !borrowing constraint
        y=valueFn
        do i=1,capitalCount
            do j=1,shockCount
                if (valuefn(i,j)<0) then
                    y(i,j) = 0
                end if
            end do
        end do
#else
        call sub_modelstop("Undefined model calling sub_evaluate")
#endif
#endif
    end FUNCTION sub_evaluate


    SUBROUTINE cashInHand(outputVars)
            !
            ! This subroutine defines tomorrow's "cash in hand" using model parameters.
            !
            ! OUTPUTS: outputVars - a m-by-n matrix representing the cash in hand function
        real(DP), dimension(:,:), INTENT(OUT) :: outputVars
        real(DP), allocatable, dimension(:) :: newInput
        integer :: index_k, index_z

        IF ( (size(outputVars,dim=1)/=capitalCount) .and. (size(outputVars,dim=2)/=shockCount) ) THEN
            PRINT '(a,i3)', 'cashInHand: outputVars must be a matrix of size ',capitalCount,&
                & '-by-',shockCount
            call sub_modelstop('program terminated in cash in hand')
        END IF

        allocate(newInput(2))
        do index_k = 1,capitalCount
            newInput(2)=grid_k(index_k)
            do index_z = 1,shockCount
                newInput(1)=shocks(index_z)
                outputVars(index_k,index_z)=sub_modelOutput(newInput)
            end do
        end do
        deallocate(newInput)

    END SUBROUTINE cashInHand

    FUNCTION getCStar(fnToSolve) RESULT(y)
        ! Get optimal consumption given the value function
        ! Note: The last step is from the envelope condition
        REAL(DP), DIMENSION(capitalCount,shockCount), INTENT(IN) :: fnToSolve
        REAL(DP), DIMENSION(capitalCount,shockCount) :: y

#ifdef NEOCLASSICAL
        REAL(DP), DIMENSION(capitalCount,shockCount) :: D
        INTEGER :: index_k

        !use slope between points as an estimate of the derivative at that point
        D(1,:)=(fnToSolve(2,:)-fnToSolve(1,:))/stepVec(1)
        D(capitalCount,:)=(fnToSolve(capitalCount,:)-fnToSolve(capitalCount-1,:))/stepVec(capitalCount-1)
        do index_k = 2,capitalCount-1
            D(index_k,:)=(fnToSolve(index_k+1,:)-fnToSolve(index_k-1,:))/(stepVec(index_k)+stepVec(index_k-1))
        enddo
        y=D**(-1/tau)
#else
        call sub_modelstop("Undefined model calling getCStar")
        y=fnToSolve
#endif
    END FUNCTION getCStar

    FUNCTION getVStar(Cstar,valuefn) RESULT(y)
        REAL(DP), DIMENSION(capitalCount,shockCount), INTENT(IN) :: valuefn, Cstar
        REAL(DP), DIMENSION(capitalCount,shockCount) :: y

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
    real (DP), parameter :: tau    = 2.d0                           ! risk aversion (for consumption)
    real (DP), parameter :: beta   = 0.9896                         ! Discount factor
    real (DP), parameter :: theta   = 0.357                         ! weight on consumption (if we have labour)
    real (DP), parameter :: rho    = 0.95                           ! Autoregressive
    real (DP), parameter :: sigma  = 0.007                          ! Variance

    !----------------------------------------------------------------
    ! 2. Numerical parameters
    !----------------------------------------------------------------
#ifdef NEOCLASSICAL
    integer :: n_tauchen = 41                                   ! Number of shocks
#endif
#ifdef AIYAGARI
    integer :: n_tauchen = 3                                   ! Number of shocks
#endif

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
        ! OUTPUT: y is the probability states
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
        i=size(stepVec)

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
        !         f_x: values for each x value in "x" input
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
    subroutine sub_endogenize(fnToSolve, grid_k, y, transition, g_k)
            ! The generalized endogenous grid method
            ! INPUTS: fnToSolve - The initial guess of the value function, size m-by-n
            !         grid_k - grid of possible capital values (..,ks,..), size m
            !         y - the possible shocks - size n
            !         transition - the transition matrix - size n-by-n
            ! OUTPUT: fnToSolve - the final value function - size m-by-n
            !         g_k - the solved capital levels

        implicit none

        real(DP), dimension(:,:), intent(inout) :: fnToSolve
        real(DP), dimension(:), intent(in) :: grid_k
        real(DP), dimension(:), intent(in) :: y
        real(DP), dimension(:,:), intent(in) :: transition
        real(DP), dimension(:,:), intent(out) :: g_k

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
        do index_z = 1,n
            kend(:,index_z)=grid_k
        end do
        do while((diff>toler) .and. (iter<MAX_ITER))
            iter=iter+1

            call getNextGuess(fnToSolve,newguess,oncapital,targetcapital)
            do index_z = 1,n
                call sub_myinterp1(oncapital(:,index_z),newguess(:,index_z), targetcapital(:,index_z),&
                    & length_grid_k, Vend1(:,index_z))
            enddo

            value1=sub_evaluate(Vend1)
            diff = maxval(abs(value1-fnToSolve))
            fnToSolve=value1
            if (mod(iter,100)==0) then
                print *, 'subvalue: Iteration: ', iter, 'Tolerance ', diff
                flush(6)
            end if
        enddo

#ifdef NEOCLASSICAL
        ! At this point the program recovers the endogenous grid for capital.
        !Aside: If we have endogenized capital, this function should have no effect
        do index_k = 1,length_grid_k
            do index_z = 1,n_tauchen
                call sub_kendogenousnewton(kend(index_k,index_z),oncapital(index_k,index_z),y(index_z),grid_k(index_k))
            enddo
        enddo

        ! Here we interpolate capital on the grid so that we can compare the results to the standard algorithm.
        !Aside: If we have endogenized capital, the outputs of these to calls should be the same
        do index_z = 1,n_tauchen
            call sub_myinterp1(kend(:,index_z),newguess(:,index_z), grid_k, length_grid_k, fnToSolve(:,index_z))
            call sub_myinterp1(kend(:,index_z),grid_k, grid_k, length_grid_k, g_k(:,index_z))
        !            call sub_myinterp1(kend(:,index_z),Cstar(:,index_z), grid_k, length_grid_k, g_c(:,index_z))
        enddo
#else
        g_k=fnToSolve
#endif

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
#ifdef AIYAGARI
    real(DP):: a_min, a_max
#endif

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
#ifdef NEOCLASSICAL
    call sub_tauchen(y,transition,transitioncdf,rho,sigma,cover_tauchen)
    print *, ' '
    print *, 'For the stochastic shock, Tauchen`s approximation method with', n_tauchen , 'points will be used.'
    print *, ' '
    flush(6)
#endif

#ifdef AIYAGARI
        transition(:,1) = (/.66_dp,.28_dp,.07_dp/)
        transition(:,2)=(/.27_dp,.44_dp,.27_dp/)
        transition(:,3)=(/.07_dp,.28_dp,.66_dp/)
        y= (/.78_dp,1.0_dp,1.27_dp/)
#endif

    cover = 0.25d0         ! Coverage of the grid (+- of steady state k)
    length_grid_k = 1000
    step1 = length_grid_k

    !----------------------------------------------------------------
    ! 1.5 Initialize modules
    !----------------------------------------------------------------

    call initialize(length_grid_k)
    call sub_modelInitialize(length_grid_k,n_tauchen,1)

    !----------------------------------------------------------------
    ! 2. Computation of the Steady State
    !----------------------------------------------------------------
#ifdef NEOCLASSICAL
    allocate(outputVars(3))

    call steadyState(outputVars)

    y_ss = outputVars(1)
    c_ss = outputVars(2)
    k_ss = outputVars(3)

    deallocate(outputVars)
#endif

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

#ifdef NEOCLASSICAL
    call sub_grid_generation(grid_k, k_ss, cover, 2.D0)
#endif
#ifdef AIYAGARI
    a_min=0.0
    a_max=50
    call sub_grid_generation(grid_k, a_min, a_max, 1.D0,1)
#endif

    call initGuess(valueFn)

    !Here we call the subroutine to carry out the endogenous grid
    call sub_endogenize(valuefn,grid_k,y,transition,g_k)

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
