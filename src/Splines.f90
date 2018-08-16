!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                          Futility Development Group                          !
!                             All rights reserved.                             !
!                                                                              !
! Futility is a jointly-maintained, open-source project between the University !
! of Michigan and Oak Ridge National Laboratory.  The copyright and license    !
! can be found in LICENSE.txt in the head directory of this repository.        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief A module interpolating 1D data using splines. Both natural cubic and
!>        monotone cubic splines are available.
!>
!> @par Module Dependencies
!>  - @ref IntrType "IntrType": @copybrief IntrType
!>  - @ref ParameterLists "ParameterLists": @copybrief ParameterLists
!>  - @ref Strings "Strings": @copybrief Strings
!>  - @ref ExceptionHandler "ExceptionHandler": @copybrief ExceptionHandler 
!>
!> @author Gavin Ridley
!>    @date 8/8/2018
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE Splines
  USE IntrType
  USE ParameterLists
  USE Strings
  USE ExceptionHandler

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: newSpline, Spline

  TYPE Spline
    REAL(SRK), ALLOCATABLE, DIMENSION(:) :: knots  ! knots, i.e. ordinates.
    REAL(SRK), ALLOCATABLE, DIMENSION(:) :: values ! f(x) at knot
    REAL(SRK), ALLOCATABLE, DIMENSION(:) :: steps ! delta x of knot
    REAL(SRK), ALLOCATABLE, DIMENSION(:) :: zi ! second derivative at knot. 
    INTEGER(SIK) :: length ! number of knots in spline
    TYPE(StringType) :: splineType
  contains
    procedure :: sample ! gets spline value at x=<arg>
    procedure :: clear ! clear out internal data
    procedure :: solve ! find where spline is equal to some value
  ENDTYPE

  !> module-wide exception handler
  TYPE(ExceptionHandlerType), SAVE :: eSplines

  CHARACTER(LEN=*), PARAMETER :: modname = 'SPLINES'

!===============================================================================
  CONTAINS
!-------------------------------------------------------------------------------

  !> @brief Base subroutine to call for initializing splines.
  !> @param datax ordinates data that become spline knots
  !> @param datay values to be interpolated
  !> @param splineOut Spline object where data is stored.
  !> @param params ParameterList specifying the spline type.
  SUBROUTINE newSpline(datax, datay, splineOut, params)
    USE Strings
    CHARACTER(LEN=*), PARAMETER :: myname = 'newSpline'
    TYPE(ParamType), INTENT(IN) :: params
    REAL(SRK), DIMENSION(:), INTENT(IN) :: datax, datay
    TYPE(Spline), INTENT(INOUT) :: splineOut
    TYPE(StringType) ::  splinetype
    integer(SIK) :: i

    !> Check if spline was already allocated. If so, clear it out.
    if (allocated(splineOut%knots)) then
        call splineOut%clear()
    endif

    !> Allocate spline internal data, make sure length of interpolated
    !> and ordinate value arrays match.
    splineOut%length = size(datax)
    allocate(splineOut%knots(splineOut%length))
    allocate(splineOut%values(splineOut%length))
    allocate(splineOut%steps(splineOut%length-1))

    splineOut%knots = datax
    splineOut%values = datay

    if (size(datay) /= splineOut%length) then
        call eSplines%raiseError(modname//'::'//myname// & 
              'Number of values and knots does not match.')
    endif

    !> ensure knots are strictly increasing
    call twoArraySorter(splineOut%knots, splineOut%values)

    do i = 1, splineOut%length-1
      splineOut%steps(i) = splineOut%knots(i+1) - splineOut%knots(i)
    enddo

    !> get spline type to manufacture
    if (params%has('Spline->type')) then
      call params%get('Spline->type', splinetype)

      if (splinetype == 'naturalcubic') then
        call makenaturalcubic(splineOut)

      else if (splinetype == 'monotonecubic') then
        call makefritsch(splineOut)
      else
        call eSplines%raiseError(modname//'::'//myname//' unknown spline type: '//splinetype)
      endif
    endif
  ENDSUBROUTINE newSpline

  !> Implementation of the cubic splines algorithm in Kincaid and Cheney's
  !> "Numerical Analysis". Amounts to solving a tridiagonal system.
  SUBROUTINE makenaturalcubic(self)
    CLASS(Spline), INTENT(INOUT) :: self
    REAL(SRK), DIMENSION(:), ALLOCATABLE :: bi, ui, vi ! work array
    INTEGER(SIK) :: i

    allocate(bi(self%length), ui(self%length), vi(self%length))
    allocate(self%zi(self%length))

    !> "Steps" make up the sub and superdiagonals.
    do i = 1, self%length-1
      bi(i) = 6.0_SRK / self%steps(i) * (self%values(i+1) - self%values(i))
    enddo

    ! ui terms make up the diagonal
    ui(2) = 2.0_SRK * (self%steps(2)+self%steps(1))
    vi(2) = bi(2) - bi(1)

    ! sweep forward:
    do i = 3, self%length-1
      ui(i) = 2.0_SRK * (self%steps(i)+self%steps(i-1)) - &
              self%steps(i-1)**2/ui(i-1)
      vi(i) = bi(i) - bi(i-1) - self%steps(i-1) * vi(i-1) / ui(i-1)
    enddo

    self%zi = [( 0.0_SRK, i = 1, self%length )]

    ! sweep backward:
    do i = self%length-1, 2, -1
      self%zi(i) = (vi(i) - self%steps(i) * self%zi(i+1)) / ui(i)
    enddo

    deallocate(bi, ui, vi)

    self%splineType = 'naturalcubic'
  ENDSUBROUTINE makenaturalcubic

  !> @brief Follows the algorithm of Moler in "Numerical Computing with MATLAB"
  !>    for forming PCHIP (piecewise cubic hermite interpolating polynomial)
  !>    splines. These have a monotone, aesthetically appealing characteristic.
  !>    Adheres to Fritsch and Carlson's monotonicity conditions described in
  !>    "Monotone Piecewise Cubic Interpolation", 1980
  SUBROUTINE makefritsch(self)
    CLASS(Spline), INTENT(INOUT) :: self
    REAL(SRK), DIMENSION(:), ALLOCATABLE :: slopes
    REAL(SRK) :: alpha, beta, tau, w1, w2
    INTEGER(SIK) :: i

    !> Values of spline's first derivative at each point.
    allocate(self%zi(self%length))

    !> values of secant line derivatives
    allocate(slopes(self%length-1))

    !> set secant slopes
    do i = 1, self%length-1
      slopes(i) = (self%values(i+1)-self%values(i)) / self%steps(i)
    enddo

    !> The method for determining tangent slope values used in
    !> Moler's "Numerical Computing with Matlab"
    do i = 2, self%length-1
      w1 = 2.0_SRK * self%steps(i) + self%steps(i-1)
      w2 = self%steps(i) + 2.0_SRK * self%steps(i-1)
      self%zi(i) = (w1 / slopes(i) + w2 / slopes(i-1)) / (w1+w2)
      self%zi(i) = self%zi(i) ** (-1)
    enddo
    self%zi(1) = slopes(1)
    self%zi(self%length) = slopes(self%length-1)

    !> Loop over subintervals, adjust values of slope of the spline
    !> itself so that the parameters alpha and beta as defined in the
    !> above mentioned paper fall in the region of monotonicity.
    do i = 1, self%length - 1
      if (abs(slopes(i)) < epsilon(alpha)) then
          self%zi(i+1) = 0.0_SRK
      else
        alpha = self%zi(i) / slopes(i)
        beta = self%zi(i+1) / slopes(i)
        if (alpha**2 + beta**2 < 9.0_SRK) then
          continue
        else
          ! Adjust both slopes
          tau = 3.0_SRK / sqrt(alpha**2 + beta**2)
          self%zi(i) = tau * alpha * slopes(i)
          self%zi(i+1) = tau * beta * slopes(i)
        endif
      endif
    enddo

    self%splinetype = 'monotonecubic'
    deallocate(slopes)
  ENDSUBROUTINE makefritsch

  !> @brief Bubble sort two arrays of reals so that the first array is
  !>        strictly increasing. Not the best, but easy to code!
  !> @param arr1 Array to sort so it's increasing only
  !> @param arr2 When value x_i and x_j get swapped in arr1, y_i and y_j in
  !>        arr2 get swapped too
  SUBROUTINE twoArraySorter(arr1, arr2)
    CHARACTER(LEN=*), PARAMETER :: myname = 'twoArraySorter'
    REAL(SRK), DIMENSION(:), INTENT(INOUT) :: arr1, arr2
    integer(SIK) :: leng, i, incr
    leng = size(arr1)
    if (leng /= size(arr2)) then
        call eSplines%raiseError(modname//'::'//myname//' input array length must match.')
    endif

    incr = 1
    i = 1
    do
      !> reached end?
      if (i == leng) exit

      if (arr1(i) < arr1(i+1)) then
        !> ordered as desired?
        i = i + 1
        cycle
      else if (arr1(i) > arr1(i+1)) then
        !> Carry this entry back.
        incr = -1
        do while (i /= 0 .and. arr1(i) > arr1(i+1))
          call swap2(arr1, arr2, i)
          i = i + incr
        enddo
      else
        !> Implies these two knot values are exactly equal, which
        !> is the single unacceptable edge case. The spline matrix
        !> becomes singular here.
        call eSplines%raiseError(modname//'::'//myname//' no two knots may be the same')
      endif
    enddo
  ENDSUBROUTINE twoArraySorter

  !> @brief Swaps value i in arrays arr1 and arr2 with value i+1 in each array.
  !>        Makes the bubble sort's code more concise.
  SUBROUTINE swap2(arr1, arr2, i)
    REAL(SRK), DIMENSION(:), INTENT(INOUT) :: arr1, arr2
    INTEGER(SIK), INTENT(IN) :: i
    REAL(SRK) :: work

    ! array 1
    work = arr1(i+1)
    arr1(i+1) = arr1(i)
    arr1(i) = work

    ! array 2
    work = arr2(i+1)
    arr2(i+1) = arr2(i)
    arr2(i) = work
  ENDSUBROUTINE swap2

  FUNCTION sample(self, xval) RESULT(val)
    REAL(SRK), INTENT(IN) :: xval
    CLASS(Spline), INTENT(IN) :: self 
    REAL(SRK) :: val, c, b, a
    INTEGER(SIK) :: i

    !> Step through subintervals until it's
    !> one step too far, then go back one.
    do i = 1, size(self%steps)
      if (xval - self%knots(i) <= 0.0_SRK) exit
    enddo
    i = i - 1
    
    !> Correct guess on interval the first time?
    if (i == 0) i = 1

    !> Evaluates a cubic polynomial in a segment with both values and
    !> second derivatives given on each side
    if (self%splinetype == 'naturalcubic') then
      a = 1.0_SRK / 6.0_SRK / self%steps(i) * (self%zi(i+1)-self%zi(i))
      b = self%zi(i) / 2.0_SRK
      c = -self%steps(i) / 6.0_SRK * self%zi(i+1) - self%steps(i) / 3.0_SRK &
          * self%zi(i) + 1.0_SRK / self%steps(i) * &
          (self%values(i+1)-self%values(i))
      val = self%values(i) + (xval - self%knots(i)) * ( c +  &
             (xval-self%knots(i))*(b+(xval-self%knots(i))*a))

    !> Evaluates a cubic polynomial with both values and first derivatives
    !> defined on each side.
    else if (self%splinetype == 'monotonecubic') then
      a = (xval - self%knots(i)) / self%steps(i)
      val = (2.0_SRK*a**3 - 3._SRK*a**2 + 1.0_SRK) * self%values(i)
      val = val + (a**3-2._SRK*a**2+a) * self%zi(i) * self%steps(i)
      val = val + (-2._SRK * a**3 + 3._SRK * a**2) * self%values(i+1)
      val = val + (a**3 - a**2) * self%zi(i+1) * self%steps(i)
    endif
  ENDFUNCTION sample

  !> Clear internal allocatable arrays
  SUBROUTINE clear(self)
    CLASS(Spline), INTENT(INOUT) :: self
    deallocate(self%knots, self%values, self%steps, self%zi)
  ENDSUBROUTINE clear

  !> @brief Find where spline is equal to "rhs"
  !>   i.e., find root of spline(x)-rhs = 0
  !>   Uses bisection to start newton method off on a good guess
  !> @param rhs - value to match spline to
  FUNCTION solve(self, rhs) RESULT(root)
    CLASS(Spline), INTENT(IN) :: self
    REAL(SRK), INTENT(IN) :: rhs
    CHARACTER(LEN=*), PARAMETER :: myname = 'solve'

    ! fractional tolerance in root
    REAL(SRK), PARAMETER :: fracdel_toler = 1e-10
    REAL(SRK) :: root, guess1, fracdel, temp

    ! used in bisection
    INTEGER(SIK) :: loc1, loc2, loc3
    LOGICAL(SBK) :: inleft, inright 

    if (self%splinetype /= 'monotonecubic') then
      call eSplines%raiseError(modname//'::'//myname//' Only monotone cubic ' &
        // 'splines are guaranteed to have one root.') 
    endif

    ! Use bisection search to pinpoint root interval:
    loc1 = 1_SIK
    loc2 = self%length
    loc3 = loc2 / 2 ! intentional integer divide
    do
      !> XOR may only be in GNU fortran... probably need some preprocessor
      !> directive to define XOR if other compilers don't have it.
      inleft = XOR(boolsign(self%values(loc1)-rhs), boolsign(self%values(loc3)-rhs))
      inright = XOR(boolsign(self%values(loc2)-rhs), boolsign(self%values(loc3)-rhs))
      if (inleft .AND. inright) then
        call eSplines%raiseError(modname//'::'//myname//'Non-unique root ' &
        // 'detected in solve.')
      else if (inleft) then
        loc2 = loc3
        loc3 = (loc3 + loc1) / 2
      else if (inright) then
        loc1 = loc3
        loc3 = (loc3 + loc2) / 2
      else
        call eSplines%raiseError(modname//'::'//myname//'WTF??')
      endif

      ! Check if bisection search can stop
      if (loc1 == loc3 .or. loc2 == loc3 .or. loc1 == loc2) then
        exit
      endif
    enddo

    ! loc1 now contains the interval the root lies in.
    ! Now just run secant method to get the root.
    guess1 = self%knots(loc1)
    root = self%knots(loc1+1)
    fracdel = abs(guess1-root)/root
    do while (fracdel > fracdel_toler)
      temp = root
      root = (guess1 * (self%sample(root)-rhs) - root * (self%sample(guess1)- &
             rhs)) / (self%sample(root) - self%sample(guess1))
      guess1 = temp
      fracdel = abs(guess1-root)/root
    enddo
  ENDFUNCTION solve

  !> @brief Return true if number of kind SRK is greater than
  !>        zero. Used in bisection method to cut out need for a
  !>        floating point multiply.
  PURE FUNCTION boolsign(num) result(ispositive)
    REAL(SRK), INTENT(IN) :: num
    LOGICAL(SBK) :: ispositive
    if (num > 0.0_SRK) then
      ispositive = .true.
    else
      ispositive = .false.
    endif
  ENDFUNCTION boolsign
     

ENDMODULE Splines
