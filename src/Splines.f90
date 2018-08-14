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
!>  - @ref ParameterLists "ParameterLists": @copybrief ParameterLists
!>
!> @author Gavin Ridley
!>    @date 8/8/2018
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE Splines
  USE IntrType
  USE ParameterLists
  USE Strings

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: newSpline, Spline

  TYPE Spline
    real(SRK), allocatable, dimension(:) :: knots  ! knots, i.e. ordinates.
    real(SRK), allocatable, dimension(:) :: values ! f(x) at knot
    real(SRK), allocatable, dimension(:) :: steps ! delta x of knot
    real(SRK), allocatable, dimension(:) :: zi ! second derivative at knot. 
    logical :: success = .false. ! whether spline has been successfully made.
    type(StringType) :: splineType
  contains
    procedure :: sample ! gets spline value at x=<arg>
    ! procedure :: matches_val ! find root of f(x) - <arg>
    procedure, private :: allocAll ! allocating all arrays in spline
  ENDTYPE
!
!===============================================================================
  CONTAINS
!-------------------------------------------------------------------------------
!> @brief
!> @param datax ordinates data that become spline knots
!> @param datay values to be interpolated
!> @param splineOut Spline object where data is stored.
!> @param params ParameterList specifying the spline type.
  SUBROUTINE newSpline(datax, datay, splineOut, params)
    USE Strings
    IMPLICIT NONE
    TYPE(ParamType), INTENT(IN) :: params
    REAL(SRK), DIMENSION(:), INTENT(IN) :: datax, datay
    REAL(SRK), DIMENSION(:), ALLOCATABLE :: bi, ui, vi ! work array
    TYPE(Spline), INTENT(INOUT) :: splineOut
    TYPE(StringType) ::  splinetype
    integer(SIK) :: i, length

    ! Allocate spline internal data, make sure length of interpolated
    ! and ordinate value arrays match.
    length = size(datax)
    if (size(datay) /= length) then
      print *, 'wtf'
    else
      call splineOut%allocAll(length)
      allocate(bi(length), ui(length), vi(length))
    endif

    splineOut%knots = datax
    splineOut%values = datay

    ! ensure knots are strictly increasing
    call twoArraySorter(splineOut%knots, splineOut%values)

    ! get spline type to manufacture
    if (params%has('Spline->type')) then
      call params%get('Spline->type', splinetype)
      if (splinetype == 'naturalcubic') then
        ! call makenaturalcubic(val, val, ...)
        print *, 'making nat cubic spline'
      endif
    endif

    do i = 1, length-1
      splineOut%steps(i) = splineOut%knots(i+1) - splineOut%knots(i)
      bi(i) = 6.0_SRK * (splineOut%values(i+1) - splineOut%values(i))
    enddo

    ui(2) = 2.0_SRK * (splineOut%steps(2)+splineOut%steps(1))
    vi(2) = bi(2) - bi(1)

    do i = 3, length-1
      ui(i) = 2.0_SRK * (splineOut%steps(i)+splineOut%steps(i-1)) - &
              splineOut%steps(i-1)**2/ui(i-1)
      vi(i) = bi(i) - bi(i-1) - splineOut%steps(i-1) * vi(i-1) / ui(i-1)
    enddo

    splineOut%zi = [( 0.0_SRK, i = 1, length )]

    do i = length-1, 2, -1
      splineOut%zi(i) = (vi(i) - splineOut%steps(i) * splineOut%zi(i+1)) / ui(i)
    enddo

    deallocate(bi, ui, vi)

    splineOut%splineType = 'naturalcubic'
    splineOut%success = .true.

    print *, splineOut%zi
  ENDSUBROUTINE newSpline

  ! Made this a separate private subroutine since I'm unsure about
  ! whether the futility allocation module should be used.
  SUBROUTINE allocAll(self, length)
    IMPLICIT NONE
    CLASS(Spline), INTENT(INOUT) :: self
    INTEGER(SIK) :: length
    allocate(self%knots(length))
    allocate(self%values(length))
    allocate(self%steps(length-1))
    allocate(self%zi(length))
  ENDSUBROUTINE allocAll

  ! Bubble sort an array of reals so that the first array
  ! is strictly increasing, and the second gets swapped
  ! around so that initially corresponding values in each
  ! array still correspond after the sort.
  SUBROUTINE twoArraySorter(arr1, arr2)
    IMPLICIT NONE
    REAL(SRK), DIMENSION(:), INTENT(INOUT) :: arr1, arr2
    integer(SIK) :: leng, i, incr
    leng = size(arr1)
    if (leng /= size(arr2)) then
      print *, 'wtf'
    endif

    incr = 1
    i = 1
    do
      ! reached end?
      if (i == leng) exit

      if (arr1(i) < arr1(i+1)) then
        ! ordered as desired?
        i = i + 1
        cycle
      else if (arr1(i) > arr1(i+1)) then
        ! Carry this entry back.
        incr = -1
        do while (i /= 0 .and. arr1(i) < arr1(i+1))
          call swap2(arr1, arr2, i)
          i = i + incr
        enddo
      else
        ! Implies these two knot values are exactly equal, which
        ! is the single unacceptable edge case. The spline matrix
        ! becomes singular here.
        print *, 'wtf'
      endif
    enddo
  ENDSUBROUTINE twoArraySorter

  ! Swaps value i in arrays arr1 and arr2 with value i+1 in each array.
  SUBROUTINE swap2(arr1, arr2, i)
    IMPLICIT NONE
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
    IMPLICIT NONE
    REAL(SRK), INTENT(IN) :: xval
    CLASS(Spline), INTENT(IN) :: self 
    REAL(SRK) :: val, c, b, a
    INTEGER(SIK) :: i

    ! Step through subintervals until it's
    ! one step too far, then go back one.
    do i = 1, size(self%steps)
      if (xval - self%knots(i) <= 0.0_SRK) exit
    enddo
    i = i - 1
    

    ! Correct guess on interval the first time?
    if (i == 0) i = 1

    a = 1.0_SRK / 6.0_SRK / self%steps(i) * (self%zi(i+1)-self%zi(i))
    b = self%zi(i) / 2.0_SRK
    c = -self%steps(i) / 6.0_SRK * self%zi(i+1) - self%steps(i) / 3.0_SRK &
        * self%zi(i) + 1.0_SRK / self%steps(i) * &
        (self%values(i+1)-self%values(i))
    val = self%values(i) + (xval - self%knots(i)) * ( c +  &
           (xval-self%knots(i))*(b+(xval-self%knots(i))*a))
  ENDFUNCTION sample

ENDMODULE Splines
