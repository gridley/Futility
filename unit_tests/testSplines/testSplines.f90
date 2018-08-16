!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                          Futility Development Group                          !
!                             All rights reserved.                             !
!                                                                              !
! Futility is a jointly-maintained, open-source project between the University !
! of Michigan and Oak Ridge National Laboratory.  The copyright and license    !
! can be found in LICENSE.txt in the head directory of this repository.        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PROGRAM testSplines
#include "UnitTest.h"
  USE UnitTest
  USE ParameterLists
  USE IntrType
  USE Splines

  IMPLICIT NONE

  TYPE(ParamType) :: params
  TYPE(Spline) :: testspline
  REAL(SRK), DIMENSION(4) :: knots, values
  REAL(SRK) :: work, lastvalue
  INTEGER(SIK) :: i

  CALL params%add('Spline->type', 'naturalcubic')

  knots = [ (real(i, SRK)*3._SRK, i = 1, 4) ]
  values = [ (real(i**2, SRK), i = 1, 4) ]
  values(3) = 0.0_SRK

  CREATE_TEST('SPLINES')

  ! make the spline
  ! Just has to not throw any errors here.
  ! -----------------------------------------------------------------
  COMPONENT_TEST('Natural cubic spline construction')
  call newSpline(knots, values, testspline, params)

  ! Spot check a few expected natural spline values:
  ! -----------------------------------------------------------------
  COMPONENT_TEST('Natural cubic spline values')
  ASSERT(testspline%sample(5.5_SRK) .APPROXEQ. 4.3148148148148149_SRK, 'nat. cubic check at 5.5')
  ASSERT(testspline%sample(10.0_SRK) .APPROXEQ. 3.1851851851851851_SRK, 'nat. cubic check at 10.0)')
  ASSERT(testspline%sample(3.0_SRK) .APPROXEQ. 1.0_SRK, 'testspline%sample(1.)')

  ! PCHIP spline, should be monotonic even with a large jump in the data
  ! -----------------------------------------------------------------
  COMPONENT_TEST('Monotone cubic spline construction')
  knots = [ (real(i, SRK)*3._SRK, i = 1, 4) ]
  values = [ 1.0_SRK, 1.0_SRK, 10.0_SRK, 11.00_SRK ]

  CALL params%clear()
  CALL params%add('Spline->type', 'monotonecubic')
  CALL newSpline(knots, values, testspline, params)

  ! spot-check value
  ASSERT(testspline%sample(10.0_SRK) .APPROXEQ. 10.451851851851851_SRK, 'PCHIP spot check')

  !> Sample many values from the spline, check for monotonicity
  work = 0.0_SRK
  do i = 0, 90
    lastvalue = work
    work = real(i, KIND=SRK) * .10_SRK + 3.0_SRK
    ASSERT(work >= lastvalue, 'PCHIP monotonicity')
  enddo

  ! -----------------------------------------------------------------
  COMPONENT_TEST('Root finding on splines')
  work = testspline%solve(5.0_SRK)
  ASSERT(work .APPROXEQ. 7.4366543036262742_SRK, 'where monotone spline matches 5.0')

  ! -----------------------------------------------------------------
  COMPONENT_TEST('Forming cubic spline with nonincreasing input data')
  knots = [ 3.0_SRK, 9.0_SRK, 6.0_SRK, 12.0_SRK ]
  values = [ 1.0_SRK, 0.0_SRK, 4.0_SRK, 16.0_SRK ]
  CALL params%clear()
  CALL params%add('Spline->type', 'naturalcubic')
  call newSpline(knots, values, testspline, params)
  ASSERT(testspline%sample(5.5_SRK) .APPROXEQ. 4.3148148148148149_SRK, 'nat. cubic check at 5.5')
  ASSERT(testspline%sample(10.0_SRK) .APPROXEQ. 3.1851851851851851_SRK, 'nat. cubic check at 10.0)')
  ASSERT(testspline%sample(3.0_SRK) .APPROXEQ. 1.0_SRK, 'testspline%sample(1.)')

  FINALIZE_TEST()
ENDPROGRAM testSplines 
