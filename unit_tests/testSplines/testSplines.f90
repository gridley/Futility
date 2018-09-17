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

  REAL(SRK), PARAMETER, DIMENSION(90) :: pchipders = [0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.0000000000000000_SRK, \
  0.54199999999999815_SRK, \
  1.0480000000000007_SRK, \
  1.5180000000000031_SRK, \
  1.9520000000000015_SRK, \
  2.3500000000000001_SRK, \
  2.7119999999999984_SRK, \
  3.0380000000000003_SRK, \
  3.3280000000000025_SRK, \
  3.5820000000000007_SRK, \
  3.8000000000000007_SRK, \
  3.9820000000000011_SRK, \
  4.1279999999999992_SRK, \
  4.2379999999999995_SRK, \
  4.3120000000000003_SRK, \
  4.3500000000000005_SRK, \
  4.3519999999999994_SRK, \
  4.3179999999999987_SRK, \
  4.2480000000000002_SRK, \
  4.1420000000000012_SRK, \
  4.0000000000000009_SRK, \
  3.8219999999999956_SRK, \
  3.6080000000000041_SRK, \
  3.3579999999999974_SRK, \
  3.0719999999999974_SRK, \
  2.7499999999999996_SRK, \
  2.3919999999999955_SRK, \
  1.9980000000000022_SRK, \
  1.5679999999999970_SRK, \
  1.1020000000000001_SRK, \
  0.59999999999999998_SRK, \
  0.56533333333333280_SRK, \
  0.53244444444444461_SRK, \
  0.50133333333333330_SRK, \
  0.47199999999999998_SRK, \
  0.44444444444444464_SRK, \
  0.41866666666666646_SRK, \
  0.39466666666666689_SRK, \
  0.37244444444444497_SRK, \
  0.35200000000000031_SRK, \
  0.33333333333333298_SRK, \
  0.31644444444444397_SRK, \
  0.30133333333333340_SRK, \
  0.28800000000000009_SRK, \
  0.27644444444444399_SRK, \
  0.26666666666666689_SRK, \
  0.25866666666666599_SRK, \
  0.25244444444444492_SRK, \
  0.24800000000000036_SRK, \
  0.24533333333333332_SRK, \
  0.24444444444444416_SRK, \
  0.24533333333333332_SRK, \
  0.24799999999999969_SRK, \
  0.25244444444444453_SRK, \
  0.25866666666666666_SRK, \
  0.26666666666666711_SRK, \
  0.27644444444444427_SRK, \
  0.28800000000000009_SRK, \
  0.30133333333333345_SRK, \
  0.31644444444444436_SRK, \
  0.33333333333333331_SRK]

  CALL params%add('Spline->type', 'naturalcubic')

  knots = [ (real(i, SRK)*3._SRK, i = 1, 4) ]
  values = [ (real(i**2, SRK), i = 1, 4) ]
  values(3) = 0.0_SRK

  CREATE_TEST('SPLINES')

  ! make the spline
  ! Just has to not throw any errors here.
  ! -----------------------------------------------------------------
  COMPONENT_TEST('Natural cubic spline construction')
  testspline = Spline(knots, values, params)

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
  testspline = Spline(knots, values, params)

  ! spot-check value
  ASSERT(testspline%sample(10.0_SRK) .APPROXEQ. 10.451851851851851_SRK, 'PCHIP spot check')

  !> Sample many values from the spline, check for monotonicity
  work = 0.0_SRK
  do i = 0, 90
    lastvalue = work
    work = real(i, KIND=SRK) * .10_SRK + 3.0_SRK
    ASSERT(testspline%sample(work) >= testspline%sample(lastvalue), 'PCHIP monotonicity')
  enddo

  !> Check derivatives of the monotone cubic splines
  COMPONENT_TEST('PCHIP derivatives')
  do i = 1, 90
    work = real(i, KIND=SRK) * 0.1_SRK + 3.0_SRK
    ASSERT( testspline%puresample(work, 1) .APPROXEQ. pchipders(i), 'PCHIP derivative value check')
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
  testspline = Spline(knots, values, params)
  ASSERT(testspline%sample(5.5_SRK) .APPROXEQ. 4.3148148148148149_SRK, 'nat. cubic check at 5.5')
  ASSERT(testspline%sample(10.0_SRK) .APPROXEQ. 3.1851851851851851_SRK, 'nat. cubic check at 10.0)')
  ASSERT(testspline%sample(3.0_SRK) .APPROXEQ. 1.0_SRK, 'testspline%sample(1.)')


  FINALIZE_TEST()
ENDPROGRAM testSplines 
