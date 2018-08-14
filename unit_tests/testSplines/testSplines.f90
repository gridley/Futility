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
  REAL(SRK) :: work
  INTEGER(SIK) :: i

  CALL params%add('Spline->type', 'naturalcubic')

  knots = [ (real(i, SRK)*3._SRK, i = 1, 4) ]
  values = [ (real(i**2, SRK), i = 1, 4) ]
  values(3) = 0.0_SRK

  ! dbg
  open(9, file='out')

  CREATE_TEST('TEST SPLINES')

  ! make the spline
  call newSpline(knots, values, testspline, params)

  ! Spot check a few expected natural spline values:
  COMPONENT_TEST('Natural cubic spline construction')

  do i = 0, 90
    work = real(i, KIND=SRK) * .10_SRK + 3.0_SRK
    write(9,*) work, testspline%sample(work)
  enddo

  ASSERT(1.0_SRK .APPROXEQ. 1.0_SRK, 'testspline%sample(5.5)')

  COMPONENT_TEST('Monotone cubic spline construction')
  ASSERT(.true., 'testspline%sample(5.5) (monotone)')

  ! dbg
  close(9)

  FINALIZE_TEST()
ENDPROGRAM testSplines 
