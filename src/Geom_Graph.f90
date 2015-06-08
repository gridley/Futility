!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                              Copyright (C) 2012                              !
!                   The Regents of the University of Michigan                  !
!              MPACT Development Group and Prof. Thomas J. Downar              !
!                             All rights reserved.                             !
!                                                                              !
! Copyright is reserved to the University of Michigan for purposes of          !
! controlled dissemination, commercialization through formal licensing, or     !
! other disposition. The University of Michigan nor any of their employees,    !
! makes any warranty, express or implied, or assumes any liability or          !
! responsibility for the accuracy, completeness, or usefulness of any          !
! information, apparatus, product, or process disclosed, or represents that    !
! its use would not infringe privately owned rights. Reference herein to any   !
! specific commercial products, process, or service by trade name, trademark,  !
! manufacturer, or otherwise, does not necessarily constitute or imply its     !
! endorsement, recommendation, or favoring by the University of Michigan.      !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief A Fortran 2003 module defining a graph type.
!>
!>
!> @par Module Dependencies
!>  - @ref IntrType "IntrType": @copybrief IntrType
!>  - @ref Allocs "Allocs": @copybrief Allocs
!>
!> @author Brendan Kochunas
!>    @date 06/06/2015
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE Geom_Graph
  USE IntrType
  USE Allocs
  
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: GraphType
  
  !> @brief a Planar Graph
  TYPE :: GraphType
    !> A list of vertices in the graph.
    !> The list is sorted lexicographically
    !> e.g. V_i with (x_i,y_i) and V_i+1 with (x_i+1,y_i+1) then
    !> x_i < x_i+1 or if x_i == x_i+1 then y_i < y_i+1
    !> The insert vertex routine inserts the vertex in order. Duplicate points
    !> are not stored.
    REAL(SRK),ALLOCATABLE :: vertices(:,:)
    !> Matrix indicating connectivity of graph
    !> 0 means no connection, 1 means linear connection,
    !> -1 means quadratic connection. Diagonal is 0 and matrix is symmetric.
    !> size is (nvertices,nvertices)
    INTEGER(SIK),ALLOCATABLE :: edgeMatrix(:,:)
    !> Similar to edgeMatrix component except for entries with -1 it stores
    !> the center of rotation and radius to define the quadratic edge.
    !> size is (3,nvertices,nvertices)
    REAL(SRK),ALLOCATABLE :: quadEdges(:,:,:)
    CONTAINS
      !> @copybrief Geom_Graph::nVert_graphType
      !> @copydetails Geom_Graph::nVert_graphType
      PROCEDURE,PASS :: nVert => nVert_graphType
      !> @copybrief Geom_Graph::nEdge_graphType
      !> @copydetails Geom_Graph::nEdge_graphType
      PROCEDURE,PASS :: nEdge => nEdge_graphType
      !> @copybrief Geom_Graph::getVertIndex_graphType
      !> @copydetails Geom_Graph::getVertIndex_graphType
      PROCEDURE,PASS :: getVertIndex => getVertIndex_graphType
      !> @copybrief Geom_Graph::nAdjacent_graphType
      !> @copydetails Geom_Graph::nAdjacent_graphType
      PROCEDURE,PASS :: nAdjacent => nAdjacent_graphType
      !> @copybrief Geom_Graph::insertVertex_graphType
      !> @copydetails Geom_Graph::insertVertex_graphType
      PROCEDURE,PASS :: insertVertex => insertVertex_graphType
      !> @copybrief Geom_Graph::defineEdge_graphType
      !> @copydetails Geom_Graph::defineEdge_graphType
      PROCEDURE,PASS :: defineEdge => defineEdge_graphType
      !> @copybrief Geom_Graph::defineQuadEdge_graphType
      !> @copydetails Geom_Graph::defineQuadEdge_graphType
      PROCEDURE,PASS :: defineQuadraticEdge => defineQuadEdge_graphType
      !> @copybrief Geom_Graph::removeVertex_graphType
      !> @copydetails Geom_Graph::removeVertex_graphType
      PROCEDURE,PASS :: removeVertex => removeVertex_graphType
      !> @copybrief Geom_Graph::removeVertex_idx_graphType
      !> @copydetails Geom_Graph::removeVertex_idx_graphType
      PROCEDURE,PASS :: removeVertexI => removeVertex_idx_graphType
      !> @copybrief Geom_Graph::removeEdge_graphType
      !> @copydetails Geom_Graph::removeEdge_graphType
      PROCEDURE,PASS :: removeEdge => removeEdge_graphType
      !> @copybrief Geom_Graph::removeVertex_idx_graphType
      !> @copydetails Geom_Graph::removeVertex_idx_graphType
      PROCEDURE,PASS :: removeEdgeIJ => removeEdge_IJ_graphType
      !> @copybrief Geom_Graph::getMCB_graphType
      !> @copydetails Geom_Graph::getMCB_graphType
      PROCEDURE,PASS :: getMCB => getMCB_graphType
      !> @copybrief Geom_Graph::clear_graphType
      !> @copydetails Geom_Graph::clear_graphType
      PROCEDURE,PASS :: clear => clear_graphType
  ENDTYPE GraphType
!
!===============================================================================
  CONTAINS
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    ELEMENTAL FUNCTION nVert_graphType(thisGraph) RESULT(n)
      CLASS(GraphType),INTENT(IN) :: thisGraph
      INTEGER(SIK) :: n
      n=0
      IF(ALLOCATED(thisGraph%vertices)) n=SIZE(thisGraph%vertices,DIM=2)
    ENDFUNCTION nVert_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    ELEMENTAL FUNCTION nEdge_graphType(thisGraph) RESULT(n)
      CLASS(GraphType),INTENT(IN) :: thisGraph
      INTEGER(SIK) :: n
      n=0
      IF(ALLOCATED(thisGraph%edgeMatrix)) n=SUM(ABS(thisGraph%edgeMatrix))/2
    ENDFUNCTION nEdge_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    PURE FUNCTION getVertIndex_graphType(thisGraph,coord) RESULT(idx)
      CLASS(GraphType),INTENT(IN) :: thisGraph
      REAL(SRK),INTENT(IN) :: coord(2)
      INTEGER(SIK) :: idx
      INTEGER(SIK) :: i,j,n
      idx=-1
      n=nVert_graphType(thisGraph)
      DO i=1,n
        IF(coord(1) .APPROXEQA. thisGraph%vertices(1,i)) THEN
          IF(coord(2) .APPROXEQA. thisGraph%vertices(2,i)) THEN
            idx=i
          ELSE
            DO j=i+1,n
              IF(coord(2) .APPROXEQA. thisGraph%vertices(2,j)) THEN
                idx=j
                EXIT
              ENDIF
            ENDDO
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDFUNCTION getVertIndex_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    ELEMENTAL FUNCTION nAdjacent_graphType(thisGraph,i) RESULT(n)
      CLASS(GraphType),INTENT(IN) :: thisGraph
      INTEGER(SIK),INTENT(IN) :: i
      INTEGER(SIK) :: n
      n=0
      IF(ALLOCATED(thisGraph%edgeMatrix)) THEN
        IF(0 < i .AND. i < SIZE(thisGraph%edgeMatrix,DIM=2)+1) &
          n=SUM(ABS(thisGraph%edgeMatrix(:,i)))
      ENDIF
    ENDFUNCTION nAdjacent_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    ELEMENTAL FUNCTION isFilament_graphType(thisGraph) RESULT(bool)
      CLASS(GraphType),INTENT(IN) :: thisGraph
      LOGICAL(SBK) :: bool
      bool=.FALSE.
    ENDFUNCTION isFilament_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    SUBROUTINE insertVertex_graphType(thisGraph,coord)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      REAL(SRK),INTENT(IN) :: coord(2)
      INTEGER(SIK) :: i,j,n,k
      INTEGER(SIK),ALLOCATABLE :: tmpE(:,:)
      REAL(SRK),ALLOCATABLE :: tmpVertices(:,:),tmpQE(:,:,:)
      IF(ALLOCATED(thisGraph%vertices)) THEN
        n=SIZE(thisGraph%vertices,DIM=2)
        CALL dmallocA(tmpVertices,2,n+1)
        j=0
        DO i=1,n
          IF(coord(1) <= thisGraph%vertices(1,i)) THEN
            IF(coord(1) .APPROXEQA. thisGraph%vertices(1,i)) THEN
              IF(coord(2) .APPROXEQA. thisGraph%vertices(2,i)) THEN
                k=-1 !Duplicate vertex
              ELSEIF(coord(2) < thisGraph%vertices(2,i)) THEN                
                k=i !Before i
              ELSE
                !After i
                DO j=i+1,n
                  !Find index for end of sequence with same x value
                  IF(.NOT.(coord(1) .APPROXEQA. thisGraph%vertices(1,j))) EXIT
                ENDDO
                !Search on y through sequence of same x
                DO k=i+1,j-1
                  IF(coord(2) < thisGraph%vertices(2,k)) EXIT
                ENDDO
              ENDIF
            ELSE
              k=i !Before i
            ENDIF
            EXIT
          ENDIF
        ENDDO
        IF(j /= 0) i=j
        IF(i == n+1) THEN
          k=n+1 !Last point
          i=n
        ENDIF
        IF(k > 0) THEN
          IF(k > 1) tmpVertices(:,1:k-1)=thisGraph%vertices(:,1:i-1)
          tmpVertices(:,k)=coord
          tmpVertices(:,k+1:n+1)=thisGraph%vertices(:,i:n)
          CALL demallocA(thisGraph%vertices)
          CALL MOVE_ALLOC(tmpVertices,thisGraph%vertices)
          
          !Expand Edge Matrices
          CALL dmallocA(tmpE,n+1,n+1)
          CALL dmallocA(tmpQE,3,n+1,n+1)
          DO j=1,k-1
            DO i=1,k-1
              tmpE(i,j)=thisGraph%edgeMatrix(i,j)
              tmpE(j,i)=thisGraph%edgeMatrix(j,i)
              tmpQE(:,i,j)=thisGraph%quadEdges(:,i,j)
              tmpQE(:,j,i)=thisGraph%quadEdges(:,j,i)
            ENDDO
            DO i=k+1,n+1
              tmpE(i,j)=thisGraph%edgeMatrix(i-1,j)
              tmpE(j,i)=thisGraph%edgeMatrix(j,i-1)
              tmpQE(:,i,j)=thisGraph%quadEdges(:,i-1,j)
              tmpQE(:,j,i)=thisGraph%quadEdges(:,j,i-1)
            ENDDO
          ENDDO
          DO j=k+1,n+1
            DO i=1,k-1
              tmpE(i,j)=thisGraph%edgeMatrix(i,j-1)
              tmpE(j,i)=thisGraph%edgeMatrix(j-1,i)
              tmpQE(:,i,j)=thisGraph%quadEdges(:,i,j-1)
              tmpQE(:,j,i)=thisGraph%quadEdges(:,j-1,i)
            ENDDO
            DO i=k+1,n+1
              tmpE(i,j)=thisGraph%edgeMatrix(i-1,j-1)
              tmpE(j,i)=thisGraph%edgeMatrix(j-1,i-1)
              tmpQE(:,i,j)=thisGraph%quadEdges(:,i-1,j-1)
              tmpQE(:,j,i)=thisGraph%quadEdges(:,j-1,i-1)
            ENDDO
          ENDDO
          CALL demallocA(thisGraph%edgeMatrix)
          CALL MOVE_ALLOC(tmpE,thisGraph%edgeMatrix)
          CALL demallocA(thisGraph%quadEdges)
          CALL MOVE_ALLOC(tmpQE,thisGraph%quadEdges)
        ENDIF
      ELSE
        CALL dmallocA(thisGraph%vertices,2,1)
        thisGraph%vertices(:,1)=coord
        CALL dmallocA(thisGraph%edgeMatrix,1,1)
        CALL dmallocA(thisGraph%quadEdges,3,1,1)
      ENDIF
    ENDSUBROUTINE insertVertex_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    PURE SUBROUTINE defineEdge_graphType(thisGraph,coord1,coord2)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      REAL(SRK),INTENT(IN) :: coord1(2)
      REAL(SRK),INTENT(IN) :: coord2(2)
      INTEGER(SIK) :: v1,v2
      v1=getVertIndex_graphType(thisGraph,coord1)
      v2=getVertIndex_graphType(thisGraph,coord2)
      IF(v1 > 0 .AND. v2 > 0 .AND. v1 /= v2) THEN
        thisGraph%edgeMatrix(v1,v2)=1
        thisGraph%edgeMatrix(v2,v1)=1
      ENDIF
    ENDSUBROUTINE defineEdge_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    PURE SUBROUTINE defineQuadEdge_graphType(thisGraph,coord1,coord2,c0,r)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      REAL(SRK),INTENT(IN) :: coord1(2)
      REAL(SRK),INTENT(IN) :: coord2(2)
      REAL(SRK),INTENT(IN) :: c0(2)
      REAL(SRK),INTENT(IN) :: r
      
      INTEGER(SIK) :: v1,v2
      REAL(SRK) :: x1,y1,x2,y2,r1,r2,rsq
      
      !Check that coord1 and coord2 exist on circle
      x1=coord1(1)-c0(1)
      y1=coord1(2)-c0(2)
      r1=x1*x1+y1*y1
      x2=coord2(1)-c0(1)
      y2=coord2(2)-c0(2)
      r2=x2*x2+y2*y2
      rsq=r*r
      IF((rsq .APPROXEQA. r1) .AND. (rsq .APPROXEQA. r2)) THEN
        v1=getVertIndex_graphType(thisGraph,coord1)
        v2=getVertIndex_graphType(thisGraph,coord2)
        IF(v1 > 0 .AND. v2 > 0 .AND. v1 /= v2) THEN
          !Update edge matrix
          thisGraph%edgeMatrix(v1,v2)=-1
          thisGraph%edgeMatrix(v2,v1)=-1

          !Store circle info in quadEdges
          thisGraph%quadEdges(1:2,v1,v2)=c0
          thisGraph%quadEdges(3,v1,v2)=r
          thisGraph%quadEdges(1:2,v2,v1)=c0
          thisGraph%quadEdges(3,v2,v1)=r
        ENDIF
      ENDIF
    ENDSUBROUTINE defineQuadEdge_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    SUBROUTINE removeVertex_graphType(thisGraph,v)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      REAL(SRK),INTENT(IN) :: v(2)
      INTEGER(SIK) :: i
      i=getVertIndex_graphType(thisGraph,v)
      CALL removeVertex_idx_graphType(thisGraph,i)
    ENDSUBROUTINE removeVertex_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    SUBROUTINE removeVertex_idx_graphType(thisGraph,idx)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      INTEGER(SIK),INTENT(IN) :: idx
      
      INTEGER(SIK) :: i,j,n
      INTEGER(SIK),ALLOCATABLE :: tmpEdge(:,:)
      REAL(SRK),ALLOCATABLE :: tmpVert(:,:),tmpQE(:,:,:)
      
      
      n=nVert_graphType(thisGraph)
      IF(0 < idx .AND. idx <= n) THEN
        CALL dmallocA(tmpVert,2,n-1)
        CALL dmallocA(tmpEdge,n-1,n-1)
        CALL dmallocA(tmpQE,3,n-1,n-1)
        
        DO i=1,idx-1
          tmpVert(:,i)=thisGraph%vertices(:,i)
          DO j=1,idx-1
            tmpEdge(j,i)=thisGraph%edgeMatrix(j,i)
            tmpQE(:,j,i)=thisGraph%quadEdges(:,j,i)
          ENDDO
          DO j=idx+1,n
            tmpEdge(j-1,i)=thisGraph%edgeMatrix(j,i)
            tmpQE(:,j-1,i)=thisGraph%quadEdges(:,j,i)
          ENDDO
        ENDDO
        
        DO i=idx+1,n
          tmpVert(:,i-1)=thisGraph%vertices(:,i)
          DO j=1,idx-1
            tmpEdge(j,i-1)=thisGraph%edgeMatrix(j,i)
            tmpQE(:,j,i-1)=thisGraph%quadEdges(:,j,i)
          ENDDO
          DO j=idx+1,n
            tmpEdge(j-1,i-1)=thisGraph%edgeMatrix(j,i)
            tmpQE(:,j-1,i-1)=thisGraph%quadEdges(:,j,i)
          ENDDO
        ENDDO
        
        CALL thisGraph%clear()
        CALL MOVE_ALLOC(tmpVert,thisGraph%vertices)
        CALL MOVE_ALLOC(tmpEdge,thisGraph%edgeMatrix)
        CALL MOVE_ALLOC(tmpQE,thisGraph%quadEdges)
      ENDIF
    ENDSUBROUTINE removeVertex_idx_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    PURE SUBROUTINE removeEdge_graphType(thisGraph,c1,c2)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      REAL(SRK),INTENT(IN) :: c1(2)
      REAL(SRK),INTENT(IN) :: c2(2)
      
      INTEGER(SIK) :: v1,v2
      
      v1=getVertIndex_graphType(thisGraph,c1)
      v2=getVertIndex_graphType(thisGraph,c2)
      IF(v1 > 0 .AND. v2 > 0) THEN
        thisGraph%edgeMatrix(v1,v2)=0
        thisGraph%edgeMatrix(v2,v1)=0
        thisGraph%quadEdges(:,v1,v2)=0.0_SRK
        thisGraph%quadEdges(:,v2,v1)=0.0_SRK
      ENDIF
    ENDSUBROUTINE removeEdge_graphType
   
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    ELEMENTAL SUBROUTINE removeEdge_IJ_graphType(thisGraph,i,j)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      INTEGER(SIK),INTENT(IN) :: i
      INTEGER(SIK),INTENT(IN) :: j
      
      INTEGER(SIK) :: n
      
      n=nVert_graphType(thisGraph)+1
      IF(i > 0 .AND. j > 0 .AND. i < n .AND. j < n) THEN
        thisGraph%edgeMatrix(i,j)=0
        thisGraph%edgeMatrix(j,i)=0
        thisGraph%quadEdges(:,i,j)=0.0_SRK
        thisGraph%quadEdges(:,j,i)=0.0_SRK
      ENDIF
    ENDSUBROUTINE removeEdge_IJ_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    SUBROUTINE removeFilament_vertIdx_graphType(thisGraph,i)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      INTEGER(SIK),INTENT(IN) :: i
      
      LOGICAL(SBK) :: loop2
      INTEGER(SIK) :: j,n,nAdj,nVerts
      INTEGER(SIK),ALLOCATABLE :: filVerts(:)
      
      n=nVert_graphType(thisGraph)
      ALLOCATE(filVerts(n))
      IF(0 < i .AND. i <= n) THEN
        nVerts=0
        nAdj=nAdjacent_graphType(thisGraph,i)
        DO WHILE(nAdj == 1)
          loop2=.TRUE.
          DO j=1,i-1
            IF(thisGraph%edgeMatrix(j,i) /= 0) THEN
              loop2=.FALSE.
              nVerts=nVerts+1
              filVerts(nVerts)=j
              EXIT
            ENDIF
          ENDDO
          IF(loop2) THEN
            DO j=i+1,n
              IF(thisGraph%edgeMatrix(j,i) /= 0) THEN
                nVerts=nVerts+1
                filVerts(nVerts)=j
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        
        DO j=1,nVerts
          CALL removeVertex_idx_graphType(thisGraph,filVerts(j))
        ENDDO
      ENDIF
    ENDSUBROUTINE removeFilament_vertIdx_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    SUBROUTINE removeFilament_graph_graphType(thisGraph,subgraph)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      CLASS(GraphType),INTENT(IN) :: subgraph
      
    ENDSUBROUTINE removeFilament_graph_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    SUBROUTINE extractPrimitive_graphType(thisGraph,v0,subgraph)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      INTEGER(SIK) :: v0
      CLASS(GraphType),INTENT(IN) :: subgraph
      
    ENDSUBROUTINE extractPrimitive_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    SUBROUTINE getMCB_graphType(thisGraph,cycles)
      CLASS(GraphType),INTENT(IN) :: thisGraph
      TYPE(GraphType),ALLOCATABLE :: cycles(:)
      
      INTEGER(SIK) :: i,n,nadj,ncycles
      TYPE(GraphType) :: g,primeGraph
      TYPE(GraphType),ALLOCATABLE :: tmpCycles(:)
      
      IF(ALLOCATED(cycles)) THEN
        DO i=1,n
          CALL cycles(i)%clear()
        ENDDO
        DEALLOCATE(cycles)
      ENDIF
      SELECTTYPE(thisGraph); TYPE IS(GraphType)
        g=thisGraph
      ENDSELECT
      ncycles=0
      DO WHILE(g%nVert() > 0)
        nadj=nAdjacent_graphType(g,1)
        IF(nadj == 0) THEN
          CALL removeVertex_idx_graphType(g,1)
        ELSEIF(nadj == 1) THEN
          CALL removeFilament_vertIdx_graphType(g,1)
        ELSE
          CALL extractPrimitive_graphType(g,1,primeGraph)
          IF(isFilament_GraphType(primeGraph)) THEN
            CALL removeFilament_graph_graphType(g,primeGraph)
          ELSE
            !Found minimum cycle, so add it to basis
            ncycles=ncycles+1
            ALLOCATE(tmpCycles(ncycles))
            DO i=1,ncycles-1
              tmpCycles(i)=cycles(i)
              CALL cycles(i)%clear()
            ENDDO
            tmpCycles(ncycles)=primeGraph
            DEALLOCATE(cycles)
            CALL MOVE_ALLOC(tmpCycles,cycles)
          ENDIF
          CALL removeEdge_GraphType(g,primeGraph%vertices(:,1), &
            primeGraph%vertices(:,2))
          CALL primeGraph%clear()
        ENDIF
      ENDDO
    ENDSUBROUTINE getMCB_graphType
!
!-------------------------------------------------------------------------------
!> @brief
!> @param
!>
!>
!>
    SUBROUTINE clear_graphType(thisGraph)
      CLASS(GraphType),INTENT(INOUT) :: thisGraph
      CALL demallocA(thisGraph%vertices)
      CALL demallocA(thisGraph%edgeMatrix)
      CALL demallocA(thisGraph%quadEdges)
    ENDSUBROUTINE clear_graphType
!
ENDMODULE Geom_Graph
