!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.
!
! Adapted to HiSPEET by Joerg Stiller (2024)

module KD_Tree
  use KD_Tree__Precision, only: KDTree_RP
  use KD_Tree__Priority_Queue
  implicit none
  private

  !------------------REAL KIND PARAMETER--------------------------

  public :: KDTree_RP

  !---------------------------------------------------------------
  !
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric is supported.
  !
  !
  ! This module is identical to 'kdtree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  dat(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is dat(1:N,1:D)
  ! as in the original kdtree module.
  !
  !-------------DATA TYPE, CREATION, DELETION---------------------

  public :: KDTree
  public :: KDTree_Result
  public :: KDTree_Node
  public :: KDTree_Create
  public :: KDTree_Destroy

  !-------------------SEARCH ROUTINES-----------------------------

  public :: KDTree_N_Nearest
  public :: KDTree_N_NearestAroundPoint
  ! Return fixed number of nearest neighbors around arbitrary vector,
  ! or extant point in dataset, with decorrelation window.
  !
  public :: KDTree_R_Nearest
  public :: KDTree_R_NearestAroundPoint
  ! Return points within a fixed ball of arb vector/extant point
  !
  public :: KDTree_SortResults
  ! Sort, in order of increasing distance, results from above.
  !
  public :: KDTree_R_Count
  public :: KDTree_R_CountAroundPoint
  ! Count points within a fixed ball of arb vector/extant point
  !
  public :: KDTree_N_NearestBruteForce
  public :: KDTree_R_NearestBruteForce
  ! brute force of KDTree_[n|r]_nearest

  !----------------------------------------------------------------

  integer, parameter :: bucket_size = 12
  ! The maximum number of points to keep in a terminal node.

  type Interval
      real(KDTree_RP) :: lower,upper
  end type Interval

  type KDTree_Node
      ! an internal tree node
      private
      integer :: cut_dim
      ! the dimension to cut
      real(KDTree_RP) :: cut_val
      ! where to cut the dimension
      real(KDTree_RP) :: cut_val_left, cut_val_right
      ! improved cutoffs knowing the spread in child boxes.
      integer :: l, u
      type (KDTree_Node), pointer :: left, right
      type(Interval), pointer :: box(:) => null()
      ! child pointers
      ! Points included in this node are indexes[k] with k \in [l,u]
  end type KDTree_Node

  type KDTree
      ! Global information about the tree, one per tree
      integer :: dimen=0, n=0
      ! dimensionality and total # of points
      real(KDTree_RP), pointer :: the_data(:,:) => null()
      ! pointer to the actual data array
      !
      !  IMPORTANT NOTE:  IT IS DIMENSIONED   the_data(1:d,1:N)
      !  which may be opposite of what may be conventional.
      !  This is, because in Fortran, the memory layout is such that
      !  the first dimension is in sequential order.  Hence, with
      !  (1:d,1:N), all components of the vector will be in consecutive
      !  memory locations.  The search time is dominated by the
      !  evaluation of distances in the terminal nodes.  Putting all
      !  vector components in consecutive memory location improves
      !  memory cache locality, and hence search speed, and may enable
      !  vectorization on some processors and compilers.

      integer, pointer :: ind(:) => null()
      ! permuted index into the data, so that indexes[l..u] of some
      ! bucket represent the indexes of the actual points in that
      ! bucket.
      logical :: sort = .false.
      ! do we always sort output results?
      logical :: rearrange = .false.
      real(KDTree_RP), pointer :: rearranged_data(:,:) => null()
      ! if (rearrange .eqv. .true.) then rearranged_data has been
      ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
      ! permitting search to use more cache-friendly rearranged_data, at
      ! some initial computation and storage cost.
      type (KDTree_Node), pointer :: root => null()
      ! root pointer of the tree
  end type KDTree

  type KDTree_SearchRecord
      !
      ! One of these is created for each search.
      !
      private
      !
      ! Many fields are copied from the tree structure, in order to
      ! speed up the search.
      !
      integer :: dimen
      integer :: nn, nfound
      real(KDTree_RP) :: ballsize
      integer :: centeridx=999, correltime=9999
      ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
      integer :: nalloc     ! how much allocated for results(:)?
      logical :: rearrange  ! are the data rearranged or original?
      ! did the # of points found overflow the storage provided?
      logical :: overflow
      real(KDTree_RP), pointer :: qv(:)          ! query vector
      type(KDTree_Result), pointer :: results(:) ! results
      type(KDTree_PriorityQueue) :: KDTree_PriorityQueue
      real(KDTree_RP), pointer :: dat(:,:) ! temp pointer to data
      integer, pointer :: ind(:)           ! temp pointer to indexes
  end type KDTree_SearchRecord

  ! everything else is private.

  type(KDTree_SearchRecord), save, target :: sr   ! A GLOBAL VARIABLE for search

contains

  function KDTree_Create(input_data, dim, sort, rearrange) result (mr)
    !
    ! create the actual tree structure, given an input array of data.
    !
    ! Note, input data is input_data(1:d,1:N), NOT the other way around.
    ! THIS IS THE REVERSE OF THE PREVIOUS VERSION OF THIS MODULE.
    ! The reason for it is cache friendliness, improving performance.
    !
    ! Optional arguments:  If 'dim' is specified, then the tree
    !                      will only search the first 'dim' components
    !                      of input_data, otherwise, dim is inferred
    !                      from SIZE(input_data,1).
    !
    !                      if sort .eqv. .true. then output results
    !                      will be sorted by increasing distance.
    !                      default=.false., as it is faster to not sort.
    !
    !                      if rearrange .eqv. .true. then an internal
    !                      copy of the data, rearranged by terminal node,
    !                      will be made for cache friendliness.
    !                      default=.true., as it speeds searches, but
    !                      building takes longer, and extra memory is used.
    !
    ! .. Function Return Cut_value ..
    type (KDTree), pointer :: mr
    integer, intent(in), optional :: dim
    logical, intent(in), optional :: sort
    logical, intent(in), optional :: rearrange
    ! ..
    ! .. Array Arguments ..
    real(KDTree_RP), target :: input_data(:,:)
    !
    integer :: i
    ! ..
    allocate (mr)
    mr%the_data => input_data
    ! pointer assignment

    if (present(dim)) then
       mr%dimen = dim
    else
       mr%dimen = size(input_data,1)
    end if
    mr%n = size(input_data,2)

    if (mr%dimen > mr%n) then
       !  unlikely to be correct
       write (*,*) 'KD_TREE_TRANS: likely user error.'
       write (*,*) 'KD_TREE_TRANS: You passed in matrix with D=',mr%dimen
       write (*,*) 'KD_TREE_TRANS: and N=',mr%n
       write (*,*) 'KD_TREE_TRANS: note, that new format is dat(1:D,1:N)'
       write (*,*) 'KD_TREE_TRANS: with usually N >> D.   If N =approx= D, then a k-d tree'
       write (*,*) 'KD_TREE_TRANS: is not an appropriate data structure.'
       stop
    end if

    call BuildTree(mr)

    if (present(sort)) then
       mr%sort = sort
    else
       mr%sort = .false.
    end if

    if (present(rearrange)) then
       mr%rearrange = rearrange
    else
       mr%rearrange = .true.
    end if

    if (mr%rearrange) then
       allocate(mr%rearranged_data(mr%dimen,mr%n))
       do i=1,mr%n
          mr%rearranged_data(:,i) = mr%the_data(:, &
           mr%ind(i))
       end do
    else
       nullify(mr%rearranged_data)
    end if

  end function KDTree_Create


  subroutine BuildTree(tp)
    type (KDTree), pointer :: tp
    ! ..
    integer :: j
    type(KDTree_Node), pointer :: dummy => null()
    ! ..
    allocate (tp%ind(tp%n))
    forall (j=1:tp%n)
       tp%ind(j) = j
    end forall
    tp%root => BuildTreeForRange(tp,1,tp%n, dummy)
  end subroutine BuildTree


  recursive function BuildTreeForRange(tp, l, u ,parent) result (res)
    ! .. Function Return Cut_value ..
    type (KDTree_Node), pointer :: res
    ! ..
    ! .. Structure Arguments ..
    type (KDTree), pointer :: tp
    type (KDTree_Node),pointer           :: parent
    ! ..
    ! .. Scalar Arguments ..
    integer, intent(in)  :: l, u
    ! ..
    ! .. Local Scalars ..
    integer :: i, c, m, dimen
    logical :: recompute
    real(KDTree_RP) :: average

    ! first compute min and max
    dimen = tp%dimen
    allocate (res)
    allocate(res%box(dimen))

    ! First, compute an APPROXIMATE bounding box of all points associated with this node.
    if ( u < l ) then
       ! no points in this box
       nullify(res)
       return
    end if

    if ((u-l)<=bucket_size) then
       !
       ! always compute true bounding box for terminal nodes.
       !
       do i=1,dimen
          call SpreadInCoordinate(tp,i,l,u,res%box(i))
       end do
       res%cut_dim = 0
       res%cut_val = 0.0
       res%l = l
       res%u = u
       res%left =>null()
       res%right => null()
    else
       !
       ! modify approximate bounding box.  This will be an
       ! overestimate of the true bounding box, as we are only recomputing
       ! the bounding box for the dimension that the parent split on.
       !
       ! Going to a true bounding box computation would significantly
       ! increase the time necessary to build the tree, and usually
       ! has only a very small difference.  This box is not used
       ! for searching but only for deciding which coordinate to split on.
       !
       do i=1,dimen
          recompute=.true.
          if (associated(parent)) then
             if (i .ne. parent%cut_dim) then
                recompute=.false.
             end if
          end if
          if (recompute) then
             call SpreadInCoordinate(tp,i,l,u,res%box(i))
          else
             res%box(i) = parent%box(i)
          end if
       end do

       c = maxloc(res%box(1:dimen)%upper-res%box(1:dimen)%lower,1)
       !
       ! c is the identity of which coordinate has the greatest spread.
       !

       if (.false.) then
          ! select exact median to have fully balanced tree.
          m = (l+u)/2
          call SelectOnCoordinate(tp%the_data,tp%ind,c,m,l,u)
       else
          !
          ! select point halfway between min and max, as per A. Moore,
          ! who says this helps in some degenerate cases, or
          ! actual arithmetic average.
          !
          if (.true.) then
             ! actually compute average
             average = sum(tp%the_data(c,tp%ind(l:u))) / real(u-l+1,KDTree_RP)
          else
             average = (res%box(c)%upper + res%box(c)%lower)/2.0
          end if

          res%cut_val = average
          m = SelectOnCoordinateValue(tp%the_data,tp%ind,c,average,l,u)
       end if

       ! moves indexes around
       res%cut_dim = c
       res%l = l
       res%u = u
!      res%cut_val = tp%the_data(c,tp%ind(m))

       res%left => BuildTreeForRange(tp,l,m,res)
       res%right => BuildTreeForRange(tp,m+1,u,res)

       if (associated(res%right) .eqv. .false.) then
          res%box = res%left%box
          res%cut_val_left = res%left%box(c)%upper
          res%cut_val = res%cut_val_left
       elseif (associated(res%left) .eqv. .false.) then
          res%box = res%right%box
          res%cut_val_right = res%right%box(c)%lower
          res%cut_val = res%cut_val_right
       else
          res%cut_val_right = res%right%box(c)%lower
          res%cut_val_left = res%left%box(c)%upper
          res%cut_val = (res%cut_val_left + res%cut_val_right)/2

          ! now remake the true bounding box for self.
          ! Since we are taking unions (in effect) of a tree structure,
          ! this is much faster than doing an exhaustive
          ! search over all points
          res%box%upper = max(res%left%box%upper,res%right%box%upper)
          res%box%lower = min(res%left%box%lower,res%right%box%lower)
       end if
    end if

  end function BuildTreeForRange


  integer function SelectOnCoordinateValue(v, ind, c, alpha, li, ui) &
    result(res)
    ! Move elts of ind around between l and u, so that all points
    ! <= than alpha (in c cooordinate) are first, and then
    ! all points > alpha are second.

    !
    ! Algorithm (matt kennel).
    !
    ! Consider the list as having three parts: on the left,
    ! the points known to be <= alpha.  On the right, the points
    ! known to be > alpha, and in the middle, the currently unknown
    ! points.   The algorithm is to scan the unknown points, starting
    ! from the left, and swapping them so that they are added to
    ! the left stack or the right stack, as appropriate.
    !
    ! The algorithm finishes when the unknown stack is empty.
    !
    ! .. Scalar Arguments ..
    integer, intent(in) :: c, li, ui
    real(KDTree_RP), intent(in) :: alpha
    ! ..
    real(KDTree_RP) :: v(1:,1:)
    integer :: ind(1:)
    integer :: tmp
    ! ..
    integer :: lb, rb
    !
    ! The points known to be <= alpha are in
    ! [l,lb-1]
    !
    ! The points known to be > alpha are in
    ! [rb+1,u].
    !
    ! Therefore we add new points into lb or
    ! rb as appropriate.  When lb=rb
    ! we are done.  We return the location of the last point <= alpha.
    !
    lb = li; rb = ui

    do while (lb < rb)
       if ( v(c,ind(lb)) <= alpha ) then
          ! it is good where it is.
          lb = lb+1
       else
          ! swap it with rb.
          tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
          rb = rb-1
       end if
    end do

    ! now lb .eq. ub
    if (v(c,ind(lb)) <= alpha) then
       res = lb
    else
       res = lb-1
    end if

  end function SelectOnCoordinateValue


  subroutine SelectOnCoordinate(v, ind, c, k, li, ui)
    ! Move elts of ind around between l and u, so that the kth
    ! element
    ! is >= those below, <= those above, in the coordinate c.
    ! .. Scalar Arguments ..
    integer, intent(in)  :: c, k, li, ui
    ! ..
    integer :: i, l, m, s, t, u
    ! ..
    real(KDTree_RP) :: v(:,:)
    integer :: ind(:)
    ! ..
    l = li
    u = ui
    do while (l<u)
       t = ind(l)
       m = l
       do i = l + 1, u
          if (v(c,ind(i))<v(c,t)) then
             m = m + 1
             s = ind(m)
             ind(m) = ind(i)
             ind(i) = s
          end if
       end do
       s = ind(l)
       ind(l) = ind(m)
       ind(m) = s
       if (m<=k) l = m + 1
       if (m>=k) u = m - 1
    end do

  end subroutine SelectOnCoordinate


  subroutine SpreadInCoordinate(tp, c, l, u, interv)
    ! the spread in coordinate 'c', between l and u.
    !
    ! Return lower bound in 'smin', and upper in 'smax',
    ! ..
    ! .. Structure Arguments ..
    type (KDTree), pointer :: tp
    type(Interval), intent(out) :: interv
    ! ..
    ! .. Scalar Arguments ..
    integer, intent(in)  :: c, l, u
    ! ..
    ! .. Local Scalars ..
    real(KDTree_RP) :: last, lmax, lmin, t, smin,smax
    integer :: i, ulocal
    ! ..
    ! .. Local Arrays ..
    real(KDTree_RP), pointer :: v(:,:)
    integer, pointer :: ind(:)
    ! ..
    v => tp%the_data(1:,1:)
    ind => tp%ind(1:)
    smin = v(c,ind(l))
    smax = smin

    ulocal = u

    do i = l + 2, ulocal, 2
       lmin = v(c,ind(i-1))
       lmax = v(c,ind(i))
       if (lmin>lmax) then
          t = lmin
          lmin = lmax
          lmax = t
       end if
       if (smin>lmin) smin = lmin
       if (smax<lmax) smax = lmax
    end do
    if (i==ulocal+1) then
       last = v(c,ind(ulocal))
       if (smin>last) smin = last
       if (smax<last) smax = last
    end if

    interv%lower = smin
    interv%upper = smax

  end subroutine SpreadInCoordinate


  subroutine KDTree_Destroy(tp)
    ! Deallocates all memory for the tree, except input data matrix
    ! .. Structure Arguments ..
    type (KDTree), pointer :: tp
    ! ..
    call DestroyNode(tp%root)

    deallocate(tp%ind)
    nullify(tp%ind)

    if (tp%rearrange) then
       deallocate(tp%rearranged_data)
       nullify(tp%rearranged_data)
    end if

    deallocate(tp)

  contains

    recursive subroutine DestroyNode(np)
      ! .. Structure Arguments ..
      type (KDTree_Node), pointer :: np
      ! ..
      ! .. Intrinsic Functions ..
      intrinsic associated
      ! ..
      if (associated(np%left)) then
         call DestroyNode(np%left)
         nullify(np%left)
      end if
      if (associated(np%right)) then
         call DestroyNode(np%right)
         nullify (np%right)
      end if
      if (associated(np%box)) deallocate(np%box)
      deallocate(np)
      return

    end subroutine DestroyNode

  end subroutine KDTree_Destroy


  subroutine KDTree_N_Nearest(tp, qv, nn, results)
    ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
    ! returning their indexes and distances in 'indexes' and 'distances'
    ! arrays already allocated passed to this subroutine.
    type (KDTree), pointer              :: tp
    real(KDTree_RP), target, intent(in) :: qv(:)
    integer, intent(in)                 :: nn
    type(KDTree_Result), target         :: results(:)

    sr%ballsize = huge(1.0)
    sr%qv => qv
    sr%nn = nn
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false.

    sr%results => results

    sr%nalloc = nn   ! will be checked

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%dat => tp%rearranged_data
    else
       sr%dat => tp%the_data
    end if
    sr%dimen = tp%dimen

    call ValidateQueryStorage(nn)
    sr%KDTree_PriorityQueue = KDTree_PQ_Create(results)

    call Search(tp%root)

    if (tp%sort) then
       call KDTree_SortResults(nn, results)
    end if

  end subroutine KDTree_N_Nearest


  subroutine KDTree_N_NearestAroundPoint(tp, idxin, correltime, nn, results)
    ! Find the 'nn' vectors in the tree nearest to point 'idxin',
    ! with correlation window 'correltime', returing results in
    ! results(:), which must be pre-allocated upon entry.
    type (KDTree), pointer        :: tp
    integer, intent(in)            :: idxin, correltime, nn
    type(KDTree_Result), target   :: results(:)

    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin)       ! copy the vector
    sr%ballsize = huge(1.0_KDTree_RP)  ! the largest real number
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%nn = nn
    sr%nfound = 0

    sr%dimen = tp%dimen
    sr%nalloc = nn

    sr%results => results

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%dat => tp%rearranged_data
    else
       sr%dat => tp%the_data
    end if

    call ValidateQueryStorage(nn)
    sr%KDTree_PriorityQueue = KDTree_PQ_Create(results)

    call Search(tp%root)

    if (tp%sort) then
       call KDTree_SortResults(nn, results)
    end if
    deallocate (sr%qv)
    return

  end subroutine KDTree_N_NearestAroundPoint


  subroutine KDTree_R_Nearest(tp, qv, r2, nfound, nalloc, results)
    ! find the nearest neighbors to point 'idxin', within SQUARED
    ! Euclidean distance 'r2'.   Upon ENTRY, nalloc must be the
    ! size of memory allocated for results(1:nalloc).  Upon
    ! EXIT, nfound is the number actually found within the ball.
    !
    !  Note that if nfound .gt. nalloc then more neighbors were found
    !  than there were storage to store.  The resulting list is NOT
    !  the smallest ball inside norm r^2
    !
    ! Results are NOT sorted unless tree was created with sort option.

    type (KDTree), pointer              :: tp
    real(KDTree_RP), target, intent(in) :: qv(:)
    real(KDTree_RP), intent(in)         :: r2
    integer, intent(out)                :: nfound
    integer, intent(in)                 :: nalloc
    type(KDTree_Result), target         :: results(:)

    sr%qv => qv
    sr%ballsize = r2
    sr%nn = 0      ! flag for fixed ball search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0

    sr%results => results

    call ValidateQueryStorage(nalloc)
    sr%nalloc = nalloc
    sr%overflow = .false.
    sr%ind => tp%ind
    sr%rearrange= tp%rearrange

    if (tp%rearrange) then
       sr%dat => tp%rearranged_data
    else
       sr%dat => tp%the_data
    end if
    sr%dimen = tp%dimen

    !
    !sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call Search(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call KDTree_SortResults(nfound, results)
    end if

    if (sr%overflow) then
       write (*,*) 'KD_TREE_TRANS: warning! return from KDTree_R_Nearest found more neighbors'
       write (*,*) 'KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball'
       write (*,*) 'KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.'
    end if

  end subroutine KDTree_R_Nearest


  subroutine KDTree_R_NearestAroundPoint(tp, idxin, correltime, r2, &
                                         nfound, nalloc, results    )
    !
    ! Like KDTree_R_Nearest, but around a point 'idxin' already existing
    ! in the data set.
    !
    ! Results are NOT sorted unless tree was created with sort option.
    !
    type (KDTree), pointer      :: tp
    integer, intent(in)         :: idxin, correltime, nalloc
    real(KDTree_RP), intent(in) :: r2
    integer, intent(out)        :: nfound
    type(KDTree_Result), target :: results(:)
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic huge
    ! ..
    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin) ! copy the vector
    sr%ballsize = r2
    sr%nn = 0    ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime

    sr%results => results

    sr%nalloc = nalloc
    sr%overflow = .false.

    call ValidateQueryStorage(nalloc)

    !    sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !    sr%il = -1               ! set to invalid indexes

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (tp%rearrange) then
       sr%dat => tp%rearranged_data
    else
       sr%dat => tp%the_data
    end if
    sr%rearrange = tp%rearrange
    sr%dimen = tp%dimen

    !
    !sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !

    call Search(tp%root)
    nfound = sr%nfound
    if (tp%sort) then
       call KDTree_SortResults(nfound,results)
    end if

    if (sr%overflow) then
       write (*,*) 'KD_TREE_TRANS: warning! return from KDTree_R_Nearest found more neighbors'
       write (*,*) 'KD_TREE_TRANS: than storage was provided for.  Answer is NOT smallest ball'
       write (*,*) 'KD_TREE_TRANS: with that number of neighbors!  I.e. it is wrong.'
    end if

    deallocate (sr%qv)

  end subroutine KDTree_R_NearestAroundPoint


  function KDTree_R_Count(tp, qv, r2) result(nfound)
    ! Count the number of neighbors within square distance 'r2'.
    type (KDTree), pointer              :: tp
    real(KDTree_RP), target, intent(in) :: qv(:)
    real(KDTree_RP), intent(in)         :: r2
    integer                             :: nfound
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic huge
    ! ..
    sr%qv => qv
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0

    nullify(sr%results) ! for some reason, FTN 95 chokes on '=> null()'

    sr%nalloc = 0            ! we do not allocate any storage but that's OK
                             ! for counting.
    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%dat => tp%rearranged_data
    else
       sr%dat => tp%the_data
    end if
    sr%dimen = tp%dimen

    !
    !sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call Search(tp%root)

    nfound = sr%nfound

  end function KDTree_R_Count


  function KDTree_R_CountAroundPoint(tp, idxin, correltime, r2) result(nfound)
    ! Count the number of neighbors within square distance 'r2' around
    ! point 'idxin' with decorrelation time 'correltime'.
    !
    type (KDTree), pointer      :: tp
    integer, intent(in)         :: correltime, idxin
    real(KDTree_RP), intent(in) :: r2
    integer                     :: nfound
    ! ..
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic huge
    ! ..
    allocate (sr%qv(tp%dimen))
    sr%qv = tp%the_data(:,idxin)
    sr%ballsize = r2

    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    nullify(sr%results)

    sr%nalloc = 0            ! we do not allocate any storage but that's OK
                             ! for counting.

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange

    if (sr%rearrange) then
       sr%dat => tp%rearranged_data
    else
       sr%dat => tp%the_data
    end if
    sr%dimen = tp%dimen

    !
    !sr%dsl = huge(sr%dsl)    ! set to huge positive values
    !sr%il = -1               ! set to invalid indexes
    !
    sr%overflow = .false.

    call Search(tp%root)

    nfound = sr%nfound

  end function KDTree_R_CountAroundPoint


  subroutine ValidateQueryStorage(n)
    !
    ! make sure we have enough storage for n
    !
    integer, intent(in) :: n

    if (size(sr%results,1) .lt. n) then
       write (*,*) 'KD_TREE_TRANS:  you did not provide enough storage for results(1:n)'
       stop
       return
    end if

  end subroutine ValidateQueryStorage


  function SquareDistance(d, iv, qv) result (res)
    ! distance between iv[1:n] and qv[1:n]
    ! .. Function Return Value ..
    ! re-implemented to improve vectorization.
    real(KDTree_RP) :: res
    ! ..
    ! ..
    ! .. Scalar Arguments ..
    integer :: d
    ! ..
    ! .. Array Arguments ..
    real(KDTree_RP) :: iv(:),qv(:)
    ! ..
    ! ..
    res = sum( (iv(1:d)-qv(1:d))**2 )
  end function SquareDistance


  recursive subroutine Search(node)
    !
    ! This is the innermost core routine of the kd-tree search.  Along
    ! with "ProcessTerminalNode", it is the performance bottleneck.
    !
    ! This version uses a logically complete secondary search of
    ! "box in bounds", whether the sear
    !
    type (KDTree_Node), pointer :: node
    ! ..
    type(KDTree_Node), pointer  :: ncloser, nfarther
    !
    integer                     :: cut_dim, i
    ! ..
    real(KDTree_RP)             :: qval, dis
    real(KDTree_RP)             :: ballsize
    real(KDTree_RP), pointer    :: qv(:)
    type(Interval), pointer     :: box(:)

    if ((associated(node%left) .and. associated(node%right)) .eqv. .false.) then
       ! we are on a terminal node
       if (sr%nn .eq. 0) then
          call ProcessTerminalNodeFixedBall(node)
       else
          call ProcessTerminalNode(node)
       end if
    else
       ! we are not on a terminal node
       qv => sr%qv(1:)
       cut_dim = node%cut_dim
       qval = qv(cut_dim)

       if (qval < node%cut_val) then
          ncloser => node%left
          nfarther => node%right
          dis = (node%cut_val_right - qval)**2
!          extra = node%cut_val - qval
       else
          ncloser => node%right
          nfarther => node%left
          dis = (node%cut_val_left - qval)**2
!          extra = qval- node%cut_val_left
       end if

       if (associated(ncloser)) call Search(ncloser)

       ! we may need to search the second node.
       if (associated(nfarther)) then
          ballsize = sr%ballsize
!          dis=extra**2
          if (dis <= ballsize) then
             !
             ! we do this separately as going on the first cut dimen is often
             ! a good idea.
             ! note that if extra**2 < sr%ballsize, then the next
             ! check will also be false.
             !
             box => node%box(1:)
             do i=1,sr%dimen
                if (i .ne. cut_dim) then
                   dis = dis + Dis2FromBnd(qv(i),box(i)%lower,box(i)%upper)
                   if (dis > ballsize) then
                      return
                   end if
                end if
             end do

             !
             ! if we are still here then we need to search mroe.
             !
             call Search(nfarther)
          end if
       end if
    end if

  end subroutine Search


  real(KDTree_RP) function Dis2FromBnd(x, amin, amax) result(res)
    real(KDTree_RP), intent(in) :: x, amin,amax

    if (x > amax) then
       res = (x-amax)**2;
       return
    else
       if (x < amin) then
          res = (amin-x)**2;
          return
       else
          res = 0.0
          return
       end if
    end if

  end function Dis2FromBnd


  subroutine ProcessTerminalNode(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure.
    !
    type (KDTree_Node), pointer :: node
    !
    real(KDTree_RP), pointer    :: qv(:)
    integer, pointer            :: ind(:)
    real(KDTree_RP), pointer    :: dat(:,:)
    !
    integer          :: dimen, i, indexofi, k, centeridx, correltime
    real(KDTree_RP)  :: ballsize, sd, newpri
    logical          :: rearrange
    type(KDTree_PriorityQueue), pointer :: priority_queue
    !
    ! copy values from sr to local variables
    !
    !
    ! Notice, making local pointers with an EXPLICIT lower bound
    ! seems to generate faster code.
    ! why?  I don't know.
    qv => sr%qv(1:)
    priority_queue => sr%KDTree_PriorityQueue
    dimen = sr%dimen
    ballsize = sr%ballsize
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    dat => sr%dat(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime

    !    doing_correl = (centeridx >= 0)  ! Do we have a decorrelation window?
    !    include_point = .true.    ! by default include all points
    ! search through terminal bucket.

    mainloop: do i = node%l, node%u
       if (rearrange) then
          sd = 0.0
          do k = 1,dimen
             sd = sd + (dat(k,i) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0
          do k = 1,dimen
             sd = sd + (dat(k,indexofi) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
       end if

       if (centeridx > 0) then ! doing correlation Interval?
          if (abs(indexofi-centeridx) < correltime) cycle mainloop
       end if


       !
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity.
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor.
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one.

       if (sr%nfound .lt. sr%nn) then
          !
          ! add this point unconditionally to fill list.
          !
          sr%nfound = sr%nfound +1
          newpri = KDTree_PQ_Insert(priority_queue,sd,indexofi)
          if (sr%nfound .eq. sr%nn) ballsize = newpri
          ! we have just filled the working list.
          ! put the best square distance to the maximum value
          ! on the list, which is extractable from the PQ.

       else
          !
          ! now, if we get here,
          ! we know that the current node has a squared
          ! distance smaller than the largest one on the list, and
          ! belongs on the list.
          ! Hence we replace that with the current one.
          !
          ballsize = KDTree_PQ_ReplaceMax(priority_queue,sd,indexofi)
       end if
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%ballsize = ballsize

  end subroutine ProcessTerminalNode


  subroutine ProcessTerminalNodeFixedBall(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure, i.e.
    ! save all within a fixed ball.
    !
    type (KDTree_Node), pointer :: node
    !
    real(KDTree_RP), pointer    :: qv(:)
    integer, pointer            :: ind(:)
    real(KDTree_RP), pointer    :: dat(:,:)
    !
    integer                     :: nfound
    integer                     :: dimen, i, indexofi, k
    integer                     :: centeridx, correltime, nn
    real(KDTree_RP)             :: ballsize, sd
    logical                     :: rearrange

    !
    ! copy values from sr to local variables
    !
    qv => sr%qv(1:)
    dimen = sr%dimen
    ballsize = sr%ballsize
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    dat => sr%dat(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime
    nn = sr%nn ! number to search for
    nfound = sr%nfound

    ! search through terminal bucket.
    mainloop: do i = node%l, node%u

       !
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity.
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor.
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one.

       ! which index to the point do we use?

       if (rearrange) then
          sd = 0.0
          do k = 1,dimen
             sd = sd + (dat(k,i) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = 0.0
          do k = 1,dimen
             sd = sd + (dat(k,indexofi) - qv(k))**2
             if (sd>ballsize) cycle mainloop
          end do
       end if

       if (centeridx > 0) then ! doing correlation Interval?
          if (abs(indexofi-centeridx)<correltime) cycle mainloop
       end if

       nfound = nfound+1
       if (nfound .gt. sr%nalloc) then
          ! oh nuts, we have to add another one to the tree but
          ! there isn't enough room.
          sr%overflow = .true.
       else
          sr%results(nfound)%dis = sd
          sr%results(nfound)%idx = indexofi
       end if
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%nfound = nfound
  end subroutine ProcessTerminalNodeFixedBall


  subroutine KDTree_N_NearestBruteForce(tp, qv, nn, results)
    ! find the 'n' nearest neighbors to 'qv' by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type(KDTree), pointer       :: tp
    real(KDTree_RP), intent(in) :: qv(:)
    integer, intent(in)         :: nn
    type(KDTree_Result)         :: results(:)

    integer :: i, j, k
    real(KDTree_RP), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = SquareDistance(tp%dimen,qv,tp%the_data(:,i))
    end do
    ! now find 'n' smallest distances
    do i = 1, nn
       results(i)%dis =  huge(1.0)
       results(i)%idx = -1
    end do
    do i = 1, tp%n
       if (all_distances(i)<results(nn)%dis) then
          ! insert it somewhere on the list
          do j = 1, nn
             if (all_distances(i)<results(j)%dis) exit
          end do
          ! now we know 'j'
          do k = nn - 1, j, -1
             results(k+1) = results(k)
          end do
          results(j)%dis = all_distances(i)
          results(j)%idx = i
       end if
    end do
    deallocate (all_distances)
  end subroutine KDTree_N_NearestBruteForce


  subroutine KDTree_R_NearestBruteForce(tp,qv,r2,nfound,results)
    ! find the nearest neighbors to 'qv' with distance**2 <= r2 by exhaustive search.
    ! only use this subroutine for testing, as it is SLOW!  The
    ! whole point of a k-d tree is to avoid doing what this subroutine
    ! does.
    type(KDTree), pointer       :: tp
    real(KDTree_RP), intent(in) :: qv(:)
    real(KDTree_RP), intent(in) :: r2
    integer, intent(out)        :: nfound
    type(KDTree_Result)         :: results(:)

    integer :: i, nalloc
    real(KDTree_RP), allocatable :: all_distances(:)
    ! ..
    allocate (all_distances(tp%n))
    do i = 1, tp%n
       all_distances(i) = SquareDistance(tp%dimen,qv,tp%the_data(:,i))
    end do

    nfound = 0
    nalloc = size(results,1)

    do i = 1, tp%n
       if (all_distances(i)< r2) then
          ! insert it somewhere on the list
          if (nfound .lt. nalloc) then
             nfound = nfound+1
             results(nfound)%dis = all_distances(i)
             results(nfound)%idx = i
          end if
       end if
    end do
    deallocate (all_distances)

    call KDTree_SortResults(nfound,results)

  end subroutine KDTree_R_NearestBruteForce


  subroutine KDTree_SortResults(nfound, results)
    !  Use after search to sort results(1:nfound) in order of increasing
    !  distance.
    integer, intent(in)         :: nfound
    type(KDTree_Result), target :: results(:)
    !

    if (nfound .gt. 1) call HeapsortStruct(results, nfound)

    return
  end subroutine KDTree_SortResults


  subroutine HeapsortStruct(a, n)
    !
    ! Sort a(1:n) in ascending order
    !
    !
    integer, intent(in)                :: n
    type(KDTree_Result), intent(inout) :: a(:)

    !
    !
    type(KDTree_Result) :: value ! temporary value

    integer :: i,j
    integer :: ileft, iright

    ileft=n/2+1
    iright=n

    if(n.eq.1) return

    do
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft)
       else
          value=a(iright)
          a(iright)=a(1)
          iright=iright-1
          if (iright == 1) then
             a(1) = value
             return
          end if
       end if
       i=ileft
       j=2*ileft
       do while (j <= iright)
          if(j < iright) then
             if(a(j)%dis < a(j+1)%dis) j=j+1
          end if
          if(value%dis < a(j)%dis) then
             a(i)=a(j);
             i=j
             j=j+j
          else
             j=iright+1
          end if
       end do
       a(i)=value
    end do
  end subroutine HeapsortStruct


end module KD_Tree
