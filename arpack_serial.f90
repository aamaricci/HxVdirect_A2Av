subroutine lanczos_arpack_d(MatVec,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose)
  !Interface to Matrix-Vector routine:
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  !Arguments
  integer                   :: Neigen
  integer                   :: Nblock
  integer                   :: Nitermax
  real(8)                   :: eval(Neigen)
  real(8)                   :: evec(:,:) !nloc,Neigen
  character(len=2),optional :: which
  real(8),optional          :: v0(size(evec,1))
  real(8),optional          :: tol
  logical,optional          :: iverbose
  !Dimensions:
  integer                   :: ldv
  integer                   :: n,nconv,ncv,nev
  !Arrays:
  real(8),allocatable       :: ax(:),d(:,:)
  real(8),allocatable       :: resid(:)
  real(8),allocatable       :: workl(:),workd(:)
  real(8),allocatable       :: v(:,:)
  logical,allocatable       :: select(:)
  integer                   :: iparam(11)
  integer                   :: ipntr(11)
  !Control Vars:
  integer                   :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
  logical                   :: rvec,verb
  integer                   :: i
  real(8)                   :: sigma
  real(8)                   :: tol_
  character                 :: bmat  
  character(len=2)          :: which_
  real(8),external          :: dnrm2
  !
  which_='SA';if(present(which))which_=which
  tol_=0d0;if(present(tol))tol_=tol
  verb=.false.;if(present(iverbose))verb=iverbose
  !
  ldv    = size(evec,1)
  !
  n      = ldv
  nev    = Neigen
  ncv    = Nblock
  !
  if(ncv>ldv)stop "ARPACK error: Ncv > Ldv"
  !
  bmat   = 'I'
  maxitr = Nitermax
  ! 
  allocate(ax(n))
  allocate(d(ncv,2))
  allocate(resid(n))
  allocate(workl(ncv*(ncv+8)))
  allocate(workd(3*n))
  allocate(v(ldv,ncv))
  allocate(select(ncv))
  ax     =0d0
  d      =0d0
  resid  =0d0
  workl  =0d0
  workd  =0d0
  v      =0d0
  lworkl = ncv*(ncv+8)
  info   = 1
  ido    = 0
  ishfts    = 1
  mode1     = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode1
  if(present(v0))then
     resid=v0
  else
     call mt_random(resid)
  endif
  resid=resid/sqrt(dot_product(resid,resid))
  !
  !MAIN LOOP, REVERSE COMMUNICATION
  do
     call dsaupd(ido,bmat,n,which_,nev,tol_,resid,ncv,v,ldv,&
          iparam,ipntr,workd,workl,lworkl,info)
     if(ido/=-1.AND.ido/=1)then
        exit
     end if
     !  Perform matrix vector multiplication
     !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
     call MatVec(N,workd(ipntr(1)),workd(ipntr(2)) )
  end do
  !
  !POST PROCESSING:
  if(info<0)then
     write(*,'(a,i6)')'Error in DSAUPD, info = ', info
     stop
  else
     rvec = .true.
     call dseupd(rvec,'All',select,d,v,ldv,sigma,bmat,n,which_,&
          nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
     !
     do j=1,neigen
        eval(j)=d(j,1)
     enddo
     !
     if(any(shape(evec)/=[ldv,Neigen]))stop "arpack_mpi ERROR: evec has wrong dimension. Reduce=T"
     !
     do j=1,neigen
        do i=1,ldv
           evec(i,j)=v(i,j)
        enddo
     enddo
     !
     !  Compute the residual norm
     !    ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and 
     !  eigenvectors.  (iparam(5) indicates how many are 
     !  accurate to the requested tolerance)
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with DSEUPD (get Evec), ierr = ',ierr
        ! else
        !    nconv =  iparam(5)
        !    do j = 1, nconv
        !       call MatVec(n, v(1,j), ax )
        !       call daxpy( n, -d(j,1), v(1,j), 1, ax, 1 )
        !       d(j,2) = dnrm2(n,ax,1)
        !       d(j,2) = d(j,2) / abs ( d(j,1) )
        !    end do
        !    if(verb)call dmout(6,nconv,2,d,maxncv,-6,'Ritz values and relative residuals')
     end if
     !
     if(info==1) then
        write(*,'(a)' ) ' '
        write(*,'(a)' ) '  Maximum number of iterations reached.'
     elseif(info==3) then
        write(*,'(a)' ) ' '
        write(*,'(a)' ) '  No shifts could be applied during implicit '&
             //'Arnoldi update, try increasing NCV.'
     end if
  endif
end subroutine lanczos_arpack_d
