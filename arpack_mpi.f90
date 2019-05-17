subroutine lanczos_parpack_d(MpiComm,MatVec,Neigen,Nblock,Nitermax,eval,evec,which,v0,tol,iverbose)
  !Arguments
  integer                      :: MpiComm
  !Interface to Matrix-Vector routine:
  interface
     subroutine MatVec(nchunk,vin,vout)
       integer                   :: nchunk
       real(8),dimension(nchunk) :: vin,vout
     end subroutine MatVec
  end interface
  !Arguments
  integer                      :: Neigen
  integer                      :: Nblock
  integer                      :: Nitermax
  real(8)                      :: eval(Neigen)
  real(8)                      :: evec(:,:) ![Nloc,Neigen]
  character(len=2),optional    :: which
  real(8),optional             :: v0(size(evec,1))
  real(8),optional             :: tol
  logical,optional             :: iverbose
  !Dimensions:
  integer                      :: ldv
  integer                      :: n,nconv,ncv,nev
  !Arrays:
  real(8),allocatable          :: evec_tmp(:) ![Nloc] see above
  real(8),allocatable          :: ax(:),d(:,:)
  real(8),allocatable          :: resid(:),vec(:)
  real(8),allocatable          :: workl(:),workd(:)
  real(8),allocatable          :: v(:,:)
  logical,allocatable          :: select(:)
  integer                      :: iparam(11)
  integer                      :: ipntr(11)
  !Control Vars:
  integer                      :: ido,ierr,info,ishfts,j,lworkl,maxitr,mode1
  logical                      :: verb
  integer                      :: i
  real(8)                      :: sigma
  real(8)                      :: tol_,norm_tmp,norm
  character                    :: bmat  
  character(len=2)             :: which_
  real(8),external             :: dnrm2
  !MPI
  integer                      :: irank,mpiRank,mpiSize,seed  
  real(8),allocatable          :: MpiSeed(:)
  logical                      :: mpi_master
  !
  which_ = 'SA'   ; if(present(which))which_=which
  tol_   = 0d0    ; if(present(tol))tol_=tol
  verb   = .false.; if(present(iverbose))verb=iverbose
  !
  !MPI setup:
  mpiRank   = get_Rank_MPI(MpiComm)
  mpiSize   = get_Size_MPI(MpiComm)
  mpi_master= Get_master_MPI(MpiComm)
  !
  ! Setup distribution of data to nodes:
  ldv = size(evec,1)            !ldv is the SMALL dimension
  !
  n   = ldv
  nev = Neigen
  ncv = Nblock
  !
  if(ncv>ldv) stop "PARPACK error: Ncv > Ldv"
  !
  bmat   = 'I'
  maxitr = Nitermax
  !=========================================================================
  !
  allocate(ax(ldv))
  allocate(resid(ldv))
  allocate(workd(3*ldv))
  allocate(v(ldv,ncv))
  allocate(d(ncv,2))
  allocate(workl(ncv*(ncv+8)))
  allocate(select(ncv))
  !
  ax     =0.d0
  d      =0.d0
  resid  =0.d0
  workl  =0.d0
  workd  =0.d0
  v      =0.d0
  lworkl = ncv*(ncv+8)
  info   = 1
  ido    = 0
  !
  if(present(v0))then
     resid = v0
  else
     allocate(MpiSeed(MpiSize))
     call random_number(MpiSeed)
     do irank=0,MpiSize-1
        if(irank==MpiRank)seed = int(123456789*MpiSeed(irank+1))
     end do

     call mersenne_init(seed)
     call mt_random(resid)
     norm_tmp=dot_product(resid,resid)
     call AllReduce_MPI(MpiComm,norm_tmp,norm)
     resid=resid/sqrt(norm)
     ! resid = 1d0/sqrt(dble(Ns)) !assume the starting vector is [1,1,1,1,....]
  endif
  !
  ishfts    = 1
  mode1     = 1
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode1
  !
  !  MAIN LOOP (Reverse communication loop)
  do
     call pdsaupd(MpiComm,ido,bmat,ldv,which_,nev,tol_,resid,&
          ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info )
     if(ido/=-1.AND.ido/=1)exit
     !  Perform matrix vector multiplication
     !    y <--- OP*x ; workd(ipntr(1))=input, workd(ipntr(2))=output
     call MatVec(ldv,workd(ipntr(1)),workd(ipntr(2)))
  end do
  !
  if(info<0)then
     if(mpi_master)write(*,'(a,i6)')'Fatal Error in PDSAUPD, info = ', info
     stop
  else
     call pdseupd (MpiComm,.true.,'All',select,d,v,ldv,sigma,bmat,&
          ldv,which_,nev,tol_,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,ierr)
     do j=1,neigen
        eval(j)=d(j,1)
     enddo
     !
     if(any(shape(evec)/=[ldv,Neigen]))stop "arpack_mpi ERROR: evec has wrong dimension. Reduce=T"
     evec=0d0
     do j=1,neigen
        do i=1,ldv
           evec(i,j)=v(i,j)
        enddo
     enddo
     !
     !
     !=========================================================================
     !  Compute the residual norm
     !    ||  A*x - lambda*x ||
     !  for the NCONV accurately computed eigenvalues and 
     !  eigenvectors.  (iparam(5) indicates how many are 
     !  accurate to the requested tolerance)
     if(ierr/=0)then
        write(*,'(a,i6)')'Error with PDSEUPD, IERR = ',ierr
        write(*,'(a)')'Check the documentation of SSEUPD.'
     else
        nconv =  iparam(5)
        !
        ! do j = 1, nconv
        !    call hprod( n, v(1,j), ax )
        !    call zaxpy( n, -d(j), v(1,j), 1, ax, 1 )
        !    rd(j,1) = dble (d(j))
        !    rd(j,2) = dimag (d(j))
        !    rd(j,3) = dznrm2 (n, ax, 1)
        !    rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
        ! end do
        ! if(mpi_master.AND.verb)call dmout(6,nconv,3,rd,maxncv,-6,'Ritz values and relative residuals')
     end if

     if(mpi_master.AND.verb)then
        write(*,'(a)') ''
        write(*,'(a)') 'ARPACK::'
        write(*,'(a)') ''
        write(*,'(a,i12)') '  Size of the matrix is:                      ', n
        write(*,'(a,i6)') '  Number of Ritz values requested is:         ', nev
        write(*,'(a,i6)') '  Number of Arnoldi vectors generated is:     ', ncv
        write(*,'(a)')    '  Portion of the spectrum:                        '//trim(which_)
        write(*,'(a,i6)') '  Number of converged Ritz values is:         ', nconv
        write(*,'(a,i6)') '  Number of Implicit Arnoldi iterations is:   ', iparam(3)
        write(*,'(a,i6)') '  Number of OP*x is:                          ', iparam(9)
        write(*,'(a,ES14.6)') '  The convergence criterion is:           ', tol
     end if

     if(mpi_master)then       
        if(info==1) then
           write(*,'(a)' ) ' '
           write(*,'(a)' ) '  Maximum number of iterations reached.'
        elseif(info==3) then
           write(*,'(a)' ) ' '
           write(*,'(a)' ) '  No shifts could be applied during implicit '&
                //'Arnoldi update, try increasing NCV.'
        end if
     endif
  endif
  deallocate(ax,resid,workd,v,d,workl,select)
end subroutine lanczos_parpack_d









