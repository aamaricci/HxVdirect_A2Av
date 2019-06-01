module SF_SP_LINALG_
  USE SF_LINALG, only: eye,eigh
  USE SF_RANDOM, only: mt_random,mersenne_init
  USE SF_MPI
  USE MPI
  implicit none
  private 

  interface sp_eigh
     module procedure :: lanczos_parpack_d
  end interface sp_eigh


  interface sp_lanc_eigh
     module procedure :: mpi_lanczos_eigh_d
  end interface sp_lanc_eigh




  complex(8),parameter              :: zero=(0d0,0d0)
  complex(8),parameter              :: one=(1d0,0d0)
  complex(8),parameter              :: xi=(0d0,1d0)
  integer,allocatable               :: seed_random(:)
  integer                           :: nrandom
  logical                           :: verb=.false.
  real(8)                           :: threshold_=1.d-12
  integer                           :: ncheck_=10
  integer                           :: nrestart_=10



  !****************************************************************************
  !                                      PUBLIC 
  !****************************************************************************
  public :: sp_eigh
  public :: sp_lanc_eigh
  !****************************************************************************


contains


  !##################################################################
  ! ARPACK METHOD for LOWEST part of the spectrum of a of a SPARSE
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
    integer :: count=0
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
       ! allocate(MpiSeed(MpiSize))
       ! call random_number(MpiSeed)
       do irank=0,MpiSize-1
          ! if(irank==MpiRank)seed = int(123456789*MpiSeed(irank+1))
          if(irank==MpiRank)seed = int(123456780 + irank)
       end do
       ! deallocate(MpiSeed)
       call mersenne_init(seed)
       call mt_random(resid)
       norm_tmp=dot_product(resid,resid)
       call AllReduce_MPI(MpiComm,norm_tmp,norm)
       resid=resid/sqrt(norm)
       !
       ! resid = 1d0
       ! norm_tmp=dble(size(resid))
       ! call AllReduce_MPI(MpiComm,norm_tmp,norm)
       ! resid=resid/sqrt(norm)
       !resid = 1d0/sqrt(dble(Ns)) !assume the starting vector is [1,1,1,1,....]
       !
       if(verb.AND.mpi_master)write(*,*)"MPI_ARPACK_EIGH: random initial vector generated"
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














  !##################################################################
  ! LANCZOS METHOD for LOWEST EigenSolution OR tri-diagonalization
  !---------------------------------------------------------------------
  !Purpose: use plain lanczos to get the groundstate energy
  !---------------------------------------------------------------------
  subroutine mpi_lanczos_eigh_d(MpiComm,MatVec,Nitermax,Egs,Vect,iverbose,threshold,ncheck)
    integer                              :: MpiComm
    interface 
       subroutine MatVec(Nloc,vin,vout)
         integer                         :: Nloc
         real(8)                         :: vin(Nloc)
         real(8)                         :: vout(Nloc)
       end subroutine MatVec
    end interface
    integer                              :: Nitermax
    integer                              :: Nloc
    real(8)                              :: egs,aa,bb
    real(8),dimension(:)                 :: vect !Nloc
    real(8),optional                     :: threshold
    integer,optional                     :: ncheck
    logical,optional                     :: iverbose
    !
    real(8),dimension(size(vect))        :: vin,vout
    real(8),dimension(:),allocatable     :: MpiSeed
    integer                              :: iter,nlanc
    real(8),dimension(Nitermax+1)        :: alanc,blanc
    real(8),dimension(Nitermax,Nitermax) :: Z
    real(8),dimension(Nitermax)          :: diag,subdiag
    real(8)                              :: esave,esave_prev,diff
    real(8)                              :: a_,b_,norm,norm_tmp
    integer                              :: i,ierr,irank,seed
    !
    integer :: MpiRank,MpiSize
    logical                              :: mpi_master
    !
    MpiRank = get_Rank_mpi(MpiComm)
    MpiSize = get_Size_mpi(MpiComm)
    mpi_master=get_master_MPI(MpiComm)
    !
    Nloc = size(vect)
    !
    if(present(iverbose))verb=iverbose
    if(present(threshold))threshold_=threshold
    if(present(ncheck))ncheck_=ncheck
    !
    norm_tmp=dot_product(vect,vect)
    call AllReduce_MPI(MpiComm,norm_tmp,norm)
    !
    if(norm==0d0)then
       ! allocate(MpiSeed(MpiSize))
       ! call random_number(MpiSeed)
       do irank=0,MpiSize-1
          ! if(irank==MpiRank)seed = int(123456789*MpiSeed(irank+1))
          if(irank==MpiRank)seed = int(123456789 + irank)
       end do
       ! deallocate(MpiSeed)
       call mersenne_init(seed)
       call mt_random(vect)
       norm_tmp=dot_product(vect,vect)
       call AllReduce_MPI(MpiComm,norm_tmp,norm)
       vect=vect/sqrt(norm)
       if(verb.AND.mpi_master)write(*,*)"MPI_LANCZOS_EIGH: random initial vector generated"
    endif
    !
    !============= LANCZOS LOOP =====================
    !
    Vin   = Vect                   !save input vector for Eigenvector calculation:
    Vout  = 0d0
    esave = 0d0
    diff  = 0d0
    alanc = 0d0
    blanc = 0d0
    nlanc = 0
    lanc_loop: do iter=1,Nitermax
       call mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,a_,b_)
       if(abs(b_)<threshold_)exit lanc_loop
       nlanc=nlanc+1
       alanc(iter) = a_ ; blanc(iter+1) = b_
       !
       diag(1:Nlanc)    = alanc(1:Nlanc)
       subdiag(2:Nlanc) = blanc(2:Nlanc)
       call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
       !
       esave_prev = esave
       esave      = diag(1)
       diff       = esave - esave_prev
       !
       if( nlanc >= Ncheck_ )then
          if(verb.AND.mpi_master)write(*,*)'Iter, E0, deltaE = ',iter,diag(1),diff
          if(abs(diff) <= threshold_)exit lanc_loop
       endif
       !
    enddo lanc_loop
    if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
    !
    !
    !============== END LANCZOS LOOP ======================
    !
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    !Get the Eigenvalues:
    egs = diag(1)
    !
    !Get the Eigenvector:
    Vin = Vect
    Vout= 0d0
    Vect= 0d0
    do iter=1,nlanc
       call mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,alanc(iter),blanc(iter))
       Vect = Vect + vin*Z(iter,1)
    end do
    norm_tmp = dot_product(vect,vect)
    call Allreduce_MPI(MpiComm,norm_tmp,norm)
    vect=vect/sqrt(norm)
    if(verb)then
       call MatVec(Nloc,vect,vout)
       if(mpi_master)write(*,*)"|H*v-E*v|=",sum(abs(vout-egs*vect))
    endif
    Nitermax=Nlanc
    !
  end subroutine mpi_lanczos_eigh_d









  !---------------------------------------------------------------------
  !Purpose: plain homebrew lanczos iteration (no orthogonalization)
  !note: the a,b variables are real, even in the complex matrix case
  !to understand why check out the Gollub-Van Loan textbook.
  !a it is easy: hermiticity->diag\in\RRR
  !b: is fixed by requiring |b|^2 = <v,v> thus you can only fix the 
  !the absolute value. A lemma shows that the phase can be chosen 
  !identically zero
  !MPI VERSION
  !---------------------------------------------------------------------
  subroutine mpi_lanczos_iteration_d(MpiComm,MatVec,iter,vin,vout,alfa,beta)
    integer                                    :: MpiComm
    interface
       subroutine MatVec(Nloc,vin,vout)
         integer                 :: Nloc
         real(8),dimension(Nloc) :: vin
         real(8),dimension(Nloc) :: vout
       end subroutine MatVec
    end interface
    real(8),dimension(:),intent(inout)         :: vin !Nloc
    real(8),dimension(size(vin)),intent(inout) :: vout
    real(8),dimension(size(vin))               :: tmp
    real(8),intent(inout)                      :: alfa,beta
    real(8)                                    :: atmp,btmp
    integer                                    :: iter,nloc
    real(8)                                    :: norm,norm_tmp
    !
    logical                                    :: mpi_master
    !
    nloc=size(vin)
    !
    mpi_master=get_master_MPI(MpiComm)
    !
    if(iter==1)then
       norm_tmp=dot_product(vin,vin)
       norm = 0d0 ; call AllReduce_MPI(MpiComm,norm_tmp,norm)
       if(mpi_master.AND.norm==0d0)stop "MPI_LANCZOS_ITERATION_D: norm = 0!!"
       vin=vin/sqrt(norm)
    else
       tmp = vin
       vin = vout/beta
       vout= -beta*tmp
    endif
    call MatVec(nloc,vin,tmp)
    vout = vout + tmp
    atmp = dot_product(vin,vout) ; alfa = 0d0; call AllReduce_MPI(MpiComm,atmp,alfa)
    vout = vout - alfa*vin
    btmp = dot_product(vout,vout); beta = 0d0; call AllReduce_MPI(MpiComm,btmp,beta)
    beta = sqrt(beta)
  end subroutine mpi_lanczos_iteration_d



end module SF_SP_LINALG_
