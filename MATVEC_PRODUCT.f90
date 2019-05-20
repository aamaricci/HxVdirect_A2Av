MODULE MATVEC_PRODUCT
  USE COMMON_VARS
  USE SETUP
  USE SF_LINALG
  USE SF_RANDOM
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none

  private

  public :: Setup_HxV
  public :: Build_HxV

  public :: HxVdirect
#ifdef _MPI
  public :: mpi_HxVdirect
#endif




  integer                       :: Dim

  !Hamiltonian hard coded paramters:
  !Hamiltonian hard coded paramters:
  real(8)                       :: e0,xmu
  real(8)                       :: Uloc
  !
  type(sector_map),dimension(2) :: Hs
  !
  integer                       :: mpi_rank
  integer                       :: mpi_size
  integer                       :: mpi_Q
  integer                       :: mpi_R
  integer                       :: mpi_Istart
  integer                       :: mpi_Iend
  integer                       :: mpi_Ishift

contains




  !####################################################################
  !             SETUP THE MPI ENVIRONMENT
  !####################################################################
  subroutine Setup_HxV(Nups,Ndws)
    integer :: Nups,Ndws,ie
    real(8) :: Deps
    integer :: seed
    !
    call build_sector(Nups,Ndws,Hs)
    !
    DimUp = get_sector_dimension(Ns,Nups)
    DimDw = get_sector_dimension(Ns,Ndws)
    Dim=DimUp*DimDw
    !
    !
    !The MPI version can allocate directly from the total dimension,
    !evaluating the chunks independently.
#ifdef _MPI
    mpi_rank = get_rank_MPI(MpiCOmm)
    mpi_size = get_size_MPI(MpiComm)
    !
    !Dw split:
    mpiQdw = DimDw/mpi_size
    mpiRdw = mod(DimDw,mpi_size)
    if(mpi_rank < mod(DimDw,mpi_size) ) then
       mpiRdw = 0
       MpiQdw = MpiQdw+1
    endif
    !Wait here to parse/read Dim
    mpi_Q = DimUp*mpiQdw
    mpi_R = DimUp*mpiRdw
    mpi_Istart = 1 + mpi_rank*mpi_Q+mpi_R
    mpi_Iend   = (mpi_rank+1)*mpi_Q+mpi_R
    mpi_Ishift = mpi_Rank*mpi_Q+mpi_R
    !
#else
    mpi_Istart = 1
    mpi_Iend   = Dim
#endif
  end subroutine Setup_HxV



  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine Build_HxV()
    !
    ! e0  = 0.2d0
    ! xmu = 0.5d0
    ! ts  = 1d0/sqrt(2d0)
    ! Uloc = 2d0
    !
    !< generate Hamiltonian parameters:
    e0  = -0.3d0
    !
    Uloc = 2d0
    !
  end subroutine Build_HxV







  subroutine HxVdirect(N,v,Hv)
    integer               :: N
    real(8),dimension(N)  :: v
    real(8),dimension(N)  :: Hv
    integer,dimension(Ns) :: nup,ndw
    integer               :: i,iup,idw
    integer               :: j,jup,jdw
    integer               :: m,mup,mdw
    integer               :: ms
    integer               :: impi,impi_up,impi_dw
    integer               :: iorb,jorb,ispin,jspin,ibath
    integer               :: kp,k1,k2,k3,k4
    integer               :: alfa,beta
    real(8)               :: sg1,sg2,sg3,sg4
    real(8)               :: htmp
    !
    Hv=0d0
    !
    !
    !LOCAL PART
    do i=mpi_Istart,mpi_Iend
       iup = iup_index(i,DimUp)
       idw = idw_index(i,DimUp)
       !
       nup  = bdecomp(Hs(1)%map(iup),Ns)
       ndw  = bdecomp(Hs(2)%map(idw),Ns)
       !
       !Diagonal Elements, i.e. local part
       htmp = 0d0
       htmp = htmp + e0*(nup(1)+ndw(1))
       htmp = htmp + Uloc*(nup(1)*ndw(1))
       htmp = htmp - 0.5d0*Uloc*(nup(1)+ndw(1)) + 0.25d0*uloc
       !
       !<
       !H_Bath: local bath energy contribution. Assume this is zero.
       ! do kp=1,Ns-1
       !    alfa = 1 + kp
       !    htmp =htmp + eps(kp)*nup(alfa) !UP
       !    htmp =htmp + eps(kp)*ndw(alfa) !DW
       ! enddo
       !
       Hv(i) = Hv(i) + htmp*v(i)
    enddo
    !
    ! IMP UP <--> BATH UP
    ! do iup=1,DimUp!first_state_up,last_state_up
    do idw=1,DimDw
       !
       do iup=1,DimUp
          mup  = Hs(1)%map(iup)
          nup  = bdecomp(mup,Ns)
          !
          i = iup + (idw-1)*DimUp
          !
          do kp=1,Ns-1
             alfa = 1 + kp
             if( (nup(1)==1) .AND. (nup(alfa)==0) )then              
                call c(1,mup,k1,sg1)
                call cdg(alfa,k1,k2,sg2)
                jup = binary_search(Hs(1)%map,k2)
                htmp = vps*sg1*sg2
                j = jup + (idw-1)*DimUp
                Hv(i) = Hv(i) + htmp*V(j)
             endif
             if( (nup(1)==0) .AND. (nup(alfa)==1) )then
                call c(alfa,mup,k1,sg1)
                call cdg(1,k1,k2,sg2)
                jup=binary_search(Hs(1)%map,k2)
                htmp = vps*sg1*sg2
                j = jup + (idw-1)*DimUp
                Hv(i) = Hv(i) + htmp*V(j)
             endif
          enddo
       enddo
       !
    enddo
    !
    !IMP DW <--> BATH DW
    do iup=1,DimUp
       do idw=1,DimDw
          mdw  = Hs(2)%map(idw)
          ndw  = bdecomp(mdw,Ns)
          !
          i = iup + (idw-1)*DimUp
          !
          do kp=1,Ns-1
             alfa = 1 + kp
             if( (ndw(1)==1) .AND. (ndw(alfa)==0) )then
                call c(1,mdw,k1,sg1)
                call cdg(alfa,k1,k2,sg2)
                jdw=binary_search(Hs(2)%map,k2)
                htmp=vps*sg1*sg2
                j = iup +  (jdw-1)*DimUp
                Hv(i) = Hv(i) + htmp*V(j)
             endif
             !
             !
             if( (ndw(1)==0) .AND. (ndw(alfa)==1) )then
                call c(alfa,mdw,k1,sg1)
                call cdg(1,k1,k2,sg2)
                jdw=binary_search(Hs(2)%map,k2)
                htmp=vps*sg1*sg2
                j = iup +  (jdw-1)*DimUp
                Hv(i) = Hv(i) + htmp*V(j)
             endif
          enddo
       enddo
    enddo
    !
    return
  end subroutine HxVdirect






  !THESE ARE THE *MPI* MATVEC PRODUCTS USING "small" VECTORS
#ifdef _MPI
  subroutine mpi_HxVdirect(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v,Hv
    !
    integer                          :: N
    real(8),dimension(:),allocatable :: vt,Hvt
    !local MPI
    integer                          :: irank
    integer                          :: MpiQup,MpiQdw
    integer,dimension(Ns)            :: nup,ndw
    integer                          :: i,iup,idw
    integer                          :: j,jup,jdw
    integer                          :: m,mup,mdw
    integer                          :: ms
    integer                          :: impi,impi_up,impi_dw
    integer                          :: iorb,jorb,ispin,jspin,ibath
    integer                          :: kp,k1,k2,k3,k4
    integer                          :: alfa,beta
    real(8)                          :: sg1,sg2,sg3,sg4
    real(8)                          :: htmp
    !
    if(MpiComm==MPI_UNDEFINED)stop "mpi_HxV_d ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "mpi_HxV_d ERRROR: MpiStatus = F"   
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=0d0
    !
    !LOCAL PART
    do i=1,Nloc
       iup = iup_index(i+mpi_Ishift,DimUp)
       idw = idw_index(i+mpi_Ishift,DimUp)
       !
       nup  = bdecomp(Hs(1)%map(iup),Ns)
       ndw  = bdecomp(Hs(2)%map(idw),Ns)
       !
       !Diagonal Elements, i.e. local part
       htmp = 0d0
       htmp = htmp + e0*(nup(1)+ndw(1))
       htmp = htmp + Uloc*(nup(1)*ndw(1))
       htmp = htmp - 0.5d0*Uloc*(nup(1)+ndw(1)) + 0.25d0*uloc
       !
       !> H_Bath: local bath energy contribution.
       !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
       ! do kp=1,Ns-1
       !    alfa = 1 + kp
       !    htmp =htmp + eps(kp)*nup(alfa) !UP
       !    htmp =htmp + eps(kp)*ndw(alfa) !DW
       ! enddo
       !
       Hv(i) = Hv(i) + htmp*v(i)
    enddo
    !
    !
    mpiQdw=DimDw/MpiSize
    if(MpiRank<mod(DimDw,MpiSize))MpiQdw=MpiQdw+1
    !
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    !ADD SOME TESTS HERE:
    !
    !Non-local terms.
    !UP part: contiguous in memory.
    ! IMP UP <--> BATH UP
    do idw=1,MpiQdw
       do iup=1,DimUp
          mup  = Hs(1)%map(iup)
          nup  = bdecomp(mup,Ns)
          !
          i = iup + (idw-1)*DimUp
          !
          do kp=1,Ns-1
             alfa = 1 + kp
             if( (nup(1)==1) .AND. (nup(alfa)==0) )then              
                call c(1,mup,k1,sg1)
                call cdg(alfa,k1,k2,sg2)
                jup = binary_search(Hs(1)%map,k2)
                htmp = vps*sg1*sg2
                j = jup + (idw-1)*DimUp
                Hv(i) = Hv(i) + htmp*V(j)
             endif
             if( (nup(1)==0) .AND. (nup(alfa)==1) )then
                call c(alfa,mup,k1,sg1)
                call cdg(1,k1,k2,sg2)
                jup=binary_search(Hs(1)%map,k2)
                htmp = vps*sg1*sg2
                j = jup + (idw-1)*DimUp
                Hv(i) = Hv(i) + htmp*V(j)
             endif
          enddo
       enddo
    enddo
    !






    ! !
    !DW part: non-contiguous in memory -> MPI transposition
    !Transpose the input vector as a whole:
    allocate(vt(mpiQup*DimDw)) ;vt=0d0
    allocate(Hvt(mpiQup*DimDw));Hvt=0d0    
    call vector_transpose_MPI(DimUp,MpiQdw,v,DimDw,MpiQup,vt)
    !
    Hvt=0d0
    do iup=1,MpiQup
       do idw=1,DimDw
          mdw  = Hs(2)%map(idw)
          ndw  = bdecomp(mdw,Ns)
          !
          i = idw + (iup-1)*DimDw
          !
          do kp=1,Ns-1
             alfa = 1 + kp
             if( (ndw(1)==1) .AND. (ndw(alfa)==0) )then
                call c(1,mdw,k1,sg1)
                call cdg(alfa,k1,k2,sg2)
                jdw=binary_search(Hs(2)%map,k2)
                htmp=vps*sg1*sg2
                j = jdw + (iup-1)*DimDw
                Hvt(i) = Hvt(i) + htmp*vt(j)
             endif
             !
             !
             if( (ndw(1)==0) .AND. (ndw(alfa)==1) )then
                call c(alfa,mdw,k1,sg1)
                call cdg(1,k1,k2,sg2)
                jdw=binary_search(Hs(2)%map,k2)
                htmp=vps*sg1*sg2
                j = jdw + (iup-1)*DimDw
                Hvt(i) = Hvt(i) + htmp*vt(j)
             endif
          enddo
       enddo
    enddo
    !
    deallocate(vt)
    allocate(vt(DimUp*mpiQdw));vt=0d0
    call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt)
    Hv = Hv + vt
    !
    return
  end subroutine Mpi_HxVdirect


  subroutine vector_transpose_MPI(nrow,qcol,a,ncol,qrow,b)    
    integer                            :: nrow,ncol,qrow,qcol
    real(8)                            :: a(nrow,qcol)
    real(8)                            :: b(ncol,qrow)
    integer,allocatable,dimension(:,:) :: send_counts,send_offset
    integer,allocatable,dimension(:,:) :: recv_counts,recv_offset
    integer                            :: counts,Ntot
    integer :: i,j,irank,ierr
    !
    counts = Nrow/MpiSize
    Ntot   = Ncol/MpiSize
    if(mod(Ncol,MpiSize)/=0)Ntot=Ntot+1
    !
    allocate(send_counts(0:MpiSize-1,Ntot));send_counts=0
    allocate(send_offset(0:MpiSize-1,Ntot));send_offset=0
    allocate(recv_counts(0:MpiSize-1,Ntot));recv_counts=0
    allocate(recv_offset(0:MpiSize-1,Ntot));recv_offset=0
    !
    do i=1,qcol
       do irank=0,MpiSize-1
          if(irank < mod(Nrow,MpiSize))then
             send_counts(irank,i) = counts+1
          else
             send_counts(irank,i) = counts
          endif
       enddo
    enddo
    !
    do i=1,Ntot
       call MPI_AllToAll(&
            send_counts(:,i),1,MPI_INTEGER,&
            recv_counts(:,i),1,MPI_INTEGER,&
            MpiComm,ierr)
    enddo
    !
    do i=1,Ntot
       do irank=1,MpiSize-1
          send_offset(irank,i) = send_counts(irank-1,i) + send_offset(irank-1,i)
       enddo
    enddo
    !
    !Get the irank=0 elements, i.e. first entries:
    recv_offset(0,1) = 0
    do i=2,Ntot
       recv_offset(0,i) = sum(recv_counts(0,:i-1))
    enddo
    !the rest of the entries:
    do i=1,Ntot
       do irank=1,MpiSize-1
          recv_offset(irank,i) = recv_offset(irank-1,i) + sum(recv_counts(irank-1,:))
       enddo
    enddo
    !
    !
    do j=1,Ntot
       call MPI_AllToAllV(&
            A(:,j),send_counts(:,j),send_offset(:,j),MPI_DOUBLE_PRECISION,&
            B(:,:),recv_counts(:,j),recv_offset(:,j),MPI_DOUBLE_PRECISION,&
            MpiComm,ierr)
    enddo
    !
    call local_transpose(b,ncol,qrow)
    !
    return
  end subroutine vector_transpose_MPI


  subroutine local_transpose(mat,nrow,ncol)
    integer                      :: nrow,ncol
    real(8),dimension(Nrow,Ncol) :: mat
    mat = transpose(reshape(mat,[Ncol,Nrow]))
  end subroutine local_transpose
#endif

end module MATVEC_PRODUCT







