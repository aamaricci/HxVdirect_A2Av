MODULE MATVEC_PRODUCT
  USE COMMON_VARS
  USE SETUP
  USE SCIFOR
  USE MPI
  implicit none

  private

  public :: Setup_HxV
  public :: mpi_HxVdirect




  integer                       :: Dim
  !
  type(sector_map),dimension(2) :: Hs
  !
  integer                       :: mpi_rank
  integer                       :: mpi_size
  integer                       :: mpi_Q
  integer                       :: mpi_R
  integer                       :: mpi_Istart
  integer                       :: mpi_Iend
  integer                       :: mpi_Ishift,mpiIdw


contains




  !####################################################################
  !             SETUP THE MPI ENVIRONMENT
  !####################################################################
  subroutine Setup_HxV(Nups,Ndws)
    integer :: Nups,Ndws,iorb,irank
    real(8) :: Deps
    integer :: seed,i
    integer :: iup,idw
    integer :: impi_up,impi_dw
    !
    !
    ! eps = 0d0
    ! eps = linspace(-2d0-mersenne(),2d0+mersenne(),Nbath)
    !
    ! call mt_random(vps)
    ! vps = vps/sqrt(5d0)
    vps = 1d0/sqrt(2d0)
    !
    Mh  = 1d0
    !
    !Set some Hamiltonian parameters:
    select case(Norb)
    case default
       stop "Norb not in (1,5)"
    case (1)
       E0(1) = 0d0
    case (2)
       E0(1) =  Mh
       E0(2) = -Mh
    case (3)
       E0(1) =  Mh
       E0(2) =  0d0
       E0(3) = -Mh
    case (4)
       E0(1) =  2*Mh
       E0(2) =  Mh
       E0(3) = -Mh
       E0(4) = -2*Mh
    case (5)
       E0(1) =  2*Mh
       E0(2) =  Mh
       E0(3) =  0d0
       E0(4) = -Mh
       E0(5) = -2*Mh
    end select
    !
    Uloc = 1d0
    Jh   = Uloc/5d0
    Ust  = Uloc - 2d0*Jh

    !
    if(MpiMaster)then
       ! print*,"Eps:",eps
       print*,"Vps:",vps
       print*,"E0 :",E0
       print*,"U  :",Uloc
       print*,"U` :",Ust
       print*,"U``:",Ust-Jh
       print*,"Jh :",Jh
       print*,"Jx :",Jh
       print*,"Jp :",Jh
    endif
    !
    call build_sector(Nups,Ndws,Hs)
    !
    DimUp = get_sector_dimension(Ns,Nups)
    DimDw = get_sector_dimension(Ns,Ndws)
    !
    !
    !The MPI version can allocate directly from the total dimension,
    !evaluating the chunks independently.
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
    mpiIdw = mpi_Rank*MpiQdw + mpiRdw
    !Wait here to parse/read Dim
    mpi_Q = DimUp*mpiQdw
    mpi_R = DimUp*mpiRdw
    mpi_Istart = 1 + mpi_rank*mpi_Q+mpi_R
    mpi_Iend   = (mpi_rank+1)*mpi_Q+mpi_R
    mpi_Ishift = mpi_Rank*mpi_Q+mpi_R

    call Barrier_MPI(MpiComm)
    if(MpiMaster)write(*,"(9A20)")"mpiRank","mpi_Q","mpi_R",&
         "mpi_Qdw","mpiR_dw","mpi_Ishift","mpi_Istart","mpi_Iend","mpi_Iend-mpi_Istart"
    call Barrier_MPI(MpiComm)
    do irank=0,mpi_size-1
       call Barrier_MPI(MpiComm)
       if(mpi_rank==irank)write(*,"(9I20)")Mpi_Rank,Mpi_Q,Mpi_R,mpiQdw,MpiRdw,Mpi_Ishift,Mpi_Istart,Mpi_Iend,Mpi_Iend-Mpi_Istart+1
    enddo
    call Barrier_MPI(MpiComm)
    !
    ! do i=1,MpiQdw*DimUp
    !    write(700+Mpi_Rank,*)i
    ! enddo
    !
    ! do idw=1,MpiQdw
    !    impi_dw = idw + mpiIdw
    !    do iup=1,DimUp
    !       write(800+Mpi_Rank,*)iup+(impi_dw-1)*DimUp
    !       write(300+Mpi_Rank,*)iup + (idw-1)*DimUp
    !    enddo
    ! enddo
    !
  end subroutine Setup_HxV




  !THESE ARE THE *MPI* MATVEC PRODUCTS USING "small" VECTORS
  subroutine mpi_HxVdirect(Nloc,v,Hv)
    integer                          :: Nloc
    real(8),dimension(Nloc)          :: v,Hv
    !
    integer                          :: N
    real(8),dimension(:),allocatable :: vt,Hvt
    !local MPI
    integer                          :: MpiQup_,MpiQdw_
    !
    integer,dimension(Ns)            :: nup,ndw
    integer                          :: i,iup,idw
    integer                          :: j,jup,jdw
    integer                          :: m,mup,mdw
    integer                          :: impi,impi_up,impi_dw
    integer                          :: iorb,jorb,ibath
    integer                          :: kp,k1,k2,k3,k4
    integer                          :: alfa
    real(8)                          :: sg1,sg2,sg3,sg4
    real(8)                          :: htmp
    logical                          :: Jcondition
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    Hv=0d0
    !LOCAL PART: do a suitable double loop here to avoid having too large integers
    !            as in sectors which are as large as 2**31
    do impi_dw = 1,MpiQdw
       idw = impi_dw + mpiIdw 
       do iup = 1,DimUp
          i = iup + (impi_dw-1)*DimUp
          !
          nup  = bdecomp(Hs(1)%map(iup),Ns)
          ndw  = bdecomp(Hs(2)%map(idw),Ns)
          !
          !Diagonal Elements, i.e. local part
          htmp = 0d0
          do iorb=1,Norb
             htmp = htmp + e0(iorb)*(nup(iorb)+ndw(iorb))
          enddo
          !
          !> H_Int: Kanamori interaction part. non-local S-E and P-H terms commented below.
          !
          ! density-density interaction: \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
          !                              + contributions of hartree terms
          do iorb=1,Norb
             htmp = htmp + Uloc*Nup(iorb)*Ndw(iorb)
             htmp = htmp - 0.5d0*Uloc*(Nup(iorb)+Ndw(iorb)) + Norb*0.25d0*uloc
          enddo
          !
          !
          if(Norb>1)then
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   !density-density interaction: different orbitals, opposite spins:
                   ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
                   ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
                   htmp = htmp + Ust*(Nup(iorb)*Ndw(jorb) + Nup(jorb)*Ndw(iorb))
                   !
                   !density-density interaction: different orbitals, parallel spins
                   ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
                   ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
                   htmp = htmp + (Ust-Jh)*(Nup(iorb)*Nup(jorb) + Ndw(iorb)*Ndw(jorb))
                   !
                   !using the Hartree-shifted chemical potential: mu=0 for half-filling
                   htmp=htmp-0.5d0*Ust*(Nup(iorb)+Ndw(iorb)+Nup(jorb)+Ndw(jorb))+0.25d0*Ust
                   htmp=htmp-0.5d0*(Ust-Jh)*(Nup(iorb)+Ndw(iorb)+Nup(jorb)+Ndw(jorb))+0.25d0*(Ust-Jh)
                enddo
             enddo
          endif
          !
          !> H_Bath: local bath energy contribution. Assume zero
          ! do iorb=1,Norb
          !    do kp=1,Nbath
          !       alfa = get_bath_index(iorb,kp)
          !       htmp =htmp + eps(kp)*(Nup(alfa)+Ndw(alfa))
          !    enddo
          ! enddo
          Hv(i) = Hv(i) + htmp*v(i)
          !
       enddo
    enddo
    !
    !
    !
    !< Local splitting of the Dw and Up sectors
    mpiQdw_=DimDw/MpiSize
    if(MpiRank<mod(DimDw,MpiSize))MpiQdw_=MpiQdw_+1
    !
    mpiQup_=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup_=MpiQup_+1
    !
    !
    !Non-local terms.
    !UP PART: CONTIGUOUS IN MEMORY.
    do idw=1,MpiQdw_               !for any DW state
       !
       do iup=1,DimUp             !loop over UP states
          mup  = Hs(1)%map(iup)   !map up_sector_state 2 up_fock_state 
          nup  = bdecomp(mup,Ns)  !get binary decomposition=occupation
          !
          i = iup + (idw-1)*DimUp !get overall fock state
          !
          do iorb=1,Norb
             do kp=1,Nbath        !loop over bath levels
                !
                !< get index of the bath sites in the fock ket 
                alfa = get_bath_index(iorb,kp)
                !
                if( (nup(iorb)==1) .AND. (nup(alfa)==0) )then              
                   call c(iorb,mup,k1,sg1)
                   call cdg(alfa,k1,k2,sg2)
                   jup = binary_search(Hs(1)%map,k2)
                   htmp = vps(kp)*sg1*sg2
                   j = jup + (idw-1)*DimUp
                   Hv(i) = Hv(i) + htmp*V(j)
                   Hv(j) = Hv(j) + htmp*V(i)
                endif
             enddo
          enddo
          !
          !
       enddo
    enddo
    !
    !
    !
    !DW PART: NON-CONTIGUOUS IN MEMORY -> MPI TRANSPOSITION
    !Transpose the input vector as a whole:
    allocate(vt(mpiQup_*DimDw)) ;vt=0d0
    allocate(Hvt(mpiQup_*DimDw));Hvt=0d0    
    call vector_transpose_MPI(DimUp,MpiQdw_,v,DimDw,MpiQup_,vt)
    do iup=1,MpiQup_
       do idw=1,DimDw
          mdw  = Hs(2)%map(idw)
          ndw  = bdecomp(mdw,Ns)
          i = idw + (iup-1)*DimDw
          !
          do iorb=1,Norb
             do kp=1,Nbath
                !
                alfa = get_bath_index(iorb,kp)
                !
                if( (ndw(iorb)==1) .AND. (ndw(alfa)==0) )then
                   call c(iorb,mdw,k1,sg1)
                   call cdg(alfa,k1,k2,sg2)
                   jdw=binary_search(Hs(2)%map,k2)
                   htmp=vps(kp)*sg1*sg2
                   j = jdw + (iup-1)*DimDw
                   Hvt(i) = Hvt(i) + htmp*vt(j)
                   Hvt(j) = Hvt(j) + htmp*vt(i)
                endif
             enddo
          enddo
          !
          !
       enddo
    enddo
    deallocate(vt) ; allocate(vt(DimUp*mpiQdw_)) ; vt=0d0
    call vector_transpose_MPI(DimDw,mpiQup_,Hvt,DimUp,mpiQdw_,vt)
    Hv = Hv + vt
    deallocate(vt)
    !
    !
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       !
       allocate(vt(N)) ; vt = 0d0
       call allgather_vector_MPI(MpiComm,v,vt)
       !
       ! do i=1,Nloc
       !    iup = iup_index(i+mpi_Ishift,DimUp)
       !    idw = idw_index(i+mpi_Ishift,DimUp)
       do impi_dw = 1,MpiQdw
          idw = impi_dw + mpiIdw 
          do iup = 1,DimUp
             i = iup + (impi_dw-1)*DimUp
             !
             mup = Hs(1)%map(iup)
             mdw = Hs(2)%map(idw)
             !
             nup  = bdecomp(mup,Ns)
             ndw  = bdecomp(mdw,Ns)
             !
             ! SPIN-EXCHANGE (S-E) 
             !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
             !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
             do iorb=1,Norb
                do jorb=1,Norb
                   Jcondition=(&
                        (iorb/=jorb).AND.&
                        (nup(jorb)==1).AND.&
                        (ndw(iorb)==1).AND.&
                        (ndw(jorb)==0).AND.&
                        (nup(iorb)==0))
                   if(Jcondition)then
                      call c(jorb,mup,k1,sg1)  !UP
                      call cdg(iorb,k1,k2,sg2) !UP
                      jup=binary_search(Hs(1)%map,k2)
                      call c(iorb,mdw,k3,sg3)  !DW
                      call cdg(jorb,k3,k4,sg4) !DW
                      jdw=binary_search(Hs(2)%map,k4)
                      htmp = Jh*sg1*sg2*sg3*sg4
                      j = jup + (jdw-1)*dimup
                      !
                      Hv(i) = Hv(i) + htmp*vt(j)
                      !
                   endif
                enddo
             enddo
             !
             ! PAIR-HOPPING (P-H) TERMS
             !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
             !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
             do iorb=1,Norb
                do jorb=1,Norb
                   Jcondition=(&
                        (iorb/=jorb).AND.&
                        (nup(jorb)==1).AND.&
                        (ndw(jorb)==1).AND.&
                        (ndw(iorb)==0).AND.&
                        (nup(iorb)==0))
                   if(Jcondition)then
                      call c(jorb,mup,k1,sg1)       !c_jorb_up
                      call cdg(iorb,k1,k2,sg2)      !c^+_iorb_up
                      jup = binary_search(Hs(1)%map,k2)
                      call c(jorb,mdw,k3,sg3)       !c_jorb_dw
                      call cdg(iorb,k3,k4,sg4)      !c^+_iorb_dw
                      jdw = binary_search(Hs(2)%map,k4)
                      htmp = Jh*sg1*sg2*sg3*sg4
                      j = jup + (jdw-1)*dimup
                      !
                      Hv(i) = Hv(i) + htmp*vt(j)
                      !
                   endif
                enddo
             enddo
             !     
          enddo
       enddo
       !
       deallocate(Vt)
    endif
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



  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine allgather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
    integer                          :: i,irank,Nloc,N
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N    = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "allgather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    V = 0d0
    call MPI_AllGatherv(Vloc,Nloc,MPI_DOUBLE_PRECISION,V,Counts,Offset,MPI_DOUBLE_PRECISION,MpiComm,MpiIerr)
    !
    return
  end subroutine Allgather_vector_MPI




end module MATVEC_PRODUCT







