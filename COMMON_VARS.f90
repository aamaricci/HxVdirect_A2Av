MODULE COMMON_VARS
  USE SF_MPI
  USE MPI
  implicit none



  !This is the local Mpi Communicator, shared across all this module.
  !The user MUST use the set/del routines to set/del the communicator
  !accordingly in order to use parallel lanczos procedure.
  !Although less elegant this covers quite all the cases without
  !introducing a distinction with respect to the serial MatVec procedures.
  integer :: MpiComm=MPI_UNDEFINED  
  logical :: MpiStatus=.false.
  integer :: MpiRank=0
  integer :: MpiSize=1
  logical :: MpiMaster=.true.

  integer :: mpiQup,mpiRup
  integer :: mpiQdw,mpiRdw


  integer :: DimUp
  integer :: DimDw
  integer :: Ns       !Number of levels per spin
  integer :: Norb     !Number of orbitals
  integer :: Nbath    !Number of bath levels per orbital and spin
  !
  real(8) :: t_start,t_end

  !Hamiltonian hard coded paramters:
  real(8),allocatable :: Vps(:)
  real(8),allocatable :: Eps(:)
  real(8),allocatable :: E0(:)
  real(8)             :: Mh
  real(8)             :: Uloc,Ust,Jh
  logical             :: jhflag=.false.
  integer,allocatable :: get_bath_index(:,:)
  
contains


  subroutine set_EDparameters(Na,Nb)
    integer :: Na,Nb
    integer :: iorb,ibath
    Norb  = Na
    Nbath = Nb
    Ns    = (Nbath+1)*Norb
    !
    if(allocated(Eps))deallocate(Eps)
    if(allocated(Vps))deallocate(Vps)
    if(allocated(E0))deallocate(E0)
    if(allocated(get_bath_index))deallocate(get_bath_index)
    !
    allocate(Eps(Nbath))
    allocate(Vps(Nbath))
    allocate(E0(Norb))
    allocate(get_bath_index(Norb,Nbath))
    do ibath=1,Nbath
       do iorb=1,Norb
          get_bath_index(iorb,ibath) = Norb + (iorb-1)*Nbath + ibath
       enddo
    enddo
  end subroutine set_EDparameters


  subroutine set_MpiComm(comm_)
    integer :: comm_
    MpiComm   = comm_
    MpiRank   = get_rank_MPI(MpiComm)
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    MpiStatus = .true.
  end subroutine set_MpiComm


  subroutine del_MpiComm()
    MpiComm = MPI_UNDEFINED
    MpiRank   = 0
    MpiSize   = 1
    MpiMaster = .false.
    MpiStatus = .true.
  end subroutine del_MpiComm





end module COMMON_VARS

