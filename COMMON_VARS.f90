MODULE COMMON_VARS
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none



  !This is the local Mpi Communicator, shared across all this module.
  !The user MUST use the set/del routines to set/del the communicator
  !accordingly in order to use parallel lanczos procedure.
  !Although less elegant this covers quite all the cases without
  !introducing a distinction with respect to the serial MatVec procedures.
#ifdef _MPI
  integer :: MpiComm=MPI_UNDEFINED  
#endif
  logical :: MpiStatus=.false.
  integer :: MpiRank=0
  integer :: MpiSize=1
  logical :: MpiMaster=.true.

  integer :: mpiQup,mpiRup
  integer :: mpiQdw,mpiRdw


  integer :: DimUp
  integer :: DimDw
  integer :: Ns       !Number of levels per spin

  real(8) :: t_start,t_end
  
contains


  subroutine set_MpiComm(comm_)
#ifdef _MPI
    integer :: comm_
    MpiComm   = comm_
    MpiRank   = get_rank_MPI(MpiComm)
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    MpiStatus = .true.
#else
    integer,optional :: comm_
#endif
  end subroutine set_MpiComm


  subroutine del_MpiComm()
#ifdef _MPI
    MpiComm = MPI_UNDEFINED
    MpiRank   = 0
    MpiSize   = 1
    MpiMaster = .false.
    MpiStatus = .true.
#endif
  end subroutine del_MpiComm





end module COMMON_VARS

