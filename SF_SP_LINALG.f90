module SF_SP_LINALG_
  USE SF_LINALG, only: eye,eigh
  USE SF_RANDOM, only: mt_random,mersenne_init
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private 

  interface sp_eigh
     module procedure :: lanczos_arpack_d
#ifdef _MPI
     module procedure :: lanczos_parpack_d
#endif
  end interface sp_eigh


  interface sp_lanc_eigh
     module procedure :: lanczos_eigh_d
#ifdef _MPI
     module procedure :: mpi_lanczos_eigh_d
#endif
  end interface sp_lanc_eigh




  complex(8),parameter              :: zero=(0d0,0d0)
  complex(8),parameter              :: one=(1d0,0d0)
  complex(8),parameter              :: xi=(0d0,1d0)
  integer,allocatable               :: seed_random(:)
  integer                           :: nrandom
  logical                           :: verb=.false.
  real(8)                           :: threshold_=1.d-12
  integer                           :: ncheck_=1



  !****************************************************************************************
  !                                      PUBLIC 
  !****************************************************************************************
  public :: sp_eigh
  public :: sp_lanc_eigh
  !****************************************************************************************




contains



  !##################################################################
  ! ARPACK METHOD for LOWEST part of the spectrum of a of a SPARSE
  !   MATRIX (defined via H*v)
  ! - DBLE and CMPLX versions included
  ! - SERIAL and PARALLEL-MPI versions included
  !##################################################################
  include "arpack_serial.f90"
#ifdef _MPI
  include "arpack_mpi.f90"
#endif


  !##################################################################
  ! LANCZOS METHOD for LOWEST EigenSolution OR tri-diagonalization
  !    of a SPARSE MATRIX (defined via H*v)
  ! - DBLE and CMPLX versions included
  ! - SERIAL and PARALLEL-MPI versions included
  !##################################################################
  include "lanczos_serial.f90"
#ifdef _MPI
  include "lanczos_mpi.f90"
#endif






end module SF_SP_LINALG_
