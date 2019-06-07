program test_MPI_LANCZOS_D
  USE COMMON_VARS
  USE SETUP
  USE MATVEC_PRODUCT
  !
  USE SF_PARSE_INPUT
  USE SF_LINALG
  USE SF_TIMER
  !
  USE SF_SP_LINALG_
  !
  USE SF_MPI
  USE MPI
  implicit none
  !Dimension
  integer             :: nup,ndw,Dim
  integer             :: Ntot
  !Matrix:
  real(8),allocatable :: Hmat(:,:)
  integer             :: N,Nloc
  real(8),allocatable :: Eval(:),Eval_old(:)
  real(8),allocatable :: Evec(:,:)
  !Lanczos
  integer             :: neigen,nitermax,niter,nblock
  integer             :: i,j,k,col,irank
  integer             :: iup,idw,jup,jdw,ncv
  !MPI:
  integer             :: comm
  integer             :: mpi_rank
  integer             :: mpi_size
  integer             :: mpi_ierr
  logical             :: mpi_master,bool
  integer             :: mpi_Q,Q,mpi_Istart
  integer             :: mpi_R,R,mpi_Iend
  integer             :: mpi_Chunk,mpi_Ishift
  !
  character(len=100)  :: finput

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)

  mpi_rank   = get_rank_MPI(comm)
  mpi_size   = get_size_MPI(comm)
  mpi_master = get_master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT","inputHxV.conf")
  call parse_input_variable(norb,"NORB",finput,default=1)
  call parse_input_variable(nbath,"NBATH",finput,default=9)
  call parse_input_variable(jhflag,"JHFLAG",finput,default=.false.)
  call parse_input_variable(Neigen,"NEIGEN",finput,default=1)
  call parse_input_variable(NCV,"NCV",finput,default=10)
  call parse_input_variable(Nvps,"NVPS",finput,default=1)
  if(mpi_master)then
     call print_input()
     call save_input(trim(finput))
  endif
  !
  !
  call set_EDparameters(Norb,Nbath)
  Nup = Ns/2
  Ndw = Ns-Nup
  !
  DimUp  = binomial(Ns,Nup)
  DimDw  = binomial(Ns,Ndw)
  !
  !
  if(mpi_master)then
     print*,"Ns         =",Ns
     print*,"Norb       =",Norb
     print*,"Nbath      =",Nbath
     print*,"Nup,Ndw    =",Nup,Ndw
     print*,"Dimensions =",DimUp,DimDw
  endif
  !
  call set_MpiComm(comm)
  !
  !Vectors have dimension: DimUp * mpiQdw (each threads takes mpiQdw columns of full len DimUp)
  call Setup_HxV(Nup,Ndw)
  !
  !
  if(mpi_master)print*,"ARPACK DIAG D - MPI:"
  !
  allocate(Eval(Neigen),Evec(DimUp*mpiQdw,Neigen))
  Eval=0d0
  Evec=0d0
  !
  Nitermax=512
  Nblock=ncv*Neigen
  !
  if(mpi_master)call start_timer
  t_start = MPI_Wtime()
  call sp_eigh(comm,mpi_HxVdirect,Neigen,Nblock,Nitermax,Eval,Evec,iverbose=.true.)
  if(mpi_master)call stop_timer
  t_end = MPI_Wtime()
  !
  if(mpi_master)then
     do i=1,Neigen
        write(*,*)Eval(i)
     end do
     open(100,file="mpi_arpack_Eval_D.dat",position="append")
     do i=1,Neigen
        write(100,*)Eval(i)
     end do
     close(100)
     open(100,file="mpi_arpack_time.dat",position="append")
     write(100,*)t_end-t_start
     close(100)
     print*,""
  endif
  deallocate(Eval,Evec)
  call Finalize_MPI()


end program test_MPI_LANCZOS_D








