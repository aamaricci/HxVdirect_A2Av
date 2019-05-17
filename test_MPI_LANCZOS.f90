program test_MPI_LANCZOS_D
  USE COMMON_VARS
  USE SETUP
  USE MATVEC_PRODUCT
  !
  USE SF_PARSE_INPUT
  USE SF_TIMER
  !
  USE SF_SP_LINALG_
  !
  USE SF_MPI
  USE MPI
  implicit none
  !Dimension
  integer                          :: nup,ndw,Dim
  integer                          :: Ntot,Nsparse,Ncycles
  !Matrix:
  real(8),allocatable              :: Hmat(:,:)
  integer                          :: Nprint,N,Nloc
  real(8),allocatable              :: Eval(:),Eval_old(:)
  real(8),allocatable              :: Evec(:,:)
  !Lanczos
  integer                          :: neigen,nitermax,niter,nblock,Nlanc
  integer                          :: i,j,k,col,irank
  integer                          :: iup,idw,jup,jdw,ncv,icycle
  !MPI:
  integer                          :: comm
  integer                          :: mpi_rank
  integer                          :: mpi_size
  integer                          :: mpi_ierr
  logical                          :: mpi_master,bool
  integer                          :: mpi_Q,Q,mpi_Istart
  integer                          :: mpi_R,R,mpi_Iend
  integer                          :: mpi_Chunk,mpi_Ishift

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)

  mpi_rank   = get_rank_MPI(comm)
  mpi_size   = get_size_MPI(comm)
  mpi_master = get_master_MPI(comm)

  call parse_cmd_variable(ntot,"ntot",default=12)
  call parse_cmd_variable(nup,"nup",default=6)
  call parse_cmd_variable(ndw,"ndw",default=6)
  call parse_cmd_variable(Neigen,"NEIGEN",default=1)
  call parse_cmd_variable(NCV,"NCV",default=10)
  call parse_cmd_variable(Ncycles,"NCYCLES",default=1)
  !
  !
  Ns = Ntot
  !
  DimUp  = binomial(Ns,Nup)
  DimDw  = binomial(Ns,Ndw)
  Nsparse= Ns
  if(mpi_master)print*,"Dimensions =",DimUp,DimDw
  !
  Nprint=Neigen
  !
  call set_MpiComm(comm)
  call Setup_HxV(Nup,Ndw)  
  allocate(Eval(1),Evec(DimUp*mpiQdw,1))

  if(mpi_master)print*,"LANCZOS DIAG D - MPI:"


  do icycle=1,Ncycles
     !
     ! build symmetric sparse matrix
     !
     call Build_HxV()
     !
     ! LANCZOS diagonalization
     !
     Nitermax=512!min(Dim,512)
     Eval=0d0
     Evec=0d0  
     Nlanc=Nitermax     
     if(mpi_master)call start_timer
     t_start = MPI_Wtime()
     call sp_lanc_eigh(comm,mpi_HxVdirect,Nlanc,Eval(1),Evec(:,1),iverbose=.true.)
     if(mpi_master)call stop_timer
     t_end = MPI_Wtime()
     if(mpi_master)then
        do i=1,Nprint
           write(*,"(2F28.15,ES28.15)")Eval(i)
        end do
        open(100,file="mpi_lanczos_Eval_D.dat",position="append")
        do i=1,Nprint
           write(100,*)Eval(i)
        end do
        close(100)
        close(100)
        open(100,file="mpi_lanczos_time.dat",position="append")
        write(100,*)t_end-t_start
        close(100)
        print*,""
     endif
  enddo
  deallocate(Eval,Evec)  
  call Finalize_MPI()


end program test_MPI_LANCZOS_D








