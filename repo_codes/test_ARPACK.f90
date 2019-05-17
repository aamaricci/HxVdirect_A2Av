program testLANCZOS_D
  USE COMMON_VARS
  USE SETUP
  USE MATVEC_PRODUCT
  !
  !
  USE SF_PARSE_INPUT
  USE SF_LINALG
  USE SF_TIMER
  !
  USE SF_SP_LINALG_
  !
  implicit none
  !Dimension
  integer                        :: nup,ndw,Dim
  integer                        :: Ntot,Nsparse,Ncycles
  !Matrix:
  real(8),allocatable            :: Hmat(:,:)
  integer                        :: Nprint,N,Nloc
  real(8),allocatable            :: Eval(:),Eval_old(:)
  real(8),allocatable            :: Evec(:,:)
  !Lanczos
  integer                        :: neigen,nitermax,niter,nblock,Nlanc
  integer                        :: i,j,k
  integer                        :: iup,idw,jup,jdw,ncv

  !
  call parse_cmd_variable(ntot,"ntot",default=10)
  call parse_cmd_variable(nup,"nup",default=5)
  call parse_cmd_variable(ndw,"ndw",default=5)
  call parse_cmd_variable(Neigen,"NEIGEN",default=1)
  call parse_cmd_variable(NCV,"NCV",default=10)
  call parse_cmd_variable(Ncycles,"NCYCLES",default=1)
  !
  Ns = Ntot
  !
  DimUp  = binomial(Ns,Nup)
  DimDw  = binomial(Ns,Ndw)
  Nsparse= Ns
  print*,"Dimensions =",DimUp,DimDw


  Nprint=Neigen

  !
  ! build symmetric sparse matrix
  !
  print*,"SETUP MATRIX BOUND:"
  call Setup_HxV(Nup,Ndw)
  !


  !
  ! ARPACK LANCZOS diagonalization
  !
  print*,"ARPACK DIAG :"
  Nitermax=512!min(Dim,512)
  Nblock=ncv*Neigen !min(ncv*Neigen,Dim)
  allocate(Eval(Neigen),Evec(DimUp*DimDw,Neigen))
  Eval=0d0
  Evec=0d0
  call start_timer
  call cpu_time(t_start)
  call sp_eigh(HxVdirect,Neigen,Nblock,Nitermax,Eval,Evec)
  call cpu_time(t_end)
  call stop_timer
  !
  do i=1,Nprint
     write(*,"(2F28.15,ES28.15)")Eval(i)
  end do
  open(100,file="arpack_Eval_D.dat")
  do i=1,Nprint
     write(100,*)Eval(i)
  end do
  close(100)
  !
  open(100,file="arpack_time.dat")
  write(100,*)t_end-t_start
  close(100)
  !
  deallocate(Eval,Evec)
  print*,""
  print*,""

  !
  ! LANCZOS diagonalization
  !
  print*,"LANCZOS DIAG :"
  Nitermax=512!min(Dim,512)
  allocate(Eval(1),Evec(DimUp*DimDw,1))
  Eval=0d0
  Evec=0d0
  Nlanc=Nitermax
  call start_timer
  call cpu_time(t_start)
  call sp_lanc_eigh(HxVdirect,Nlanc,Eval(1),Evec(:,1))
  call cpu_time(t_end)
  call stop_timer
  !
  do i=1,Nprint
     write(*,"(2F28.15,ES28.15)")Eval(i)
  end do
  open(100,file="lanczos_Eval_D.dat")
  do i=1,Nprint
     write(100,*)Eval(i)
  end do
  close(100)
  !
  !
  open(100,file="lanczos_time.dat")
  write(100,*)t_end-t_start
  close(100)
  !
  deallocate(Eval,Evec)  
  print*,""
  print*,""


end program testLANCZOS_D






