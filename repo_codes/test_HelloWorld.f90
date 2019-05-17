program test_MPI_LANCZOS_D
  USE SCIFOR
  USE MPI
  implicit none
  !MPI:
  integer                          :: comm,irank
  integer                          :: mpi_rank
  integer                          :: mpi_size
  integer                          :: mpi_ierr
  logical                          :: mpi_master,bool
  integer                          :: mpi_Q,Q,mpi_Istart
  integer                          :: mpi_R,R,mpi_Iend
  integer                          :: mpi_Chunk,mpi_Ishift
  real(8) :: t_start,t_end
  real(8) :: Matrix(1024,1024),eval(1024)

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)

  mpi_rank   = get_rank_MPI(comm)
  mpi_size   = get_size_MPI(comm)
  mpi_master = get_master_MPI(comm)

  call Barrier_MPI(comm)

  do irank=0,mpi_size-1
     call Barrier_MPI(comm)
     if(irank==mpi_rank)then
        write(*,*)"Node",mpi_rank," of ",mpi_size,": active"     
        call mt_random(matrix)
        call eigh(matrix,eval)
        write(*,*)"Node",mpi_rank,":",eval(1)
     end if
  enddo
  call Barrier_MPI(comm)

  
  call Finalize_MPI()
end program test_MPI_LANCZOS_D
