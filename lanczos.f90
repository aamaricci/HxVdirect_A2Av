
subroutine eigensolve(MpiComm,MatVec,lambda,V)                                        !
  integer             :: MpiComm
  interface 
     subroutine MatVec(Nloc,vin,vout)
       integer        :: Nloc
       real(8)        :: vin(Nloc)
       real(8)        :: vout(Nloc)
     end subroutine MatVec
  end interface
  real(8),intent(out) :: lambda(:), V(:,:)

  integer             :: i, info, n
  real(8),allocatable :: T(:,:), Q(:,:), work(:)

  n = size(lambda)
  allocate(T(3, n), Q(n, n), work(2 * n - 2))

  call lanczos(A, T, V)

  call dstev('V', n, T(2, 1:n), T(3, 1:n-1), Q, n, work, info)

  V = matmul(V, Q)

  do i = 1, n
     V(:,i) = V(1,i) / dabs(V(1, i)) * V(:,i)
  enddo

  lambda = T(2, 1:n)

end subroutine eigensolve




!--------------------------------------------------------------------------!
subroutine lanczos(A, T, Q)                                                !
  !--------------------------------------------------------------------------!
  ! Perform n steps of the Lanczos process on the matrix A, putting the      !
  ! coefficients into the 3 x n array T and the vectors into the array Q.    !
  !--------------------------------------------------------------------------!
  ! input/output variables
  class(sparse_matrix_interface), intent(in) :: A
  real(dp), intent(out) :: T(:,:), Q(:,:)
  ! local variables
  integer :: i, k, n
  real(dp) :: w(A%nrow), alpha, beta

  n = size(T, 2)

  T = 0.0_dp
  Q = 0.0_dp
  w = 0.0_dp

  !-----------------------------------
  ! First step of the Lanczos process
  call init_seed()

  ! Make a random unit vector
  call random_number(Q(:,1))
  Q(:,1) = 2 * Q(:,1) - 1
  Q(:,1) = Q(:,1) / dsqrt(sum(Q(:,1) * Q(:,1)))

  ! Fill out the entries in T, Q
  call A%matvec(Q(:,1), w)
  alpha = sum(Q(:,1) * w)
  w = w - alpha * Q(:,1)
  beta = dsqrt(sum(w * w))
  Q(:,2) = w / beta
  T(2,1) = alpha
  T(3,1) = beta
  T(1,1) = beta


  !-----------------------------------------
  ! Subsequent steps of the Lanczos process
  do i = 2, n - 1
     call A%matvec(Q(:,i), w)
     alpha = sum(Q(:,i) * w)
     w = w - alpha * Q(:,i) - beta * Q(:,i - 1)

     ! Re-orthogonalize the current vector against previous Lanczos vectors
     ! to avoid floating point goofs.
     do k = 1, i - 2
        w = w - sum(Q(:,k) * w) * Q(:,k)
     enddo

     beta = dsqrt(sum(w * w))
     Q(:,i+1) = w / beta
     T(2,i) = alpha
     T(3,i) = beta
     T(1,i) = beta
  enddo

  !---------------------------------------
  ! Last iteration of the Lanczos process
  call A%matvec(Q(:,n) , w)
  T(2, n) = sum(Q(:,n) * w)

end subroutine lanczos
