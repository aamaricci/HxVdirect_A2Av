!---------------------------------------------------------------------
!Purpose: use plain lanczos to get the groundstate energy
!---------------------------------------------------------------------
subroutine lanczos_eigh_d(MatVec,Nitermax,Egs,Vect,iverbose,threshold,ncheck)
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  integer                              :: Nitermax
  real(8)                              :: egs
  real(8),dimension(:)                 :: vect
  real(8),optional                     :: threshold
  integer,optional                     :: ncheck
  logical,optional                     :: iverbose
  !
  integer                              :: Ndim
  real(8),dimension(size(vect))        :: vin,vout
  integer                              :: iter,nlanc
  real(8),dimension(Nitermax+1)        :: alanc,blanc
  real(8),dimension(Nitermax,Nitermax) :: Z
  real(8),dimension(Nitermax)          :: diag,subdiag,esave
  real(8)                              :: a_,b_
  real(8)                              :: norm,diff
  integer                              :: i,ierr
  !
  if(present(iverbose))verb=iverbose
  if(present(threshold))threshold_=threshold
  if(present(ncheck))ncheck_=ncheck
  !
  Ndim = size(vect)
  !
  norm=dot_product(vect,vect)
  if(norm==0d0)then
     call mt_random(vect)
     vect=vect/sqrt(dot_product(vect,vect))
     if(verb)write(*,*)"LANCZOS_EIGH_D: random initial vector generated:"
  endif
  !
  !============= LANCZOS LOOP =====================
  !
  vin = vect
  vout= 0.d0
  alanc=0.d0
  blanc=0.d0
  nlanc=0
  !
  lanc_loop: do iter=1,Nitermax
     call lanczos_iteration_d(MatVec,iter,vin,vout,a_,b_)
     if(abs(b_)<threshold_)exit lanc_loop
     !
     nlanc=nlanc+1
     !
     alanc(iter) = a_ ; blanc(iter+1) = b_
     diag    = 0d0
     subdiag = 0.d0
     Z       = eye(Nlanc)
     diag(1:Nlanc)    = alanc(1:Nlanc)
     subdiag(2:Nlanc) = blanc(2:Nlanc)
     call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
     !
     if(nlanc >= Ncheck_)then
        esave(nlanc-(Ncheck_-1))=diag(1)
        if(nlanc >= (Ncheck_+1))then
           diff=esave(Nlanc-(Ncheck_-1))-esave(Nlanc-(Ncheck_-1)-1)
           if(verb)write(*,*)'iter, E0, deltaE = ',iter,diag(1),diff
           if(abs(diff).le.threshold_)exit lanc_loop
        endif
     endif
     !
  enddo lanc_loop
  if(nlanc==nitermax)print*,"LANCZOS_SIMPLE: reach Nitermax"
  !
  !============== END LANCZOS LOOP ======================
  !
  diag    = 0d0
  subdiag = 0.d0
  Z       = eye(Nlanc)
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
  !
  !Get the Eigenvalues:
  egs = diag(1)
  !
  !Get the Eigenvector:
  vin =vect
  vout=0.d0
  vect=0.d0
  do iter=1,Nlanc
     call lanczos_iteration_d(MatVec,iter,vin,vout,alanc(iter),blanc(iter))
     vect = vect + vin*Z(iter,1)
  end do
  norm=sqrt(dot_product(vect,vect))
  vect=vect/norm
  !
  if(verb)then
     call MatVec(Ndim,vect,vout)
     write(*,*)"|H*v-E*v|=",sum(abs(vout-egs*vect))/Ndim
  endif
  Nitermax=Nlanc
  !
end subroutine lanczos_eigh_d





















!---------------------------------------------------------------------
!Purpose: plain homebrew lanczos iteration (no orthogonalization)
!note: the a,b variables are real, even in the complex matrix case
!to understand why check out the Gollub-Van Loan textbook.
!a it is easy: hermiticity->diag\in\RRR
!b: is fixed by requiring |b|^2 = <v,v> thus you can only fix the 
!the absolute value. A lemma shows that the phase can be chosen 
!identically zero
!---------------------------------------------------------------------
subroutine lanczos_iteration_d(MatVec,iter,vin,vout,alfa,beta)
  interface
     subroutine MatVec(Nloc,vin,vout)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: vin
       real(8),dimension(Nloc) :: vout
     end subroutine MatVec
  end interface
  integer                                    :: iter
  real(8),dimension(:),intent(inout)         :: vin
  real(8),dimension(size(vin)),intent(inout) :: vout
  real(8),dimension(size(vin))               :: tmp
  real(8),intent(inout)                      :: alfa,beta
  integer                                    :: nloc
  real(8)                                    :: norm
  !
  nloc=size(vin)
  !
  if(iter==1)then
     norm=sqrt(dot_product(vin,vin))
     if(norm==0.d0)stop "LANCZOS_ITERATION_D: norm =0!!"
     vin    = vin/norm
  else
     tmp = vin
     vin = vout/beta
     vout= -beta*tmp
  endif
  call MatVec(nloc,vin,tmp)
  vout    = vout + tmp
  alfa = dot_product(vin,vout)
  vout    = vout - alfa*vin
  beta = sqrt(dot_product(vout,vout))
  ! if(iter==1)then
  !    norm=sqrt(dot_product(vin,vin))
  !    if(norm==0.d0)stop "LANCZOS_ITERATION_D: norm =0!!"
  !    vin=vin/norm
  !    b=0d0
  ! end if
  ! call MatVec(nloc,vin,tmp)
  ! tmp = tmp-b*vout
  ! a   = dot_product(vin,tmp)
  ! tmp = tmp-a*vin
  ! b   = sqrt(dot_product(tmp,tmp))
  ! vout= vin
  ! vin = tmp/b
end subroutine lanczos_iteration_d
