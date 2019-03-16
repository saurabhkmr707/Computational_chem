program deepwell
implicit none

integer n,i,j,k,iplusj,nrot,ifail,lwork

real*8 smatrix, hmatrix, diag, x, res, pi, eigvec
allocatable::smatrix(:,:),hmatrix(:,:),diag(:), eigvec(:,:)
real*8 l

real*8 amat,a1mat, hmat,work
allocatable:: amat(:,:),a1mat(:,:),hmat(:,:),work(:)

write (*,*) 'Give nr. of basis states'
read (*, *) n

allocate(smatrix(0:n-1,0:n-1),hmatrix(0:n-1,0:n-1),diag(0:n-1),eigvec(0:n-1,0:n-1))
allocate(amat(0:n-1,0:n-1),a1mat(0:n-1,0:n-1),hmat(0:n-1,0:n-1),work(64*n))

pi = 4.0*datan(1.d0)
l = 2. !check here

do i=0,n-1
  do j=0,n-1
    iplusj = i+j
    if ( i .eq. j) then
    	smatrix(i,j) = 1.
    	hmatrix(i,j) = -(j+1)*(j+1)*pi*pi/(l*l) - 1./l
    else
    	smatrix(i,j) = 0.
    	hmatrix(i,j) = - (sin((i+1)*pi/l-(j+1)*pi/l))/(l*((i+1)*pi/l-(j+1)*pi/l))
    end if
  ENDDO
ENDDO

!diagonalize S-mat first
lwork=64*n
call dsyev('v','u',n,smatrix,n,diag,work,lwork,ifail)


do i = 0, n-1
 do j = 0, n-1
  amat(i,j) = smatrix(i,j)/sqrt(diag(j))
 enddo
enddo





hmat=matmul(transpose(amat),matmul(hmatrix,amat))

!diagonalize h-mat 
call dsyev('v','u',n,hmat,n,diag,work,lwork,ifail)

!carry out Av=c
a1mat=matmul(amat,hmat)

! Output the variational eigenvalues to the screen  together with the exact ones
! exact = n^2*pi^2/length^2

write (6,*) 'Variational     Exact'
do i=0, n-1
  write (6,'(5F12.4)') diag(i), (i+1)*(i+1)*pi*pi/4.d0
enddo

end

