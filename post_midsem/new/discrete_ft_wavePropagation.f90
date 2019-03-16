module m
 contains
 subroutine x_grid(x)
  implicit none
  real,dimension(256) :: x
  integer i 
  x(1) = -2.
  
  do i = 2,256
   x(i) = x(i-1)+0.02 
  end do
   
 end subroutine x_grid
 
 
 subroutine psix(x,psi)
  implicit none
  real,dimension(256)::x
  complex,dimension(256) :: psi
  real alpha,pi,k0,x0
  complex io
  integer i 
  io = cmplx(0.,1.)
  alpha = 20.
  k0 = 20.
  x0 = -0.5
  pi = 3.14
  
  do i = 1,256
   psi(i) = (((2.*alpha)/pi)**(1./4.))*exp(io*k0*(x(i)-x0))*exp(-alpha*((x(i)-x0)**2.))
  end do
 end subroutine psix
 
 subroutine potential (v,x)
  implicit none
  real,dimension(256) :: v,x
  integer i
  do i = 1,256
	if (x(i)  .ge. 0. .and. x(i) .le. 0.5 ) then
	v(i) = 0.1
	else
	v(i) = 0.
	end if
  end do
 end subroutine potential
 
 subroutine k_values(k,ksq,nx)
 implicit none
 real,dimension(256) ::  k,ksq
 integer nx,i
 real pi
 pi = 3.14
 do i = 1,nx
  if (i .le. nx/2) then
   k(i) = (2*pi*(i-1))/4.
  else
   k(i) = (2*pi*(i-nx-1))/4.
  end if
 end do 
 
  do i =1,256
   ksq(i) = (k(i)*k(i))
  end do  
 end subroutine k_values

 subroutine ft(phi,k,x,temp)
  implicit none
  complex,dimension(256) :: phi,temp
  real,dimension(256) :: k,x
  integer m,j
  complex i0
  i0 = cmplx(0.,1.)
  
  do m =1,256
   phi(m) = cmplx(0.,0.)
   do j =1,256
    phi(m) = phi(m) + temp(j)*exp(-i0*k(m)*x(j))
   end do
   phi(m) = phi(m)/(sqrt(256.))
  end do
  
  
 end subroutine ft
 
 subroutine ift(psi2,k,x,phi)
  implicit none
  complex,dimension(256) :: psi2,phi
  real,dimension(256) :: k,x
  integer m,j
  complex i0
  i0 = cmplx(0.,1.)
  
  do j = 1,256
  psi2(j) = cmplx(0.,0.)
   do m = 1,256
    psi2(j) = psi2(j) + phi(m)*exp(i0*k(m)*x(j))
   end do
  psi2(j) = psi2(j)/(sqrt(256.))
  end do
  
 end subroutine ift

 
end module m


program main
 use m
 implicit none
 integer i,j,t,nx
 real,dimension(256) :: x,v,k,ksq
 complex,dimension(256) :: psi,temp,phi,psi2
 complex eye
 character (len=100):: new
 eye = cmplx(0.,1.)
 nx = 256
 
 call x_grid(x)
 call potential(v,x)
 call psix(x,psi)
 call k_values(k,ksq,nx) !most probably correct upto here

 do t = 1,5000
 write(new,1) t+10000
 1 format('',i5)
 open(unit = 2,file = new)
  do i =1,256
   temp(i) = exp(-eye*(V(i)/2.)*0.1)*psi(i)
  end do
  
  
  call fft(temp,256,1)
  
  
  do i = 1,256
   temp(i) = temp(i)*exp(-(eye*((k(i)*k(i))/(2.*14500.)))*0.1)
  end do
  
  
  
  call fft(temp,256,-1)
  
  
  do i = 1,256
   psi(i) = temp(i)*exp(-eye*(V(i)/2.)*0.1)
  end do

  do i = 1,256
   write(2,*) x(i),abs(psi(i))**2.
  end do
  close(2)
  print*,t
  
 end do 
  
 

 
end program main
