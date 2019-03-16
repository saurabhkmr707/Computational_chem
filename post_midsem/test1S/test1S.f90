
program main
 implicit none
 integer i,j,t,nx
 real,dimension(256) :: x,v,k
 complex,dimension(256) :: psi,temp,phi,psi2
 complex eye,io
 real alpha,k0,x0,pi
 character (len=100):: new
 io = cmplx(0.,1.)
  alpha = 20.
  k0 = 50.
  x0 = -0.5
  pi = 3.14
 eye = cmplx(0.,1.)
 nx = 256
 
 x(1) = -2.
  
  do i = 2,256
   x(i) = x(i-1)+0.02 
  end do
  
  
 do i = 1,256
	if ((x(i)  .ge. 0.) .and. (x(i) .le. 0.5) ) then
	v(i) = 0.1 
	else
	v(i) = 0.
	end if
  end do
  
  open(unit = 6000,file = 'potential.txt')
   do i = 1,256
    write(6000,*) x(i),v(i)
   end do 
  close(6000)
  
 do i = 1,256
   psi(i) = (((2.*20.)/pi)**(1./4.))*exp(io*k0*(x(i)-x0))*exp(-alpha*((x(i)-x0)**2.))
  end do
  
 do i = 1,256
  if (i .le. (256./2 -1.)) then
   k(i) = (2*pi*(i-1))/(265.*0.02)
  else
   k(i) = (2*pi*(i-256.-1))/((265.*0.02))
  end if
 end do 
 
 open(unit = 6001,file = 'kval.txt')
   do i = 1,256
    write(6001,*) x(i),k(i)
   end do 
  close(6001)

print*,'hello'

 do t = 1,5000
 write(new,1) t+10000
 1 format('',i5)
 open(unit = 15,file = new)
  do i =1,256
   psi(i) = exp(-eye*(v(i)/2.)*0.1)*psi(i)
  end do
  
  
  call fft(psi,256,1)
  
  
  do i = 1,256
   psi(i) = psi(i)*exp(-eye*((k(i)*k(i))/(2.*14500.))*0.1)
  end do
  
  
  
  call fft(psi,256,-1)
  
  do i =1,256
   psi(i) = psi(i)/256.
  end do
  
  
  do i = 1,256
   psi(i) = psi(i)*exp(-eye*(V(i)/2.)*0.1)
  end do

  do i = 1,256
   write(15,*) x(i),(abs(psi(i))**2.)
  end do
  
  print*,t
  
 end do 
  
 

 
end program main
