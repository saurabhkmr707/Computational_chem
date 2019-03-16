program main
 implicit none
 real,dimension(256) :: x,k,v
 complex,dimension(256) ::  psi
 complex i0
 integer i,j,t
 real k0,pi,alpha,x0
 
 k0 = 20
 pi = 3.14
 i0 = cmplx(0.,1.)
 dx = 0.02
 dt = 0.1
 alpha = 20.
 x0 = -0.5
 
 !x grid
 x(1) = -2.
 do i =2,256
  x(i) = x(i-1) + dx
 end do
 
 !potential grid
 do i =2,256
  if (x(i)  .ge. 0. .and. x(i) .le. 0.5 ) then
	v(i) = 0.
	else
	v(i) = 0.1
   end if
 end do
 
 !k grid
 do i = 1,256
  if (i .le. (256./2.-1)) then
   k(i) = (2.*pi*(i-1))/4.
  else
   k(i) = (2.*pi*(i-256.-1.))/4.
  end if
 end do 
 
 
 
 ! initial psi
 
 do i = 1,256
  psi(i) = (((2.*alpha)/pi)**(1./4.))*exp(i0*k0*(x(i)-x0))*exp(-alpha*((x(i)-x0)**2.))
 end do 
 
 
 
end program main
