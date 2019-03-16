program main
 implicit none
 real x0,alpha,p0,m,xmin,dx,dt,pi
 real,dimension(256) :: x(0:255),v(0:255),k(0:255)
 complex,dimension(256) :: psi0(0:255),psi1(0:255),psitemp(0:255),psi(0:255)
 complex i0
 integer i,j,t
 character(len = 100) new
 
 i0 = cmplx(0.,1.)
 pi = 3.14159
 alpha = 20. 
 dt = 0.1
 m = 14500.
 ! x-grid
 xmin = -2.0
 dx = 0.02
 x0 = -0.5
 p0 = 20.
 x(0) = xmin 
 do i = 1,255
 	x(i) = x(i-1) +dx
 end do	
 
 ! v-grid
 do i = 0,255
 	if ( x(i) .lt. 0.) then
 		v(i) = 0.
 	else 
 		v(i) = 0.1
 	end if
 end do 
 
 !k-grid
 do i = 0,255
 	if ( i .lt. 256./2-1.) then
 		k(i) = 2.*pi*i/(256.*dx)
 	else 
 		k(i) = 2.*pi*(i-256.-1)/(256.*dx)
 	end if
 end do 
 
 !psi calculation
  
  !inti psi
  do i =0,255
  	psi0(i) = ((2.*alpha/pi)**0.25)*exp(i0*k(i)*(x(i)-x0))*exp(-alpha*(x(i)-x0)**2.)
  	psitemp(i) = psi0(i)
  end do
  
  call dft(psitemp,k,x)
  
  do i = 0,255
  	psi1(i) = psi0(i) + i0*dt*((1./(2.*m))*psitemp(i)- v(i)*psi0(i))
  	psitemp(i) = psi1(i)
   end do 
  
  do t = 1,5000
	  write(new,1) t+10000
	  1 format('',i5)
  	call dft(psitemp,k,x)
  	do i = 0,255
  		psi(i) = psi0(i) + 2.*i0*dt*((1./(2.*m))*psitemp(i) - v(i)*psi1(i))
  	end do 
  	
  	open(unit = t,file = new)
  	do i = 0,255
  		write(t,*) x(i),abs(psi(i))**2.
  	end do 
  	close(t)
  	
  	do i =0,255
  		psi0(i)  = psi1(i)
  		psitemp(i) = psi(i)
  		psi1(i) = psi(i)
  		  
  	end do 
  end do  
  
  
 
 
 
 open (unit = 16000,  file = 'test')
 do i =  0,255
 	write(1,*) x(i),(abs(psi1(i)))**2.
 end do 	
 close(16000)
 
 contains 
 subroutine dft(psi,k,x)
 implicit none
 real,dimension(256) :: x(0:255),v(0:255),k(0:255)
 complex,dimension(256) :: psi(0:255),phi(0:255)
 complex i0
 integer i,j
 i0 = cmplx(0.,1.)
 
 
 do i = 0,255
 	phi(i) = cmplx(0.,0.)
 	do j = 0,255
 		phi(i) = phi(i) + (1./sqrt(256.))*psi(j)*exp(-i0*k(i)*x(j))
 	end do 
 end do 
 
 do i = 0,255
 	phi(i) = -phi(i)*k(i)*k(i)
 end do
 
 do i = 0,255
 	psi(i) = cmplx(0.,0.)
 	do j = 0,255
 		psi(i) = psi(i) + (1./sqrt(256.))*phi(j)*exp(i0*k(j)*x(i))
 	end do 
 end do 
 return
 end subroutine dft
end program main
