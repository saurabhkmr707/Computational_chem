program main
 implicit none
 integer i,j,t
 real*8 dx,x0,sigma,beta,v1,beta1,x1,x2,v2,v3,vlower,v5,v4,vasym,beta2,beta3,beta4,x4,x3, &
 pi,mu,dt,xmask,xmax,dxmask,p0,de
 real*8,dimension(2048) :: x,k,v1ad,v2ad,v11,v22,v12,d
 complex*16,dimension(2048) :: psi1,psi2
 complex*16,dimension(2,2048) :: psi
 real*8,dimension(250) :: e,k1,s11,s111,s112
 complex*16,dimension(2048,1000) :: psitemp
 complex*16 i0
 complex*16,dimension (2,2) :: va,vb
 character(len = 100) p
 
 
 de = 0.001
 i0 = cmplx(0.,1.)
 dx = 90./2048.
 x0 = 9.
 x1 = -4.364721325998*0.01
 sigma = 0.3 
 beta = 1./(4.*sigma*sigma)
 beta1 = 5.5
 v2 = 4.79833373*10.**(-3.)
 vasym = 3.61196179*0.1
 v3 = 9.8998917754*0.1
 v1 = 4.0167971782296*0.01 
 vlower = 0.0
 x3 = -7.6042693477*0.1
 beta3 = 2.3471780470
 v5 = 7.9781762366*0.1
 v4 = 1.122019*0.01
 x2 = 5.0012635420*0.01
 beta2 = 4.9818195151
 beta4 = 1.0487590725
 x4 = 8.1790045179*0.1
 pi = 3.1415926
 mu = 3474.057
 dt = 8.
 xmax = 45.
 
!x-grid
 x(1) = -45.
 do i =2,2048
 	x(i) = x(i-1)+dx	
 enddo
 
 !k-grid
 do i = 1,2048
	  if (i .le. 2048./2.) then
	   k(i) = (2.*pi*(i-1))/(2048.*dx)
	  else
	   k(i) = (2.*pi*(i-2048.-1))/(2048.*dx)
	  end if
 end do 
 
 
 !------------------
 !v grid
 
 do i =1,2048 !v1ad
 	v1ad(i) = (v1*exp(beta1*(x(i)-x1)))/(1.+exp(beta1*(x(i)-x1)))**2. &
 	+ (v2*exp(beta1*(x(i)-x1)))/(1.+exp(beta1*(x(i)-x1)))
 enddo
 
 do i =1,2048 !v2ad
 	v2ad(i) = vasym - (v3*exp(beta2*(x(i)-x2)))/(1.+exp(beta2*(x(i)-x2)))**2. &
 	- (v4*exp(beta2*(x(i)-x2)))/(1.+exp(beta2*(x(i)-x2))) &
 	- (v5*exp(beta3*(x(i)-x3)))/(1.+exp(beta3*(x(i)-x3)))**2. - vlower
 enddo
 
 do i =1,2048 !v11
 	v11(i) = (1-f(x(i)))*v1ad(i) + f(x(i))*v2ad(i)
 end do
 
 do i =1,2048 !v22
 	v22(i) = f(x(i))*v1ad(i) + (1.-f(x(i)))*v2ad(i)
 end do 
 
 do i =1,2048 !v12 (was correct upto here v12)
 	v12(i) = -sqrt(f(x(i))*(1-f(x(i))))*(v2ad(i)-v1ad(i))
 end do 
 
 
 
 !---------------------------------------------
 !inital wave packet
 
  !p0
 p0 = -sqrt(abs(2*mu*(0.029-v11(1230))))
 
 
 do i =1,2048 !initial psi
 	psi1(i)= ((1./(2*pi*sigma*sigma))**(0.25))*exp(-beta*(x(i)-x0)**2.)*exp(i0*p0*(x(i)-x0))
 end do  !(correct upto here)
 
!----------------------------------------------

 do i =1,2048 !2d psi
 	psi(1,i) = psi1(i)
 	psi(2,i) = cmplx(0.,0.)
 	psi2(i) = cmplx(0.,0.)
 end do 
 
 !-----------------------
 
 do i=1,2048
	if(x(i) .ge. 0 ) then
		xmask = 10.
	else
		xmask = -30.
	endif
	if(abs(x(i)).ge. abs(xmask)) then 
		dxmask = xmax - abs(xmask)
		d(i) = sin(pi*(abs(xmask)+ dxmask - abs(x(i)))/(2*dxmask))
	else
		d(i) = 1.0
	endif
	write(20,*) x(i), d(i)
 enddo

 

 do t = 1,1000
 	do i =1,2048 !multiplied va
 		va(1,1) = exp(-i0*v11(i)*dt/4.)
 		va(2,2) = exp(-i0*v22(i)*dt/4.)
 		va(1,2) = cmplx(0.,0.)
 		va(2,1) = cmplx(0.,0.)
 		psi(:,i) = matmul(va,psi(:,i))
 	end do 
 	
 	do i =1,2048 !calculating vb converting into cos and sin form
 		vb(1,1) = cos(v12(i)*dt/2.)
 		vb(1,2) = -i0*sin(v12(i)*dt/2.)
 		vb(2,1) = -i0*sin(v12(i)*dt/2.)
 		vb(2,2) = cos(v12(i)*dt/2.)
 		
 		psi(:,i) = matmul(vb,psi(:,i))
 	end do 
 	
 	do i =1,2048 !multiplied va
 		va(1,1) = exp(-i0*v11(i)*dt/4.)
 		va(2,2) = exp(-i0*v22(i)*dt/4.)
 		va(1,2) = cmplx(0.,0.)
 		va(2,1) = cmplx(0.,0.)
 		psi(:,i) = matmul(va,psi(:,i))
 	end do 
 	
 	
 	call fft(psi(1,:),2048,1)
 	call fft(psi(2,:),2048,1)
 	
 	
 	do i =1,2048 !multiplied exp(-i*T*dt)
 		psi(1,i) = psi(1,i)*exp(-i0*k(i)*k(i)*dt/(2.*mu))
 		psi(2,i) = psi(2,i)*exp(-i0*k(i)*k(i)*dt/(2.*mu))
 	end do
 	
 	
 	
 	call fft(psi(1,:),2048,-1)
 	call fft(psi(2,:),2048,-1)
 	
 	psi = psi/2048.
 	
 	
 	do i =1,2048 !multiplied va
 		va(1,1) = exp(-i0*v11(i)*dt/4.)
 		va(2,2) = exp(-i0*v22(i)*dt/4.)
 		va(1,2) = cmplx(0.,0.)
 		va(2,1) = cmplx(0.,0.)
 		psi(:,i) = matmul(va,psi(:,i))
 	end do 
 	
 	do i =1,2048 !calculating vb converting into cos and sin form
 		vb(1,1) = cos(v12(i)*dt/2.)
 		vb(1,2) = -i0*sin(v12(i)*dt/2.)
 		vb(2,1) = -i0*sin(v12(i)*dt/2.)
 		vb(2,2) = cos(v12(i)*dt/2.)
 		
 		psi(:,i) = matmul(vb,psi(:,i))
 	end do 
 	
 	do i =1,2048 !multiplied va
 		va(1,1) = exp(-i0*v11(i)*dt/4.)
 		va(2,2) = exp(-i0*v22(i)*dt/4.)
 		va(1,2) = cmplx(0.,0.)
 		va(2,1) = cmplx(0.,0.)
 		psi(:,i) = matmul(va,psi(:,i))
 		psitemp(i,t) = psi(1,i)
 	end do 

 end do 
 

 ! e-grid
 e(1) = 0.02
 do i = 2,250
 	e(i) = e(i-1) + de 
 end do 
 
 do i =1,250
 	k1 = -sqrt(2*mu*(e(i)-v11(1936)))
 end do 
 
 do i =1,250
 	s112(i) = cmplx(0.,0.)
 	do t = 1,1000
 		s112(i) = s112(i) + psitemp(1936,t)*exp(i0*e(i)*8.*t)*8.  
 	enddo	
 end do 
 
 do i = 1,250
 	s111(i) = sqrt(abs(k1(i)*k1(i)))*exp(i0*k1(i)*x(1936))/(2.*pi*mu)
 end do 
 
 do i =1,250
 	s11(i) = s111(i)*s112(i)
 enddo
 

	open(unit=3, file='s11' )
	do i=1,250
		write(3,*) e(i),abs(s11(i)**2.) 
	enddo
	close(3)
	
	open(unit = 4,file = 'v11')
	do i =1,2048
		write(4,*) x(i),v11(i)
	end do 
	close(4)
 
 
 contains 
 real*8 function f(x)
  implicit none
  real*8 x
  f = 0.5*(1-tanh(1.0487590725*(x-8.1790045179*0.1)))
 end function f
end program main

