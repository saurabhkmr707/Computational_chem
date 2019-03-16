program splitfft
implicit none

complex, dimension(128):: psix,psiy,psi
complex,dimension(128,128) :: psixy
real, dimension(128):: vx,vy,x,y
integer i,j,h
real a,x0,k0,m,xmin,dx,dt,ymin,dy,y0
complex:: iota =cmplx(0.,1.)
character*5 p
real:: pi= acos(-1.0)
real,dimension(128,128) :: vxy

a=20.0
x0=-0.5
k0=50.0
m=14500.0
xmin=-2.0
dx=0.04
dy = 0.04
dt=0.1
ymin = -2.0


y(1) = ymin
do i = 2,128
	y(i) = y(i-1) + dy
end do

x(1)=xmin
do i=2,128
	x(i) = x(i-1) + dx
enddo


do i=1,128
	if(x(i) .ge. 0.0 .and. x(i) .le. 0.5) then 
                vx(i) = 0.8
        else
                vx(i) = 0.0
        endif
enddo

do i=1,128
	if(y(i) .ge. 0.0 .and. y(i) .le. 0.5) then 
                vy(i) = 0.8
        else
                vy(i) = 0.0
        endif
enddo
open (unit = 6001,file = 'potential.txt')

do i = 1,128
  do j=1,128
	vxy(i,j) = vx(i)*vy(j)
	write(6001,*) x(i),y(j),vxy(i,j)
  end do 	
  write(6001,*)
end do


close(6001) 

do h=1,128
	psix(h) = sqrt(sqrt(2.*a/pi)) * exp(-a*(x(h)-x0)**2) * exp(iota*(k0*(x(h)-x0)))
	psiy(h) = sqrt(sqrt(2.*a/pi)) * exp(-a*(y(h)-y0)**2) * exp(iota*(k0*(y(h)-y0)))
enddo



open(1,file='10001')
do i=1,128
	do j =1,128
		psixy(i,j) = psix(i)*psiy(j)
		write(1,*) x(i),y(j), (abs(psixy(i,j))**2)
	end do
 write(1,*)
enddo
close(1)

do j=2,5000
	
	
	do i = 1,128
		do h = 1,128
			psi(h) = psixy(i,h)
			psi(h) = psi(h)*exp(-iota*vxy(i,h)*dt/2.)
		end do
		
		call fft(psi,128,+1)
		do h = 1,128
		 psi(h) = exp(-iota*dt*k(h)*k(h)/(2.*m)) * psi(h)
		end do
		
		call fft(psi,128,-1)
		
		do h = 1,128
			psi(h) = psi(h)/128.
		end do

		do h=1,128
			psi(h) = exp(-iota*vxy(i,h)*dt/2.) * psi(h)
		enddo
		
		do h = 1,128
			psixy(i,h) = psi(h)
		end do

		
	end do 
	
	
	do h = 1,128
		do i = 1,128
			psi(i) = psixy(i,h)
			psi(i) = psi(i)*exp(-iota*vxy(i,h)*dt/2.)
		end do
		
		call fft(psi,128,+1)
		do i = 1,128
		 psi(i) = exp(-iota*dt*k(i)*k(i)/(2.*m)) * psi(i)
		end do
		
		call fft(psi,128,-1)
		
		do i = 1,128
			psi(i) = psi(i)/128.
		end do

		do i=1,128
			psi(i) = exp(-iota*vxy(i,h)*dt/2.) * psi(i)
		enddo
		
		do i = 1,128
			psixy(i,h) = psi(i)
		end do

		
	end do 
	
	
	if(mod(j,100).eq.1) then
	write(p,'(I5)')10000+j
	open(unit=j, file=p )
	do i=1,128
		do h =1,128
		
		write(j,*) x(i),y(h), (abs(psixy(i,h))**2)
		end do
 	write(j,*)
	enddo
	endif
	print*,'writing file: ',j 
	
	
enddo

contains

real function k(h1)
integer h1
if (h1 .le. 64)then  
k = 2.*pi*(h1-1.)/(128.*dx)
else
k = 2*pi*(h1-1-128.)/(128.*dx)
endif
return
end function k

end
