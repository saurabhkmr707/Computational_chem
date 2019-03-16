program tdseft
implicit none

 complex(kind = 16), dimension(20,5000):: psi
real(kind = 16), dimension(20):: v,x,z
real(kind = 16), dimension(5000):: t
integer i,j
real(kind = 16) a,x0,k0,m,xmin,dx,dt
character*5 p

open(10, file='psi')

a=20.0
x0=-0.5
k0=20.0
m=14500
xmin=-2.0
dx=0.02
dt=0.1

x(1)=xmin
V(1)=0
do i=2,20
	x(i) = x(i-1) + dx
	if(x(i).ge.0) then 
                v(i) = 1
        else
                v(i) = 0
        endif
enddo



z(1) = 1.89*exp(-a*(x(1)-x0)**2)
	psi(1,1) = z(1) * exp(cmplx(0,1)*(k0*(x(1)-x0)))
	psi(1,2) = psi(1,1) + dt*cmplx(0,1)*( (1/(2*m))*( (psi(3,1) - 2*psi(2,1) + psi(1,1))/dt**2) - v(i)*psi(1,1) )
	do j=3,5000
		psi(1,j) = psi(1,j-2) - (2*dt)*cmplx(0,1) * ( -1/(2*m)*dp(1,j-1) + v(1)*psi(1,j-1) )
	enddo

do i=2,20
	z(i) = 1.89*exp(-a*(x(i)-x0)**2)
	psi(i,1) = z(i) * exp(cmplx(0,1)*(k0*(x(i)-x0)))
	psi(i,2) = psi(i,1) + dt*cmplx(0,1)*( (1/(2*m))*dp(i,1) - v(i)*psi(i,1) )
enddo


do j=3,5000
	do i=2,20
		psi(i,j) = psi(i,j-2) - (2*dt)*cmplx(0,1) * ( -1/(2*m)*dp(i,j-1) + v(i)*psi(i,j-1) )
	enddo
enddo

do j=1,5000,100
	write(p,'(I5)')20000+j
	open(unit=j, file=p ) 
	do i=1,20
		write(j,*) x(i),(real(psi(i,j))**2+aimag(psi(i,j))**2)
	enddo
	close(j)
enddo

contains

complex(kind = 16) function dp(i1,j1)
integer i1,j1,m1,h
complex(kind = 16) dd
complex, dimension(20) :: phi,k
do m1=1,20
	if(m1.le.9) then
		k(m1) = (6.28/5.12)*m1
	else
		k(m1) =  (6.28/5.12)*(m1-50)
	endif
enddo
do m1=1,20
	phi(m1) = cmplx(0,0)
	do h=1,20
		phi(m1) = phi(m1) + (1./16)*psi(h,j1)*exp(-cmplx(0,1)*k(m1)*x(h)) 
	enddo
enddo

dd = cmplx(0,0)
do m1=1,20
	dd = dd - (k(m1)**2)*(1./16)*phi(m1)*exp(cmplx(0,1)*k(m1)*x(i1))
enddo
dp=dd

end function dp

end