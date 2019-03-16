program main
implicit none
real x,v,e,a,b,maxm,minm,pavg,t
real,dimension(1000):: csq
complex,dimension(1000):: psi,test
complex psi0,psi1,iota
integer i,j,k
a = 0.00
b = 0.01
iota = (0.0,1.0)
e = 1.0
do k = 1,1000
open (unit = 2,file = 'prob-trans.txt')
open(unit = 1,file = 'potbar.txt')
x = 0.01
psi(1) = exp (-1*sqrt(2*e)*iota*0)
psi(2) = exp(-1*sqrt(2*e)*iota*0.01)
write(1,*) a,real(psi(1)),aimag(psi(1)),(real(psi(1))**2+aimag(psi(1))**2)
write(1,*) b,real(psi(2)),aimag(psi(2)),(real(psi(2))**2+aimag(psi(2))**2)
do i = 3,1000
 x = x + 0.01
 if ((x .ge. 4) .and. (x .le. 5)) then
  v = 9.0
 else
  v = 0.0
 end if
 psi(i) = 2 * ( 1 + 0.0001 * ( V - E ) ) * psi(i-1) - psi(i-2)
 write(1,*) x,real(psi(i)),aimag(psi(i)),(real(psi(i))**2)+(aimag(psi(i))**2)
 
end do
close(1)

do j=1,1000
 csq(j)=(real(psi(j))**2+aimag(psi(j))**2)
enddo

minm=csq(1)
do j=600,1000
    if (csq(j)<minm) then
        minm=csq(j)
    end if
enddo

maxm=csq(1)
do j=600,1000
    if (csq(j)>maxm) then
        maxm=csq(j)
    end if
enddo

pavg = (minm + maxm)/2
t = 2/(1+pavg)
write (2,*) e,t
e = k*0.1
end do
close(2) 
end program main
