program main
implicit none
real(kind=8),dimension(201) :: e_x,e_p,t,rk2_x,rk2_p,rk4_x,rk4_p
real(kind=8) :: h,k,m,pi,omega,period,s1,s2,s3,s4,k1,k2,k3,k4
integer :: i



k=1.
m=1.

omega = sqrt(k/m)
pi = 2.*asin(1.)

period = 2*pi/omega
h = 0.02*period

!----EULER--------------------
e_x(1) = 1.
e_p(1) = 0.
t(1) = 0.

open(unit=10,file="euler_p_vs_x.txt")
open(unit=11,file="euler_e_vs_t.txt")
open(unit=12,file="euler_x_vs_t.txt")
do i=1,200
e_p(i+1) = e_p(i) + h*f(t(i),e_x(i),e_p(i))
e_x(i+1) = e_x(i) + h*g(t(i),e_x(i),e_p(i))
t(i+1) = t(i) + h
enddo

do i=1,201
write(10,*) e_x(i),e_p(i)
write(11,*) t(i)/period,(e_x(i)**2 + e_p(i)**2)
write(12,*) t(i)/period,e_x(i)
enddo
!-----------------------------
!----------RK2----------------
rk2_x(1) = 1.
rk2_p(1) = 0.
open(unit=13,file="rk2_p_vs_x.txt")
open(unit=14,file="rk2_e_vs_t.txt")
open(unit=15,file="rk2_x_vs_t.txt")
do i=1,200

s1 = h*f(t(i),rk2_x(i),rk2_p(i))
k1 = h*g(t(i),rk2_x(i),rk2_p(i))
s2 = h*f(t(i)+h,rk2_x(i)+k1,rk2_p(i)+s1)
k2 = h*g(t(i)+h,rk2_x(i)+k1,rk2_p(i)+s1)
rk2_p(i+1) = rk2_p(i) + (s1+s2)/2.
rk2_x(i+1) = rk2_x(i) + (k1+k2)/2.
enddo

do i=1,201
write(13,*) rk2_x(i),rk2_p(i)
write(14,*) t(i)/period,(rk2_x(i)**2 + rk2_p(i)**2)
write(15,*) t(i)/period,rk2_x(i)
enddo

!-----------------------------


!----------RK4----------------
rk4_x(1) = 1.
rk4_p(1) = 0.
open(unit=16,file="rk4_p_vs_x.txt")
open(unit=17,file="rk4_e_vs_t.txt")
open(unit=18,file="rk4_x_vs_t.txt")

do i=1,200
s1 = h*f(t(i),rk4_x(i),rk4_p(i))
k1 = h*g(t(i),rk4_x(i),rk4_p(i))
s2 = h*f(t(i)+h/2.,rk4_x(i)+k1/2.,rk4_p(i)+s1/2.)
k2 = h*g(t(i)+h/2.,rk4_x(i)+k1/2.,rk4_p(i)+s1/2.)
s3 = h*f(t(i)+h/2.,rk4_x(i)+k2/2.,rk4_p(i)+s2/2.)
k3 = h*g(t(i)+h/2.,rk4_x(i)+k2/2.,rk4_p(i)+s2/2.)
s4 = h*f(t(i)+h,rk4_x(i)+k3,rk4_p(i)+s3)
k4 = h*g(t(i)+h,rk4_x(i)+k3,rk4_p(i)+s3)
rk4_p(i+1) = rk4_p(i) + (s1+2*s2+2*s3+s4)/6.
rk4_x(i+1) = rk4_x(i) + (k1+2*k2+2*k3+k4)/6.
enddo

do i=1,201
write(16,*) rk4_x(i),rk4_p(i)
write(17,'(F4.2,a2,F4.2)') t(i)/period," ",(rk4_x(i)**2 + rk4_p(i)**2)
write(18,*) t(i)/period,rk4_x(i)
enddo

!-----------------------------


contains
real(kind=8) function f(t,x,p)
implicit none
real(kind=8) :: t,x,p
f = -x
return
end function f

real(kind=8) function g(t,x,p)
implicit none
real(kind=8) :: t,x,p
g = p
return
end function g

end program main

