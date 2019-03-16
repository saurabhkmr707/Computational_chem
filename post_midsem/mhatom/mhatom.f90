program main
implicit none
real(kind=8),dimension(10000)::r11,p,R1
 integer::i,j,n
character(len=6) :: str
real(kind=8)::dr=0.0005,e=-0.6,s1,s2,s3,s4,k1,k2,k3,k4
R1(1)=0.000001
p(1)=-1000.0
r11(1)=0.0005

do n=1,20
if(n .gt. 9)then
write(str,"(A4,I2)") "file",n
else
write(str,"(A4,I1)") "file",n
endif
open(unit=2,file=str)


do i=2,10000
  r11(i)=r11(i-1)+dr
 enddo

 do i=1,9999

s1 = dr*f(r11(i),R1(i),P(i),e)
k1 = dr*g(r11(i),R1(i),P(i),e)
s2 = dr*f(r11(i)+dr/2.,R1(i)+k1/2.,P(i)+s1/2.,e)
k2 = dr*g(r11(i)+dr/2.,R1(i)+k1/2.,P(i)+s1/2.,e)
s3 = dr*f(r11(i)+dr/2.,R1(i)+k2/2.,P(i)+s2/2.,e)
k3 = dr*g(r11(i)+dr/2.,R1(i)+k2/2.,P(i)+s2/2.,e)
s4 = dr*f(r11(i)+dr,R1(i)+k3,P(i)+s3,e)
k4 = dr*g(r11(i)+dr,R1(i)+k3,P(i)+s3,e)
P(i+1) = P(i) + (s1+2*s2+2*s3+s4)/6.
R1(i+1) = R1(i) + (k1+2*k2+2*k3+k4)/6.
enddo

do i=1,10000
 write(2,*)r11(i),(abs(r11(i)*R1(i))**2),R1(i),0.0
enddo
e = e+0.01
enddo

contains
real(kind=8) function f(r11,R1,P,e)
implicit none
real(kind=8) :: r11,R1,P,e
f = -(2.*P/r11) - 2*(e+1./r11)*R1
return
end function f


real(kind=8) function g(r11,R1,P,e)
implicit none
real(kind=8) :: r11,R1,P,e
g = P
return
end function g
end


