program main
 implicit none
 real,dimension(10000) :: r,r11,p
 real,dimension(20) :: e0
 real cR
 character(len = 100) :: str
 integer i,j,e
 real dr,s1,s2,s3,s4,k1,k2,k3,k4
 
 
 dr = 0.0005
 r11(1) = 0.000001
 p(1) = -1000 
 r(1) = 0.0005
 !r-grid
 do i =2,10000
 	r(i) = r(i-1) + dr
 end do
 
 !e-grid
 e0(1) = -0.6
 do i = 2,20
 	e0(i) = e0(i-1) + 0.01 
 end do
 
 do e = 1,20
 
 	if(e .gt. 9)then
	write(str,"(A4,I2)") "file",e
	else
	write(str,"(A4,I1)") "file",e
	endif
	
	open(unit = 2,file = str)
	
 	do i =1,9999
 		s1 = dr*f(r(i),p(i),r11(i),e0(e))
 		k1 = dr*g(p(i))
 		s2 = dr*f(r(i)+dr/2.,p(i)+s1/2.,r11(i)+k1/2.,e0(e))
 		k2 = dr*g(p(i)+s1/2.)
 		s3 = dr*f(r(i)+dr/2.,p(i)+s2/2.,r11(i)+k2/2.,e0(e))
 		k3 = dr*g(p(i)+s2/2.)
 		s4 = dr*f(r(i)+dr,p(i)+s3,r11(i)+k3,e0(e))
 		k4 = dr*g(p(i)+s3)
 		
 		p(i+1) = p(i) + (s1+2.*s2+2.*s3+s4)/6.
 		r11(i+1) = r11(i) + (k1+2.*k2+2.*k3+k4)/6.
 		write(2,*) r(i),(abs(r(i)*r11(i)))**2.
 	end do
 	
 	 
 end do
 
 contains 
 real function f(r,p,r11,e)
  implicit none
  real r,p,r11,e
  f = -2.*p/r -2.*(e+1./r)*r11
 end function f
 
 real function g(p)
  implicit none
  real p
  g = p
  
 end function g
end program main
