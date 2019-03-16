!question - use rk4 to esitmate y(0.4) if yprime = x**2+y**2 with y(0) = 0 let h = 0.2

program main
 implicit none
 real yk1,yk,h,s1,s2,s3,s4,xk
 integer i,j
 
 h = 0.2
 yk = 0.
 xk = 0.
 do 
 	s1 = h*f(xk,yk)
 	s2 = h*f(xk+h/2.,yk+s1/2.)
 	s3 = h*f(xk+h/2.,yk+s2/2.)
 	s4 = h*f(xk+h,yk+s3)
 	yk1 = yk+(s1+2*s2+2*s3+s4)/6.
 	yk = yk1
 	xk = xk + h
 	if(xk .eq. 0.4) exit
 	
 enddo
 
 print*, yk1
 
 contains 
 real function f(x,y)
 implicit none
 real x,y
 f = X**2+y**2
 end function f
end program main
