module new
 contains 
 real function f(x)
  real x
  f = x**2 + x - 2
 end function f
end module new

program main
 use new
 implicit none
 real a,b,x3
 integer i
 print*,"enter a and b"
 read*,a,b
 if (f(a)*f(b)< 0) then
 x3 = (a+b)/2
  do 	 
   if (f(x3)*f(a)<0)then
    b = x3
    x3 = (a+x3)/2
   else if (f(x3)*f(b)<0) then
    a = x3
    x3 = (b+x3)/2
   end if
   if (f(x3) < 0.001 .and. f(x3)> -0.001) exit
  end do 
 else
 print*,"bisection can't be applied here."
 
 end if
 
 print*, x3, f(x3)
end program main

