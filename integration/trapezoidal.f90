program main
 implicit none
 integer i,j,n
 real a,b,h,s
 real,dimension(:),allocatable :: x,y
 print*,'enter a and b'
 read*,a,b
 print*,"enter n"
 read*,n
 h =(b-a)/n
 allocate (x(n),y(n))
 x(1) = a
 y(1)= f(x(1))
 
 s = (f(a)+f(b))/2
 do i = 1,n-1
  s = s+f(a+i*h)
 end do 
 
 s = s*h
 print*,s
  
  contains 
  real function f(x)
   real x
   f = exp(x) 
  end function f
end program 
