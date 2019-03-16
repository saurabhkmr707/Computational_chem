program main
 implicit none
 real a,b,h,s
 integer i,j,n
 real, dimension(:), allocatable:: x,y 
 print*,'enter the value of a and b: '
 read*,a,b
 print*,"enter n(even no. ): "
 read*, n
 allocate (x(n),y(n))
 x(1) = a
 y(1) = f(x(1))
 h = (b-a)/n
 do i = 2,n
  x(i) = x(i-1) +  h
  y(i) = f(x(i))
 end do
 
 s = f(a)+f(b)
 do i = 2,n,2
  s = s+4*f(x(i))
 end do
 
 do i = 1,n-1,2
  s = s+2*f(x(i))
 end do 
 
 
print*,""
print*,s*h/3
 contains
 real function f(x)
  real x
  f = exp(x)
 end function f
end program main
