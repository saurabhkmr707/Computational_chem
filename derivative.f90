program derivative
 
 real x,h1,h2,f_x
 print*,"enter the number you want to caluclate derivative at: "
 read*,x
 h1 = 0.5
 f_x = (f(x+h1)-f(x))/h1
 print*,f_x
 
 contains
 real function f(x)
  real x
  f = -0.1*x**4-0.15*x**3-0.5*x**2-0.25*x+1.2
 end function f
 
end program derivative
