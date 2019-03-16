program main
 implicit none
 real x,g
  x = 0.0
  open (unit = 1, file ="morsepotential.txt")
  do 
   x = x+0.01	
   if (x >= 10) exit
   write(1,*) x,f(x)
  end do 
 close(1)
 
 contains 
 real function f(x)
 real x
 f= 1*(1-exp(-1*(x-0.5)))**2
 
 end function f
 
end program main

