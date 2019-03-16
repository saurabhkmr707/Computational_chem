module waveFunction
contains
subroutine wave_func(n,y)
 implicit none
 real gp,f,psi
 integer n,y
 gp = 1.0/n
 f = 0.0
 open (unit = 1,file = "psiData.txt")
 do 
  f = f+gp
  if (f .ge. 1) exit
  psi = sqrt(2/1.0)*sin(((y*3.14*f))/1.0)
  write (1,*) f,psi,psi**2
 end do 
 close(1)
 
end subroutine wave_func

end module waveFunction

program main 
use waveFunction
implicit none
integer n,y
print*,"enter n: "
read*,n
print*,"enter quantum number: "
read*,y
call wave_func(n,y)
 
end program 
