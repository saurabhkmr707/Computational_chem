program main
 implicit none
 complex,dimension(-100:100,0:5000) :: psi,d2_dx2
 !complex,dimension(-100:100) :: d2_dx2
 real,dimension(-100:100,0:5000) :: psir
 real h
 integer i,j
 
 h = 0.02
 
 do i = -100,100
  do j = 0,5000
   psi(i,j) = 0
  end do 
 end do 
 
 do i = -100,100
  psi(i,0) = psi0(i*0.02)
 end do
 
 j =0
 
 d2_dx2(-100,j) = (psi(-98,j) - 2*psi(-99,j) + psi(-100,j))/(0.02**2) 
 d2_dx2(-99,j) = (psi(-97,j) - 2*psi(-98,j) + psi(-99,j))/(0.02**2)
 
 
 do i = -98,100
  d2_dx2(i,j) = d2_dx2(i-1,j)*(h**2) - psi(i-2,j) + 2*psi(i-1,j)
 end do
   
 
 
 do i = -100,100
  psi(i,1) = psi(i,0) + cmplx(0,1)*0.1*((1/(2*14500))*d2_dx2(i,0) - V(i*0.02)*psi(i,0))
 end do 
 
 do j =1,4999
 
 d2_dx2(-100,j) = (psi(-98,j) - 2*psi(-99,j) + psi(-100,j))/(0.02**2) 
 d2_dx2(-99,j) = (psi(-97,j) - 2*psi(-98,j) + psi(-99,j))/(0.02**2)
 
 
 do i = -98,100
  d2_dx2(i,j) = d2_dx2(i-1,j)*(h**2) - psi(i-2,j) + 2*psi(i-1,j)
 end do
 do i = -100,100
  psi(i,j+1) = psi(i,j-1) - (2*cmplx(0,1)*0.1*(-(1/(2*14500))*d2_dx2(i,j) + V(i*0.02)*psi(i,j)))
 end do 
 end do 
 


 
 psir = abs(psi)
 open(unit = 1,file='ha3.txt')
 do i = -100,100
  write(1,*) (psir(i,j), j=0,8), i*0.02
 end do
 close(1)
 
 contains
 complex function psi0(x)
 implicit none
  real x
  psi0 = (((2*20)/3.14)**(1/4))*exp(-(20*(x+0.5)**2))*exp( cmplx(0,1)*20*(x+0.5) )
 end function psi0
 
 real function V(x)
 implicit none
  real x
  if (x < 0) then
   V = 0
  else if (x >= 0) then
   V = 1
  end if
  
  
  
 end function V
 
end program main 
