program main
 implicit none
 complex(kind = 16),dimension(-100:100,0:5000) :: psi,d2_dx2
 integer i,j
 
 do i =-100,100
  psi(i,0) = (((2.*20.)/3.14)**(1./4.))*exp(cmplx(0.,1.)*20.*(i*0.02+0.5))*exp(-20.*(i*0.02+0.5)**2.)
 end do
 
 d2_dx2(-100,0) = (psi(-98,0) - 2.*psi(-99,0) + psi(-100,0))/(0.02**2.) 
 d2_dx2(-99,0) = (psi(-97,0) - 2.*psi(-98,0) + psi(-99,0))/(0.02**2.)
 
 do i = -98,100
  d2_dx2(i,0) = d2_dx2(i-1,0)**(0.02**2.) + 2.*psi(i-1,0) - psi(i-2,0)
 end do
 
 
 do i =-100,100
  psi(i,1) = psi(i,0) + cmplx(0.,1.)*0.5* &
   ((1./(2.*14500.))*d2_dx2(i,0) - V(i*0.02)*psi(i,0))
 end do 
 
 do j =1,4999
 
 d2_dx2(-100,j) = (psi(-98,j) - 2.*psi(-99,j) + psi(-100,j))/(0.02**2.) 
 d2_dx2(-99,j) = (psi(-97,j) - 2.*psi(-98,j) + psi(-99,j))/(0.02**2.)
 
 
 do i = -98,100
  d2_dx2(i,j) = d2_dx2(i-1,j)*(0.02**2.) - psi(i-2,j) + 2.*psi(i-1,j)
 end do
 do i = -100,100
  psi(i,j+1) = psi(i,j-1) - (2.*0.1*cmplx(0.,1.)*& 
  (-(1./(2.*14500.))*d2_dx2(i,j)*1.0+1.0*V(i*0.02)*psi(i,j)))
 end do 
 end do 
 print*,(d2_dx2(-100,100))
 open(unit = 1,file='ha31.txt')
 do i = -100,100
  write(1,*) ((real(psi(i,j)))**2+(aimag(psi(i,j)))**2, j=0,5000), i*0.02
 end do
 close(1)
 
 contains
 real(kind = 16) function V(x)
 implicit none
  real(kind = 4) x
  if (x < 0 ) then 
   v = 0.
  else 
   v = 1.
  endif
 end function V
end program main

