program main
 implicit none
 complex(kind = 16),dimension(-100:100,0:5000) :: psi,d2_dx2
 integer i,j
 
 do i =-100,100
  psi(i,0) = (((2.*20.)/3.14)**(1./4.))*exp(cmplx(0.,1.)*20.*(i*0.02+0.5))*exp(-20.*(i*0.02+0.5)**2.)
 end do
 
 psi(-100,1) = psi(-100,0) + cmplx(0.,1.)*0.5*((1./(2.*14500.))* & 
 (((psi(-98,0)-2.*psi(-99,0)+psi(-100,0)))/(0.02*0.02)) - V(-100*0.02)*psi(i,0))

 do i =-99,100
  psi(i,1) = psi(i,0) + cmplx(0.,1.)*0.5*((1./(2.*14500.))* &
   (((psi(i+1,0)-2.*psi(i,0)+psi(i-1,0)))/(0.02*0.02)) - V(i*0.02)*psi(i,0))
 end do 
 
 do j =2,5000
 psi(-100,j) = psi(-100,j-2) - 2.*cmplx(0.,1.)*0.1*(-1/(2.*14500.)* & 
 ((psi(-98,j-1)-2.*psi(-99,j-1)+psi(-100,j-1))/(0.02*0.02)) +V(-100*.02)*psi(-100,j-1))
 
 do i = -99,100
 psi(i,j) = psi(i,j-2) - 2.*cmplx(0.,1.)*0.1*(-1/(2.*14500.)* &
  ((psi(i+1,j-1)-2.*psi(i,j-1)+psi(i-1,j-1))/(0.02*0.02)) +V(i*.02)*psi(i,j-1))
 end do
 
 end do 
 print*,(d2_dx2(-100,100))
 open(unit = 1,file='ha32.txt')
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

