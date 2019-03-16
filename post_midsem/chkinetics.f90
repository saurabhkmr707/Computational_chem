program main
 implicit none
 real(kind = 8) kab,kba,kbc,kcb,kcd,dt,x,t
 integer i
 double precision,dimension(600000) :: A,B,C,D
 
 kab = 1
 kba = 3
 kbc = 4.2
 kcb = 7.3
 kcd = 0.4
 dt = 0.0001
 x = 0.0001 
 open(unit = 1,file = 'chkineticsData.txt')
 
 A(1) = 5
 B(1) = 0
 c(1) = 0
 D(1) = 0
 
 write(1,*) A(1),B(1),C(1),D(1),x
 
 do i = 1,599999
  A(i+1) = A(i) + dt*(-kab*A(i)+kba*B(i))
  B(i+1) = B(i) + dt*(kab*A(i)-kba*B(i)-kbc*B(i)+kcb*(C(i)*C(i)))
  C(i+1) = C(i) + dt*(2*kbc*B(i) - 2*kcb*(C(i)*C(i)) -kcd*C(i))
  D(i+1) = D(i) + dt*(kcd*C(i))
  t = (i+1)*0.0001
  write(1,*) A(i+1),B(i+1),C(i+1),D(i+1),t
  
 end do 
 
 close(1)
 
end program main
