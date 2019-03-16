program main 
 implicit none
 integer i,j
 integer,parameter :: N = 100, LDA = N,LWORK = 100000
 integer INFO 
 double precision :: A(LDA,N),W(N),WORK(LWORK)
 real B(N)
 real t
 t = 200
 do i = 1,N
  do j = 1,N
   A(i,j) = 2*t*delta(i,j) - t*delta(i,j-1) - t*delta(i,j+1)
  end do 
 end do  
 
 !do i = 1,5
  !print*,(A(i,j),j=1,5)
 !enddo
 
 call dsyev('V','U',N,A,LDA,W,WORK,LWORK,INFO)
 
 open (unit =2,file ="matB.txt")
 do i = 0,N-1
  b(i+1) = ((i**2)*((2*3.14)**2))/(8*1*5**2)
  write(2,*) w(i+1),(A(i+1,j=1,N)) 
 end do 
 close(2)
 
 
 contains 
 real function delta(i,j)
 integer i,j
 if (i .eq. j) then
  delta = 1
 else 
  delta = 0
 end if
 end function delta
end program 

