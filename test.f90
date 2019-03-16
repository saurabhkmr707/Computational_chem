program main 
 implicit none
 integer i,j
 integer,parameter :: N = 500, LDA = N,LWORK = 100000
 integer INFO 
 double precision :: A(LDA,N),W(N),WORK(LWORK)
 real t,g
 g = 0.05
 t = 200
 do i = 1,N
  do j = 1,N
   A(i,j) = (2*t+g*j*0.01)*delta(i,j) - t*delta(i,j-1) - t*delta(i,j+1)
  end do 
 end do  
 
 !do i = 1,5
  !print*,(A(i,j),j=1,5)
 !enddo
 
 call dsyev('V','U',N,A,LDA,W,WORK,LWORK,INFO)
 
 open (unit =2,file ="theta.txt")
 do i = 1,N
  write(2,*) w(i),A(i,1),A(i,2),a(i,3),i*0.05
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

