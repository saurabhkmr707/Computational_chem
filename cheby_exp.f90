program cheby_pol

 implicit none
 integer i,j,n
 real,dimension(:),allocatable:: x,c,t,g
 real pi,xx,f
 pi = 3.14
 
 print*,"enter degree(n): "
 read*,n
 allocate(x(n+1),c(n+1),t(n+1),g(n+1))
 
 do j = 1,n+1
  x(j) = cos(((2*(j-1)+1)*pi)/(2*(n+1)))
 end do 
 
 do i = 1,n+1
  c(i) = 0
  do j = 1,n+1
   if (i == 1) then
    c(i) = c(i) + (exp(x(j)))/(n+1)
   else
    c(i) = c(i) + 2*(exp(x(j))*(cos((i-1)*acos(x(j)))))/(n+1)
   end if
  end do 
  
 end do 
 
 open (unit = 1,file = 'chebyshev_exp.txt')
 xx = -1.0
 do 
 xx = xx +0.05
 f = 0.0
  do i = 1,n+1
   f = f + c(i)*(cos((i-1)*acos(xx)))
  end do 
  write(1,*) xx,f
  if (xx .ge. 1.0)exit
  end do 
 close (1)
 
 
 print*,(x(j),j=1,n+1)
 print*,''
 print*,(c(i),i=1,n+1)
 print*,''
 !print*,(f(i),i=1,n+1)
 
end program cheby_pol
