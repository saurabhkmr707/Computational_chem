program main
 real,dimension(:),allocatable :: x,y
 real,dimension(:,:),allocatable:: d
 real x_input,l
 integer i,j,n
 
 print*,"enter num of points:(in the newton_data.txt file it is 4) "
 read*,n
 allocate (x(n),y(n),d(n-1,n-1))
 print*,"enter x"
 read*, x_input
 
 open (unit = 1,file = 'newton_data.txt')
 do i = 1,n
  read(1,*) x(i),y(i)
 end do
 close (1)
 
 
 do i = 1,n-1
  d(i,1) = (y(i+1)-y(i))/(x(i+1)-x(i))
 end do 
 
 do j = 2,n-1
  do i = 1,n-1-(j-1)
   d(i,j) = (d(i+1,j-1)-d(i,j-1))/(x(i+j)-x(i))
  end do
 end do
 
 !do i = 1,3
  !print*, (d(i,j),j=1,3)
 !end do 
 
 output = 0
 do i = 2,n
  l = 1.0
  do j = 1,i-1
   l = l*(x_input-x(j))
  end do 
  output = output + l*d(1,i-1)
 end do
 print*,output+y(1)
end program main
