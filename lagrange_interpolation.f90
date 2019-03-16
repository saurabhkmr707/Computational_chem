module lagrange_interpolation

 contains
 subroutine lagrange(xdata,ydata)
 implicit none
real,dimension(5)::xdata,ydata
integer i,j,k
 real l,x,f,output
 x = 0
 open(unit = 2,file = "lagrange_data_points.txt")
do k = 1,50
x =  k*0.1
 output = 0
 do i = 1,5
  l = 1.0
   do j = 1,5
    if (i == j) then
     l = l*1
    else if(xdata(i)==xdata(j)) then
     l = l*1
    else
     l = l*(x-xdata(j))/(xdata(i)-xdata(j))
    end if
   end do
   output = output + l*ydata(i)
 end do
 write(2,*) x,output
 end do
 close(2)
 end subroutine lagrange
 
end module lagrange_interpolation

program main 

use lagrange_interpolation

implicit none

real,dimension(5)::xdata,fdata
real output
integer i

open (unit = 1,file = "data.txt")
 do i= 1,5
  read(1,*) xdata(i),fdata(i)
 end do
close(1)
call lagrange(xdata,fdata)
end program main
