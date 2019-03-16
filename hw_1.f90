
program particle_in_a_box
real, dimension(0:1000):: psix,psiy , V , csq
integer j,i,g
real E,k,min,max,pavg,t,x
open(10,file='data.txt')
E = 0.0
do i=1,1000
k = sqrt(2*E)

do j=0,1000
    if (j.ge.400 .and. j.le.500) then
        V(j) = 9.00
    else
        V(j) = 0.00
    end if
enddo

psix(0) = 1.00
psix(1) = cos(k*0.01)
psiy(0) = 0.00
psiy(1) = -sin(k*0.01)
do j=2,1000
 psix(j) = 2 * ( 1 + 0.0001 * ( V(j) - E ) ) * psix(j-1) - psix(j-2)
 psiy(j) = 2 * ( 1 + 0.0001 * ( V(j) - E ) ) * psiy(j-1) - psiy(j-2)
enddo
open (unit = 1,file = 'random.txt')
 x = 0.0
 do g= 0,1000
 x = x+0.01
 write(1,*) x,psix(g),psiy(g),(psix(g))**2+(psiy(g))**2
 end do 
close(1)
do j=0,1000
 csq(j)=(psix(j)*psix(j)+psiy(j)*psiy(j))
enddo

min=csq(0)
do j=600,1000
    if (csq(j)<min) then
        min=csq(j)
    end if
enddo

max=csq(0)
do j=600,1000
    if (csq(j)>max) then
        max=csq(j)
    end if

enddo

pavg = (min + max)/2
t = 2/(1+pavg)

write(10,*) e,t
E = i*0.1
enddo
end program particle_in_a_box
