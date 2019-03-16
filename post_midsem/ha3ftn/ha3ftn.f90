program TDSE

implicit none

integer n,i,j,a,tnum,time
character (len=100):: haf
real pi,delt,delx,x0,p0,m,xmin,k0
real,dimension(256) :: V,k,ksq,x
real,dimension(256)::psireal,psiimag
complex,dimension(256)::psiinit,psiatone,psitemp,psivary
complex io


io=(0.0,1.0)
delt=0.1
delx=0.02
n=256
a=20
x0=-0.5
p0=20.0
m=14500.0
xmin=-2.0
tnum=5000
pi=3.1415
k0=p0

open(unit=10, file='pot')
																									

call xgrid(x,xmin,delx,n)

call pot(x,n,V)
																									!
								
do i=1,n
	write(10,*)x(i),V(i)
enddo



call kval(n,pi,k,ksq,delx)




open(unit=11, file='psi-0')



do i=1,n

	psiinit(i) = ((2.0*a/pi)**0.25) * exp(io*k0*(x(i)-x0)) * exp((-a)*(x(i)-x0)*(x(i)-x0))
	

enddo




do i=1,n
	write(11,*)x(i),realpart(psiinit(i))**2+imagpart(psiinit(i))**2
	psitemp(i)=psiinit(i)
enddo


call furfun(psitemp,x,k,ksq,n)


do i = 1,n
		psiatone(i) = psiinit(i) + io*delt*( psitemp(i) - (V(i)*psiinit(i)))
enddo

open(unit=12, file='1')

do i=1,n
		write(12,*)x(i),realpart(psiatone(i))**2+imagpart(psiatone(i))**2
		psitemp(i)=psiatone(i)
enddo




do time= 2 , tnum
        write(haf,19)time +10000
        19 format ('',i5)
	open(unit=15, file=haf)
	call furfun(psitemp,x,k,ksq,n)
	
	do i=1,n
			psivary(i) = psiinit(i) +  (2.0*io*delt)*(psitemp(i) - (V(i)*psiatone(i)))
	
	enddo

	do i=1,n
				write(15,*)x(i),realpart(psivary(i))**2+imagpart(psivary(i))**2
	enddo

	do i=1,n
				psiinit(i) = psiatone(i)
				psitemp(i) = psivary(i)
				psiatone(i) = psivary(i)
	
	enddo



enddo



end






subroutine furfun(psitemp,x,k,ksq,n)
integer i,j,n
real pi,m
real ,dimension(256) ::x,k,ksq
complex, dimension(256) :: fi,psitemp
complex io,q,r


pi=3.1415
io=(0.0,1.0)
m=14500.0

do i=1,n
		q=(0.0,0.0)
		
		do j=1,256
				q=q+psitemp(j)*exp((-io)*k(i)*x(j))
		enddo

		fi(i)=(q/sqrt(256.0))
enddo



do i=1,n

	fi(i) = (-ksq(i)) * fi(i)

enddo


do i=1,n
		r=(0.0,0.0)
		
		do j=1,256
				r = r + fi(j)*exp((io)*k(j)*x(i))
		enddo

		psitemp(i)=(r/sqrt(256.0))
enddo

do i=1,n

	psitemp(i)=psitemp(i)/(2.0*m)

enddo
return

end


subroutine xgrid(x,xmin,delx,n)
integer i,n
real, dimension(256) :: x
real delx,xmin

do i=1,n

		x(i)=xmin + (i-1)*delx
enddo

return
end



subroutine pot(x,n,V)
integer i,n
real,dimension(256)::V,x


do i=1,n
	if (x(i) < 0) then
		V(i)=0.0
	else
		V(i)=1.0
	endif

enddo

return

end




subroutine kval(n,pi,k,ksq,delx)
integer i,n
real pi,delx
real,dimension(256)::k,ksq



do i=1,n
	if(i < (n/2.0) ) then
		k(i)=(2.0*pi*(i-1))/(n*delx)
		ksq(i)=k(i)**2
		
	else
		k(i)=(2.0*pi*(i-1-n))/(n*delx)
		ksq(i)=k(i)**2
	endif

enddo

return
end
