

program hydrogenatom
implicit none
double precision s11,s12,s13,s14,s21,s22,s23,s24
double precision r,Rr,h,p,E,p1,Rr1
integer i,j
character*12 file1,file2
h=(5.0-0.0005)/10000.0

E=-0.6

do j=1,20
Rr=0.000001
p=-1000.0
r=0.0005

write(file1,'("Rr",I2)')1+j
write(file2,'("rdist",I2)')2+j
open(1+j, file=file1)
open(2+j, file=file2)

write(1+j,*) Rr,r
write(2+j,*) ((r*Rr)**2.0),r

do i=1,10000
s11=h*f1(p)
s21=h*f2(p,r,Rr,E)

s12=h*f1(p+s11/2.0)
s22=h*f2(p+(s21/2.0),r+(h/2.0),Rr,E)

s13=h*f1(p+s12/2.0)
s23=h*f2(p+(s22/2.0),r+(h/2.0),Rr,E)

s14=h*f1(p+s13)
s24=h*f2(p+s23,r+h,Rr,E)

Rr1=Rr+((s11+(2.0*s12)+(2.0*s13)+s14)/6.0)
p1=p+((s21+(2.0*s22)+(2.0*s23)+s24)/6.0)

write(1+j,*) Rr1,r+h
write(2+j,*) ((r*Rr)**2.0),r+h

r=r+h
Rr=Rr1
p=p1
enddo

Close(1+j)
Close(2+j)

E=E+0.01
enddo
contains

double precision function f2(p,r,Rr,E)
implicit none
double precision p,r,Rr,E
f2=-(2.0/r)*p-(2.0*(E+(1.0/r))*Rr)
return
end function

double precision function f1(p)
implicit none
double precision p
f1=p
return
end function

end program hydrogenatom


