program main
 implicit none
 real kab,kba,kbc,kcb,kcd,dt,x,t,as1,as2,as3,as4,&
 bs1,bs2,bs3,bs4,cs1,cs2,cs3,cs4,s1,s2,s3,s4
 
 integer i
 real A,B,C,D,A1,B1,C1,D1
 
 kab = 1.
 kba = 3.
 kbc = 4.2
 kcb = 7.3
 kcd = 0.4
 dt = 0.0001
 open(unit = 1,file = 'chkineticsrk4Data.txt')
 
 A = 5.
 B = 0.
 c = 0.
 D = 0.
 t = 0.
 
 do 
          t = t + dt
  	  s1 = dt*fA(t,A,B,C,D)
  	  s2 = dt*fA(t+dt/2.,A+s1/2.,B+s1/2.,C,D)
  	  s3 = dt*fA(t+dt/2.,A+s2/2.,B+s2/2.,C,D)
  	  s4 = dt*fA(t+dt,A+s3,B+s3,C,D)
  	  A1 =  A + (s1+2.*s2+2.*s3+s4)/6.
  	  
  	  s1 = dt*fB(t,A,B,C,D)
  	  s2 = dt*fB(t+dt/2.,A+s1/2.,B+s1/2.,C+s1/2.,D)
  	  s3 = dt*fB(t+dt/2.,A+s2/2.,B+s2/2.,C+s2/2.,D)
  	  s4 = dt*fB(t+dt,A+s3,B+s3,C+s3,D)
  	  B1 =  B + (s1+2.*s2+2.*s3+s4)/6.
  	  
  	  s1 = dt*fC(t,A,B,C,D)
  	  s2 = dt*fC(t+dt/2.,A,B+s1/2.,C+s1/2.,D)
  	  s3 = dt*fC(t+dt/2.,A,B+s2/2.,C+s2/2.,D)
  	  s4 = dt*fC(t+dt,A,B+s3,C+s3,D)
  	  C1 =  C + (s1+2.*s2+2.*s3+s4)/6.
  	  
  	  s1 = dt*fD(t,A,B,C,D)
 	  s2 = dt*fD(t+dt/2.,A,B,C+s1/2.,D)
 	  s3 = dt*fD(t+dt/2.,A,B,C+s2/2.,D)
 	  s4 = dt*fD(t+dt,A,B,C+s3,D)
 	  D1 =  D + (s1+2.*s2+2.*s3+s4)/6.
 	  
	  
	  A = A1
	  B = B1
	  C = C1
	  D = D1
	  write(1,*) t,A,B,C,D
	  if (t .ge. 60) exit
  
 end do 
 
 
 
 contains
 real function fA(t,A,B,C,D)
  implicit none
  real t,A,B,C,D
  fA = -kab*A +kba*B
 end function fA
 
 real function fB(t,A,B,C,D)
  implicit none
  real t,A,B,C,D
  fB = kab*A -kba*B -kbc*B +kcb*C*C
 end function fB
 
 real function fC(t,A,B,C,D)
  implicit none
  real t,A,B,C,D
  fC = 2*kbc*B -2*kcb*C*C - kcd*C
 end function fC
 
 real function fD(t,A,B,C,D)
  implicit none
  real t,A,B,C,D
  fD = kcd*C
 end function fD
 
 
 
end program main
