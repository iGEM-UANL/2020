program Draining_Flow
  !THIS CODE IS MADE FOR THE SIMULATION OF THE DRAINING FLOW
implicit none
!This program is to simulate the draining flow
!First we are going to dreclare the VARIABLES and arrays
real*8 EPSILON,B1,C,T_C,H_MIN,L,A
real*8 H0,RHO,SIGMA,MIU,D,G0,A_G,D_G
real*8 CTE1,CTE2,CTE3
real*8 dx,dt,t_inicial,t_final,coeffa,coeffb,coeffc,coeffd
real*8, dimension (501,501)::H
real*8, dimension (501,501)::G
real*8, dimension (501,501)::U
real*8, dimension(501)::a1,b,c1,d1,XU
real*8 Q,HXXX,GXX,GX,PIX,UX,HX,HXX,HXXXX
real*8 BETA,CT1,CT2,gravity, B0
integer X,T,m,n,k,nsteps
!Seting the files
open(10, File="H(X,T)VSX.txt")
open(100, file="XVSG(X,T).txt")
open(200, file="XVSU(X,T).txt")
open(11, File="H.plt")
open(12, file="G.plt")
open(13, file="U.plt")
write(11,*) "plot  'H(X,T)VSX.txt' w l "
write(12,*) "plot 'XVSG(X,T).txt' w l "
write(13,*) "plot 'XVSU(X,T).txt' w l "


!INPUT PARAMETERS OF THE SYSTEM (IN CGS SYSTEM)
!ALL OF THESE PARAMETERS ARE DIMENSIONLESS
!Here is the section to use the program, for any situation the only thing to change are the next PARAMETERS
!---------------------------------------------
!---------------------------------------------
EPSILON=
C= !(3/4*EPSILON**2)
T_C= !CHARACTERISTIC TIME
H_MIN= !H ASTERISK IN THE PROBLEM
gravity=
L=!Length in cm
H0= !INITIAL THICKNESS
RHO= !DENSITY OF THE LIQUID WITHOUT SURFACTANT
SIGMA= !Surfacte tension
MIU= !VISCOSITY
D= !Diffusion constant
G0= !Initial surfactante concentration
!Another constants depending the foam
A=
A_G=(A/(EPSILON**2))
D_G=(D/((L*L)/T_C))
B0=(RHO*gravity*(L*L))/(SIGMA)
B1=B0/EPSILON
!After this is the same code for any situation
!---------------------------------------------------
!Mesh and discretization
!We need to set the mesh
t_inicial=0.0  !INITIAL TIME
t_final=  !Final time
nsteps=  !NUMBER OF STEPS IN TIME
dt=(t_final-t_inicial)/(nsteps) !INFINITISIMAL STEP
K=300
dx= !INFINITISIMAL STEP OF X
m= !This values are on the model
n=
BETA=1
!---------------------------------------------
!-----------------------------------------------
!INITIAL CONDITION AND BOUNDARY IN ARRAYS
!Those are declared on the model
do X=1,k
H(X,1)=H0
end do

do X=1,k
G(X,1)=G0
end do
!-----------------------------------------------
!FINITE DIFFERENCES METHOD

do T=1,nsteps

do X=1,K   !SOLVING U(X,T)

a1(X)=coeffa(X,T,H,dx,K)
b(X)=coeffb(X,T,H,dx,K,)    !THOMAS COEFFICENTS
c1(X)=coeffc(X,T,H,dx,K,)
d1(X)=coeffd(X,T,C,G,H_MIN,B1,dx,A_G,m,n,K,H)
end do


call solve_tridiag(a1,b,c1,d1,XU,k,T) !THOMAS ALGORITHM

do X=1,K
    U(X,T)=XU(X)
end do


    do X=1,K
    H(X,T+1)=(dt/dx)*((-Q(X+1,T,B1,dx,H_MIN,m,n,H,k,A_G,U,H0)+Q(X-1,T,B1,dx,H_MIN,m,n,H,k,A_G,U,H0)))+H(X,T)
    end do

do X=1,K
IF(X.EQ.1)THEN
GX=0
GXX=(2*G(X+1,T)-2*G(X,T))/(dx**2)
UX=(U(X+1,T)-0)/(2*dx)
G(X,T+1)=-BETA*(U(X,T)*GX)-BETA*(UX*G(X,T))+D_G*(dx/dx)*GXX+G(X,T)
ELSE
IF(X.EQ.K)THEN
GX=0
GXX=(-2*G(X,T)+2*G(X-1,T))/(dx**2)
UX=(0-U(X-1,T))/(2*dx)
G(X,T+1)=-BETA*(U(X,T)*GX)-BETA*(UX*G(X,T))+D_G*(dx/dx)*GXX+G(X,T)
ELSE
GX=(G(X+1,T)-G(X-1,T))/(2*dx)
GXX=(G(X+1,T)-2*G(X,T)+G(X-1,T))/(dx**2)
UX=(U(X+1,T)-U(X-1,T))/(2*dx)
G(X,T+1)=-BETA*(U(X,T)*GX)-BETA*(UX*G(X,T))+D_G*(dx/dx)*GXX+G(X,T)
END IF
END IF
end do
end do
do T=1,nsteps
do X=1,k
write(10,*)  H(X,T), dx*(X-1)
write(100,*) G(X,T), dx*(X-1)
write(200,*) U(X,T), dx*(X-1)
end do
end do
!Calling gnuplot for make the graphs
call system('gnuplot -p H.plt')
call system('gnuplot -p G.plt')
call system('gnuplot -p U.plt')

end program !End of the code
!Here is the function's section
!Here we calculate HXXX,Q,PI AND CONSTANTS
function Q(X,T,B1,dx,H_MIN,m,n,H,k,A_G,U,H0)
implicit none
real*8 dx,H_MIN,B1,Q,HXXX,PIX,A_G,H0,CTE1,CTE2,CTE3,dt
integer X,T,m,n,k
real*8,dimension (201,201)::H
real*8,dimension (201,201)::U
if (X.eq.0 .or. X.EQ.K+1) then  !CONDITIONS THAT WE USE IN THE BOUNDARY CONDITIONS
Q=0
else
if (X.EQ.1) then
HXXX=(H(X+2,T)-2*H(X+1,T)+2*1-(1-dx))/(2*dx**3)
CTE1=((H(X+1,T)-1)/(H(X,T)*dx))
CTE2=(m*((H_MIN/H(X,T))**m))
CTE3=(n*((H_MIN/H(X,T))**n))
PIX=CTE1*(CTE2-CTE3)
Q=U(X,T)*H(X,T)+(H(X,T)**3)*(B1+HXXX+A_G*PIX)
else
	if (X.EQ.2) then
    HXXX=((H(X+2,T)-2*H(X+1,T)+2*1-(1-dx)))/(2*dx**3)
	CTE1=((H(X+1,T)-H(X-1,T))/(H(X,T)*dx))
    CTE2=(m*((H_MIN/H(X,T))**m))
    CTE3=(n*((H_MIN/H(X,T))**n))
    PIX=CTE1*(CTE2-CTE3)
    Q=U(X,T)*H(X,T)+(H(X,T)**3)*(B1+HXXX+A_G*PIX)
else
	if(X.EQ.K-1) then
    HXXX=((1-2*H(X+1,T)+2*H(X-1,T)-H(X-2,T)))/(2*dx**3)
	CTE1=((H(X+1,T)-H(X-1,T))/(H(X,T)*dx))
    CTE2=(m*((H_MIN/H(X,T))**m))
    CTE3=(n*((H_MIN/H(X,T))**n))
    PIX=CTE1*(CTE2-CTE3)
	Q=U(X,T)*H(X,T)+(H(X,T)**3)*(B1+HXXX+A_G*PIX)
else
    if(X.EQ.K) then
	HXXX=(((1+dx)-2*1+2*H(X-1,T)-H(X-2,T)))/(2*dx**3)
	CTE1=((1-H(X-1,T))/(H(X,T)*dx))
    CTE2=(m*((H_MIN/H(X,T))**m))
    CTE3=(n*((H_MIN/H(X,T))**n))
    PIX=CTE1*(CTE2-CTE3)
    Q=U(X,T)*H(X,T)+(H(X,T)**3)*(B1+HXXX+A_G*PIX)
else
    HXXX=((H(X+2,T)-2*H(X+1,T)+2*H(X-1,T)-H(X-2,T)))/(2*dx**3)
    CTE1=((H(X+1,T)-H(X-1,T))/(H(X,T)*dx))
    CTE2=(m*((H_MIN/H(X,T))**m))
    CTE3=(n*((H_MIN/H(X,T))**n))
    PIX=CTE1*(CTE2-CTE3)
    Q=U(X,T)*H(X,T)+(H(X,T)**3)*(B1+HXXX+A_G*PIX)
end if
end if
end if
end if
end if
return
end


function coeffa(X,T,H,dx,K)
implicit none
real*8,dimension(201,201)::H
real*8 dx,coeffa
integer X,T,K
if (X.EQ.1) THEN
coeffa=H(X,T)-((H(X+1,T)-1)/(2))
ELSE
IF (X.EQ.K) THEN
coeffa=H(X,T)-((1-H(X-1,T))/(2))
ELSE
coeffa=H(X,T)-((H(X+1,T)-H(X-1,T))/(2))
END IF
END IF
return
end function

function coeffb(X,T,H,dx,K)
implicit none
real*8,dimension(201,201)::H
real*8 dx,coeffb
integer X,T,K
coeffb=-2*H(X,T)
return
end function

function coeffc(X,T,H,dx,K)
implicit none
real*8,dimension(201,201)::H
real*8 dx,coeffc
integer X,T,K
IF (X.EQ.1) THEN
coeffc=H(X,T)+((H(X+1,T)-1)/(2))
ELSE
IF (X.EQ.K) THEN
coeffc=H(X,T)+((1-H(X-1,T))/(2))
ELSE
coeffc=H(X,T)+((H(X+1,T)-H(X-1,T))/(2))
END IF
END IF
return
end function

function coeffd(X,T,C,G,H_MIN,B1,dx,A_G,m,n,K,H)
implicit none
real*8,dimension(K+1, K+1)::G
real*8,dimension(K+1,K+1)::H
real*8 coeffd
real*8 C,H_MIN,B1,dx,GX,HXXX,HX,CT1,CT2,PIX,A_G
integer X,T,K,m,n
IF(X.EQ.1) THEN
GX=0
HXXX=(H(X+2,T)-H(X+1,T)+2*1-(1-dx))/(2*dx**3)
HX=(H(X+1,T)-1)/(2*dx)
CT1=(H_MIN/H(X,T))**m
CT2=(H_MIN/H(X,T))**n
PIX=(HX/H(X,T))*(m*CT1-n*CT2)
coeffd=(2*dx**2)*C*(GX-H(X,T)*(B1+HXXX+A_G*PIX))
ELSE
IF(X.EQ.K) THEN
GX=0
HXXX=((1+dx)-1+2*H(X-1,T)-H(X-2,T))/(2*dx**3)
HX=(1-H(X-1,T))/(2*dx)
CT1=(H_MIN/H(X,T))**m
CT2=(H_MIN/H(X,T))**n
PIX=(HX/H(X,T))*(m*CT1-n*CT2)
coeffd=(2*dx**2)*C*(GX-H(X,T)*(B1+HXXX+A_G*PIX))
ELSE
IF(X.EQ.K-1)THEN
GX=(G(X+1,T)-G(X-1,T))/(2*dx)
HXXX=(1-H(X+1,T)+2*H(X-1,T)-H(X-2,T))/(2*dx**3)
HX=(H(X+1,T)-H(X-1,T))/(2*dx)
CT1=(H_MIN/H(X,T))**m
CT2=(H_MIN/H(X,T))**n
PIX=(HX/H(X,T))*(m*CT1-n*CT2)
coeffd=(2*dx**2)*C*(GX-H(X,T)*(B1+HXXX+A_G*PIX))
ELSE
GX=(G(X+1,T)-G(X-1,T))/(2*dx)
HXXX=(H(X+2,T)-H(X+1,T)+2*H(X-1,T)-H(X-2,T))/(2*dx**3)
HX=(H(X+1,T)-H(X-1,T))/(2*dx)
CT1=(H_MIN/H(X,T))**m
CT2=(H_MIN/H(X,T))**n
PIX=(HX/H(X,T))*(m*CT1-n*CT2)
coeffd=(2*dx**2)*C*(GX-H(X,T)*(B1+HXXX+A_G*PIX))
END IF
END IF
END IF
return
end function
!Solving the tridiagonal matrix
subroutine solve_tridiag(a1,b,c1,d1,XU,k,T)
implicit none
integer,intent(in) :: K
real(8),dimension(K+1),intent(in) :: a1,b,c1,d1
real(8),dimension(K+1),intent(out) :: XU
real(8),dimension(K+1) :: cp,dp
real(8) :: m
integer i,T


! initialize c-prime and d-prime
cp(1) = c1(1)/b(1)
dp(1) = d1(1)/b(1)
! solve for vectors c-prime and d-prime
do i = 2,K
m = b(i)-cp(i-1)*a1(i)
cp(i) = c1(i)/m
dp(i) = (d1(i)-dp(i-1)*a1(i))/m
enddo
! initialize x
XU(K) = dp(K)
!*, U(K,T),K
! solve for x from the vectors c-prime and d-prime
do i = K-1, 1, -1
XU(i) = dp(i)-cp(i)*XU(i+1)
end do

end subroutine solve_tridiag
