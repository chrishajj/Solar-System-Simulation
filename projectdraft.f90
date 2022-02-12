module variables
  implicit none
  integer, parameter :: N=9
  real(16), parameter :: Ttotal=30,Msun=1.989D30, Mearth=5.972D24, Mmoon=7.35D22, Dmoon=384.4D6, G=6.67D-11
  real(16), parameter ::  pi=3.1415 , wearth=sqrt(G*Mearth/Dmoon**3), period=(2*pi/wearth)
  real(16), parameter :: D1=147.84D9, M1= Mearth*Msun, w1= sqrt(G*Msun/D1**3), period1=(2*pi/w1), D2=D1+Dmoon
  real(16), parameter :: Mmars= 6.39D23, Dmars=236.56D9, Wmars=sqrt(G*Msun/Dmars**3), Tmars=(2*pi/Wmars)
  real(16), parameter :: Mmercury=0.0553*Mearth, Mvenus=0.815*Mearth, Mjupiter= 317.8284*Mearth
  real (16), parameter:: Muranus= 14.5357*Mearth, Mneptune= 17.1478*Mearth, Msaturn=95.1609*Mearth
  real (16), parameter :: Dmercury=0.3871*D1, Dvenus=0.7233*D1, Djupiter=D1*5.2028, Dsaturn=D1*9.5428
  real (16), parameter:: Duranus=D1*19.1921, Dneptune= D1*30.0689
  real(16), dimension (N,N) :: rij
  real(16), dimension (N) :: M,W,D
end module variables

program Nbody
  use variables
  implicit none
  integer ::smthg, imax,i
  real (16), dimension(4*N) :: X, dX
  real(16):: dt, t
  external :: deriv
  open (1,file='orbits.txt')

  M(1)=Msun
  M(2)=Mmercury
  M(3)=Mvenus
  M(4)=Mearth
  M(5)=Mmars
  M(6)=Mjupiter
  M(7)=Msaturn
  M(8)=Muranus
  M(9)=Mneptune
 ! M(10)=Mphobos
  
  D(1)=0
  D(2)=Dmercury
  D(3)=Dvenus
  D(4)=D1
  D(5)=Dmars
  D(6)=Djupiter
  D(7)=Dsaturn
  D(8)=Duranus
  D(9)=Dneptune
   
  W(1)=0
  
  do i=2,N
     W(i)=sqrt(G*Msun/D(i)**3)
  end do
  
  t=0
  dt=0.0001*2*pi/W(8)
  
  imax=100*period

  ! X
  
  do i=1,N
     X(i)=D(i)
  end do

  ! Y
  
  do i=N+1, 2*N
     X(i)=0.
  end do
  
  ! Vx
  do i=2*N+1, 3*N
     X(i)=0.
  end do
  
  ! Vy
  do i=1,N
     smthg=3*N+i
     X(smthg)=D(i)*W(i)
  end do

  !Vy Phobos
 ! X(4*N)=Dphobos*W(10) + Dmars*W(5)
  ! RK4 LOOP
  
  do smthg=1,1
     write (1,*) X
   !  write (1,*)
    ! write (1,*)
     call rk4(t,X,dt,4*N,deriv)
     t=t+dt
  end do
end program Nbody


! THE FINAL DERIV
subroutine deriv(dt, X, dX,Np)
  use variables
  implicit none
  integer(16) ::Np, i,j, cte, xx, yy
  real (16), dimension (4*N) :: X, dX
  real (16) :: dt

  dX=0
  
 ! MOVE VELOCITIES FROM X to dX
  do i=1,2*N
     cte=2*N + i
     dX(i)= X(cte)
  end do

 ! DEFINE Rij
  do i=1,N
     do j=1,N
        rij(i,j)= sqrt((X(i)-X(j))**2+(X(N+i)-X(N+j))**2)
     end do
  end do
  
 ! MODIFY dX 
  do i=1,N
     do j=1,N
        xx=2*N+i
        yy=3*N+i
        if (i.ne.j)then
           dX(xx)=dX(xx)-G*M(j)*(X(i)-X(j))/(rij(i,j))**3
           dX(yy)=dX(yy)-G*M(j)*(X(N+i)-X(N+j))/(rij(i,j))**3
        endif
     end do
  end do
end subroutine deriv

subroutine rk4(t,x,dt,n,deriv)
! 4th order Runge-Kutta. See Numerical Recipes p. 701ff
implicit none
integer           , intent(in)    :: n
real(16)              , intent(in)    :: t, dt
real(16), dimension(n), intent(inout) :: x
real(16)                              :: ddt
real(16), dimension(n)                :: xp, k1, k2, k3, k4

ddt = 0.5*dt
call deriv(t,x,k1,n)      ; xp = x + ddt*k1
call deriv(t+ddt,xp,k2,n) ; xp = x + ddt*k2
call deriv(t+ddt,xp,k3,n) ; xp = x +  dt*k3
call deriv(t+dt,xp,k4,n)  ; x  = x +  dt*( k1 + 2.0*k2 + 2.0*k3 + k4 )/6.0
end subroutine rk4
