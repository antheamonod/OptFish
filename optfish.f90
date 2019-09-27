!Author: Anthea Monod
!Reference: Random Effects Modeling and the Zero-Inflated Poisson Distribution
!Line-Search Optimization without computation of gradient
!To compile, use g95 optfish.f90 -o prog ~/bin/called/*.f90
!To run, use ./optfish > optfish_output to write output in a separate file

module commondata
  save
  integer, parameter :: n=1000 ! sample size
  integer, dimension(n,2) :: y,y_nz,y_z,y_1,y_2
  real, dimension(n,2) :: x,x_nz,x_z,x_1,x_2
  integer :: n_nz,n_z,n_y1,n_y2
  real :: logu!  auxiliary variable to decrease execution time
  real, dimension(n) :: logprod1_y1,logprod2_y2,logprod12_nz ! auxiliary variables to decrease execution time
end module commondata

program test_mle
  use commondata
  use random
  implicit none
  integer :: i,j,k,l,counter
  real :: f,fmax
  real, dimension(4) :: theta0,thetabest,thetatmp,thetatotal,ctm
  real :: alphamax,beta1max,beta2max,umax
  integer, parameter :: partition=20,simus=100
  real :: dummy
  theta0(1)=0.22  ! theta0(1) = alpha, true parameter values
  theta0(2)=0.5  ! theta0(2) = beta1
  theta0(3)=0.5  ! theta0(3) = beta2
  theta0(4)=1.5  ! theta0(4) = u
  alphamax = 0.9
  beta1max = 1. 
  beta2max = 1. 
  umax = 3.

  thetatotal=0.
  write(*,*) 'theta0=',theta0
  write(*,*) 'observations=',n,'simulations=',simus,'partition=',partition
  write(*,*) 'alphamax=',alphamax,'beta1max=',beta1max,'beta2max=',beta2max,'umax=',umax
  write(*,*) '-----------------------we start-----------------------'
  do counter=1,simus
     call gendata(theta0)
     call sort
     fmax=-huge(dummy) ! smallest real value

!start line search loop
     do l = 1,partition
        thetatmp(4) = (umax*l)/partition
        call aux(thetatmp(4)) ! compute expressions that do not depend upon thetatmp(1,2,3)
        do j = 1,partition
           thetatmp(2) = (beta1max*j)/partition
           do k = 1,partition
              thetatmp(3) = (beta2max*k)/partition
              do i=1,partition
                 thetatmp(1) = (alphamax*i)/partition
                 call mletheta(thetatmp,f)
                 if (f.gt.fmax) then
                    fmax = f
                    thetabest = thetatmp
                 end if
              end do
           end do
        end do
     end do
!Write results of current simulation
     write(*,*) 'thetabest=',thetabest, '   iteration=',counter,'/',simus
     thetatotal(1)=thetatotal(1)+thetabest(1)
     thetatotal(2)=thetatotal(2)+thetabest(2)
     thetatotal(3)=thetatotal(3)+thetabest(3)
     thetatotal(4)=thetatotal(4)+thetabest(4)
     ctm(1)=thetatotal(1)/counter
     ctm(2)=thetatotal(2)/counter
     ctm(3)=thetatotal(3)/counter
     ctm(4)=thetatotal(4)/counter
     write(*,*) 'current thetamean=',ctm
  end do
!Write global result (all simulations)
  write(*,*) '--------------------results---------------------------'
  write(*,*) 'theta0=',theta0
  write(*,*) 'observations=',n,'simulations=',simus,'partition=',partition
  write(*,*) 'alphamax=',alphamax,'beta1max=',beta1max,'beta2max=',beta2max,'umax=',umax
  write(*,*) 'thetamean(1)=',thetatotal(1)/simus
  write(*,*) 'thetamean(2)=',thetatotal(2)/simus
  write(*,*) 'thetamean(3)=',thetatotal(3)/simus
  write(*,*) 'thetamean(4)=',thetatotal(4)/simus
  return
end program test_mle

!Objective function
subroutine mletheta(xx,f) ! xx is proposed solution
  use commondata
  implicit none
  integer :: i
  real,dimension(4) :: xx
  real :: f

  f = 0
  do i=1,n_nz
     f = f + 2*log(1-xx(1)) + y_nz(i,1)*x_nz(i,1)*xx(2) +&
          & y_nz(i,2)*x_nz(i,2)*xx(3) + xx(4)*logu - (xx(4)+y_nz(i&
          &,1)+y_nz(i,2))*log(xx(4)+exp(x_nz(i,1)*xx(2))+exp(x_nz(i,2)&
          &*xx(3))) + logprod12_nz(i)
  end do

  do i=1,n_z
     f = f + log( xx(1)*xx(1) + xx(1)*(1-xx(1))*(&
          & exp(log(xx(4)/(xx(4)+exp(x_z(i,1)*xx(2))))*xx(4)) +&
          & exp(log(xx(4)/(xx(4)+exp(x_z(i,2)*xx(3))))*xx(4)) ) + ((1&
          &-xx(1))**2) * exp( log(xx(4)/(xx(4)+exp(x_z(i,1)*xx(2))+exp(x_z(i,2)*xx(3))))*xx(4) ))
  end do

  do i=1,n_y1
     f = f + log(1-xx(1)) + y_1(i,1)*x_1(i,1)*xx(2)  +&
          & logprod1_y1(i) +  log( xx(1)*exp(log( xx(4)&
          &/(exp(x_1(i,1)*xx(2))+xx(4)) ) *xx(4))&
          &/exp(log(exp(x_1(i,1)*xx(2))+xx(4))*y_1(i,1)) + (1&
          &-xx(1))*exp(log(  xx(4)/(exp(x_1(i,1)*xx(2))+exp(x_1(i,2)&
          &*xx(3))+xx(4)))*xx(4))&
          &/exp(log(exp(x_1(i,1)*xx(2))+exp(x_1(i,2)*xx(3))+xx(4))*y_1(i&
          &,1)))
  end do

  do i=1,n_y2
     f = f + log(1-xx(1)) + y_2(i,2)*x_2(i,2)*xx(3)  +&
          & logprod2_y2(i) + log( xx(1)*exp(log( xx(4)&
          &/(exp(x_2(i,2)*xx(3))+xx(4)) ) *xx(4))&
          &/exp(log(exp(x_2(i,2)*xx(3))+xx(4))*y_2(i,2)) + (1&
          &-xx(1))*exp(log(  xx(4)/(exp(x_2(i,1)*xx(2))+exp(x_2(i,2)&
          &*xx(3))+xx(4)))*xx(4))&
          &/exp(log(exp(x_2(i,1)*xx(2))+exp(x_2(i,2)*xx(3))+xx(4))*y_2(i&
          &,2)))
  end do
  return
end subroutine mletheta

!Function to compute (a+b-1)(a+b-2)...a
!If b.le.0, prod returns 1
function prod(a,b)
  implicit none
  real :: a,prod
  integer :: b,i
  prod = 1.
  do i = 0,b-1
     prod = prod*(a+i)
  end do
  return
end function prod

!Auxiliary subroutine
subroutine aux(u)
  use commondata
  implicit none
  real :: u, prod
  integer :: i
  logu = log(u)
  do i=1,n_nz
     logprod12_nz(i) = log(prod(u,y_nz(i,1)+y_nz(i,2)))
  end do
  do i=1,n_y1
     logprod1_y1(i) = log(prod(u,y_1(i,1)))
  end do
  do i=1,n_y2
     logprod2_y2(i) = log(prod(u,y_2(i,2)))
  end do
end subroutine aux

!Generate data, gives a matrix y of values from a Poisson distribution and zeros
subroutine gendata(theta)
  use random
  use commondata
  implicit none
  real,dimension(4),intent(in) :: theta
  double precision :: t
  real :: lam,delta,lam_re
  integer :: i

  y=0  ! by default
  do i=1,n
     x(i,1) = random_normal() ! Generates random number from N(0,1)
     x(i,2) = random_normal() ! Generates random number from N(0,1)
     call random_number(t)    ! Generates random number uniformly on [0,1]
     if (t.gt.theta(1)) then  ! Checks zero-inflation
        lam = exp(x(i,1)*theta(2)) ! Computes Poisson parameter lambda
        delta = random_gamma(theta(4),.TRUE.) ! Generates random number from Gamma distribution with u0 = shape = scale parameters
        lam_re = lam*delta ! Computes random-effects Poisson parameter
        y(i,1) = random_Poisson(lam_re,.TRUE.) ! Generates random number from Poisson distribution with random-effects parameter
     end if
     call random_number(t)
     if (t.gt.theta(1)) then ! Checks zero-inflation
        lam = exp(x(i,2)*theta(3)) ! Computes Poisson parameter lambda
        delta = random_gamma(theta(4),.TRUE.) ! Generates random number from Gamma distribution with u0 = shape = scale parameters
        lam_re = lam*delta ! Computes random-effects Poisson parameter
        y(i,2) = random_Poisson(lam_re,.TRUE.) ! Generates random number from Poisson distribution with random-effects parameter
     end if
  end do
end subroutine gendata

!Subroutine to count and sort the number of pairs of nonzero, zero, mixed y_i1 & y_i2
subroutine sort
  use commondata
  implicit none
  integer :: i

  n_nz = 0
  n_z = 0
  n_y1 = 0
  n_y2 = 0
  do i = 1,n
     if ((y(i,1).ne.0).and.(y(i,2).ne.0)) then
        n_nz = n_nz + 1
        y_nz(n_nz,1) = y(i,1)
        y_nz(n_nz,2) = y(i,2)
        x_nz(n_nz,1) = x(i,1)
        x_nz(n_nz,2) = x(i,2)
     elseif ((y(i,1).eq.0).and.(y(i,2).eq.0)) then
        n_z = n_z + 1
        y_z(n_z,1) = y(i,1)
        y_z(n_z,2) = y(i,2)
        x_z(n_z,1) = x(i,1)
        x_z(n_z,2) = x(i,2)
     elseif ((y(i,1).ne.0).and.(y(i,2).eq.0)) then
        n_y1 = n_y1 + 1
        y_1(n_y1,1) = y(i,1)
        y_1(n_y1,2) = y(i,2)
        x_1(n_y1,1) = x(i,1)
        x_1(n_y1,2) = x(i,2)
     else
        n_y2 = n_y2 + 1
        y_2(n_y2,1) = y(i,1)
        y_2(n_y2,2) = y(i,2)
        x_2(n_y2,1) = x(i,1)
        x_2(n_y2,2) = x(i,2)
     endif
  end do
end subroutine sort
