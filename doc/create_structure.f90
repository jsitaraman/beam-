!>
!! create the structural dynamic model
!! by forming element matrices and assembling them
!!
subroutine create_structure
!
use rotor
use fem
use parallelvar
!
implicit none
!
real*8 :: mmat(ndof,ndof),kmat(ndof,ndof),cmat(ndof,ndof),qmat(ndof)
real*8 :: theta,thetadot
integer :: i,j,k
!
! set global mass damping and stiffness matrices
! to zero
!
do i=1,ntotaldof
   do j=1,ntotaldof
      mm(i,j)=0.
      kk(i,j)=0.
      if (i==j) then
       dd(i,j)=str_damp  !< diagonal damping for now
      else
       dd(i,j)=0.
      endif
   enddo
   qq(i)=0.
enddo
!
cmat=0.
qmat=0.
!
! need to find theta and thetadot here
! based on rigid motions
!
theta=0.
thetadot=0.
!
do i=elstart,elend
   call m_matrix(i,mmat,theta)
   call k_matrix(i,kmat,theta)
   call q_matrix(i,qmat,theta,thetadot)
   call assemble(mmat,kmat,cmat,qmat,i)
enddo
!
!write(6,*) '===================================================='
!do i=1,7
! write(6,"(7(1X,E14.8))") (mmat(i,j),j=1,7)
!enddo
!write(6,*) '===================================================='
!
!
!call mpi_barrier(mpi_comm_world,ierr)
!
call allreduce_real8(mm,ntotaldof*ntotaldof)
call allreduce_real8(kk,ntotaldof*ntotaldof)
!call allreduce_real8(dd,ntotaldof*ntotaldof)
!
return
end subroutine create_structure
!>
!! Element mass matrix
!!
subroutine m_matrix(ielem,Mmat,theta)
!
use rotor
use fem, only : ndof
use shapefunctions, only : Hs,Hphi
use gauss
!
implicit none
!
integer, intent(in) :: ielem            !< element number
real*8, intent(in)  :: theta            !< twist angle
real*8, intent(out) :: Mmat(ndof,ndof)  !< mass matrix 
!
real*8  :: xmass,wg,ykmsq,xl,FA1,eg1,ct
integer :: i,j,n
!      
xl=dx(ielem)
xmass=xm(ielem)
eg1=eg(ielem)
!eg1=0.
ct=cos(theta)
!
do i=1,4
   do j=i,4
      Mmat(i,j)=0.
      do n=1,6
         wg=weight(n)*0.5*xl
         Mmat(i,j)=Mmat(i,j)+ &
              xmass*wg*Hs(xabs(n),i,xl)*Hs(xabs(n),j,xl)
      enddo
   enddo
enddo
!  
ykmsq=yk2(ielem)+yk1(ielem)
!
do i=1,4
  do j=5,7
    Mmat(i,j)=0.
    do n=1,6
     wg=weight(n)*0.5*xl
     Mmat(i,j)=Mmat(i,j)+wg*xmass*eg1*ct*Hs(xabs(n),i,xl)*Hphi(xabs(n),j-4)
    enddo
  enddo
enddo
!
do i=5,7
   do j=5,7
      Mmat(i,j)=0.
      do n=1,6
         wg=weight(n)*0.5*xl
         Mmat(i,j)=Mmat(i,j)+xmass*ykmsq* &
              wg*Hphi(xabs(n),i-4)*Hphi(xabs(n),j-4)
      enddo
   enddo
enddo
!
do i=1,7
   do j=1,i-1
      Mmat(i,j)=Mmat(j,i)
   enddo
enddo

return
end subroutine m_matrix
!>
!! Element stiffness matrix
!!
subroutine k_matrix(ielem,Kmat,theta)
  use rotor
  use fem, only : ndof
  use shapefunctions, only : Hs,dHs,ddHs,Hphi,dHphi
  use gauss
  !
  implicit none
  !
  ! arguments
  !
  integer , intent(in) :: ielem           !< element number
  real*8, intent(out)  :: kmat(ndof,ndof) !< stiffness matrix
  real*8, intent(in)   :: theta           !< twist of the element
  !
  ! local variables
  !
  integer :: i,j,n
  real*8 :: xi,EI1,EI2,GJ1,xl,xmass,xx
  real*8 :: ct0,ct1,st1,ct2,FA,ykmsq,FA1,term1,term2,wg
  real*8 :: omega2,eg1
  !
  xi=x(ielem)
  EI1=EIy(ielem)
  EI2=EIz(ielem)
  GJ1=GJ(ielem)
  eg1=eg(ielem)
  !eg1=0.
  xl=dx(ielem)
  xmass=XM(ielem)
  xx=xi+xl*0.5
  !
  ct0=cos(theta)
  ct1=ct0*ct0
  st1=sin(theta)**2
  ct2=cos(2*theta)
  omega2=omega*omega
  !
  FA=0.
  do i=ielem+1,nelem
     FA=FA+xm(i)*(x(i)+dx(i)*0.5)*dx(i)
  enddo
  !
  Kmat=0.
  !
  do i=1,4
     do j=i,4
        Kmat(i,j)=0.
        do n=1,6
           wg=weight(n)*0.5*xl
           xx=(1+xabs(n))*0.5*xl
           FA1=xmass*(xl-xx)*(xi+(xl+xx)*0.5)
           term1=(FA+FA1)*dHs(xabs(n),i,xl)*dHs(xabs(n),j,xl)*omega2
           term2=(EI1*ct1+EI2*st1)*ddHs(xabs(n),i,xl)*ddHs(xabs(n),j,xl)
           Kmat(i,j)=Kmat(i,j)+wg*(term1+term2)
        enddo
     enddo
  enddo
  !
  do i=1,4
   do j=5,7
    Kmat(i,j)=0.
    do n=1,6
     wg=weight(n)*0.5*xl
     xx=(1+xabs(n))*0.5*xl
     Kmat(i,j)=Kmat(i,j)+wg*(xmass*(x(ielem)+xx)*omega2*eg1*ct0*&
	                     dHs(xabs(n),i,xl)*Hphi(xabs(n),j-4))
    enddo
   enddo
  enddo
  !
  ykmsq=yk2(ielem)-yk1(ielem)
  !
  do i=5,7
     do j=5,7
        Kmat(i,j)=0.
        do n=1,6
           wg=weight(n)*0.5*xl
           Kmat(i,j)=Kmat(i,j)+xmass*omega2*ykmsq*Hphi(xabs(n),i-4)*&
                Hphi(xabs(n),j-4)*ct2*wg
           Kmat(i,j)=Kmat(i,j)+GJ1*dHphi(xabs(n),i-4,xl)*&
                dHphi(xabs(n),j-4,xl)*wg
        enddo
     enddo
  enddo
  !
  do i=1,7
     do j=1,i-1
        Kmat(i,j)=Kmat(j,i)
     enddo
  enddo
  
  return
end subroutine k_matrix
!>
!! 
!!  Forcing vector
!!
subroutine q_matrix(ielem,Qmat,theta,thetadd)
  use rotor
  use fem, only : ndof
  use shapefunctions, only : Hphi,Hs,dHs
  use gauss
  !
  implicit none
  !
  ! Arguments
  !
  integer, intent(in) :: ielem      !< element number
  real*8, intent(out) :: qmat(ndof) !< forcing vector 
  real*8, intent(in)  :: theta      !< twist
  real*8, intent(in)  :: thetadd    !< torsional acceleration
  !
  real*8 :: xi,ei1,gj1,xl,xmass,ykmsq,wg,omega2,eg1,xx
  real*8 :: st,ct
  integer :: i,n
  !
  xi=x(ielem)
  EI1=EIy(ielem)
  GJ1=GJ(ielem)
  xl=dx(ielem)
  eg1=eg(ielem)
  xmass=XM(ielem)
  ykmsq=yk2(ielem)+yk1(ielem)
  omega2=omega**2
  !
  st=sin(theta)
  ct=cos(theta)
  !
  do i=1,ndof
     Qmat(i)=0
  enddo
  !
  do i=1,4
    do n=1,6
     wg=weight(n)*0.5*xl
     xx=(1+xabs(n))*xl
     Qmat(i)=Qmat(i)&
             -wg*(xmass*omega2*(betap*(x(ielem)+xx))+thetadd*eg1*st)*Hs(xabs(n),i,xl)&
	     -wg*(xmass*omega2*eg1*st*(x(ielem)+xx))*dHs(xabs(n),i,xl)
    enddo
  enddo
  !
  do i=5,7
     do n=1,6
        wg=weight(n)*0.5*xl
        Qmat(i)=Qmat(i)-wg*(xmass*ykmsq*thetadd &
             +xmass*omega2*(yk2(ielem)-yk1(ielem)) &
             *st*ct)*Hphi(xabs(n),i-4)
     enddo
  enddo
  !
  return
end subroutine q_matrix
!>
!! Assemble elemental matrices in to the
!! global stiffness matrix
!!
subroutine assemble(Mmat,Kmat,Cmat,Qmat,ielem)
!
use rotor, only : nelem
use fem 
!
implicit none
!
integer, intent(in) :: ielem          !< element number
real*8, intent(in) :: mmat(ndof,ndof) !< element mass matrix
real*8, intent(in) :: kmat(ndof,ndof) !< element stiffness matrix
real*8, intent(in) :: cmat(ndof,ndof) !< element damping matrix
real*8, intent(in) :: qmat(ndof)      !< element forcing vector
!
integer :: i,j
integer :: ip(ndof)
!
!
! location of the dof's in the global vector
!
call getElementPos(ielem,ip,nelem,ndof)
!
do i=1,ndof
   do j=1,ndof
      mm(ip(i),ip(j))=mm(ip(i),ip(j))+Mmat(i,j)
      kk(ip(i),ip(j))=kk(ip(i),ip(j))+Kmat(i,j)
      dd(ip(i),ip(j))=dd(ip(i),ip(j))+Cmat(i,j)
   enddo
   qq(ip(i))=qq(ip(i))+qmat(i)
enddo
!
return
end subroutine assemble



