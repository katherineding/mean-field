  program main
  implicit none
  integer, parameter :: nx=16
  integer, parameter :: ny=8
  integer, parameter :: nsites=nx*ny
  integer, parameter :: ndim=2*nsites
  integer, parameter :: cmax=40
  integer, parameter :: cmin=20
  integer, parameter :: q=16
  integer i,j,count,site,flag
  integer info
  double precision hubu,sumd,dens(nsites),spin(nsites),mu
  double precision occ(ndim),ndens(nsites),nspin(nsites),fill
  double precision t(ndim,ndim),w(ndim),work(3*ndim-1),pi
  character jobz,uplo

  open(unit=1111,file='bands.dat')

  pi=2.0d0*dasin(1.0d0)
  hubu=8.0d0
  jobz='V'
  uplo='U'
  fill=0.875d0
  do i=1,nx
     do j=1,ny
        site=i+(j-1)*nx
        dens(site)=fill+0.1d0*dcos(4.0d0*pi*i/q)
        spin(site)=0.1d0*dcos(2.0d0*pi*i/q)*(-1.0d0)**(i+j)
     end do
  end do

  count=0
  
666 continue

  flag=1
  mu=0.5d0*hubu*fill

  t=0.0d0
  do i=1,nx
     do j=1,ny
        site=i+(j-1)*nx
        if(i.lt.nx)t(site,site+1)=-1.0d0
        if(i.gt.1)t(site,site-1)=-1.0d0
        if(j.lt.ny)t(site,site+nx)=-1.0d0
        if(j.gt.1)t(site,site-nx)=-1.0d0
        t(site,site)=0.5d0*hubu*(dens(site)-spin(site))-mu
     end do
  end do

  do i=1,nx
     do j=1,ny
        site=nsites+i+(j-1)*nx
        if(i.lt.nx)t(site,site+1)=1.0d0
        if(i.gt.1)t(site,site-1)=1.0d0
        if(j.lt.ny)t(site,site+nx)=1.0d0
        if(j.gt.1)t(site,site-nx)=1.0d0
        t(site,site)=-0.5d0*hubu*(dens(site-nsites)+spin(site-nsites))+mu
     end do
  end do

!  do i=1,ndim
!     write(6,*)(t(i,j),j=1,ndim)
!  end do

!  goto 999
  
  call dsyev(jobz,uplo,ndim,t,ndim,w,work,3*ndim-1,info)
        
  occ=0.0d0
  do i=1,ndim
!     write(6,*)i,w(i)
     if(w(i).le.0.0d0)then
        do j=1,nsites
           occ(j)=occ(j)+t(j,i)**2
        end do
     else
        do j=1+nsites,ndim
           occ(j)=occ(j)+t(j,i)**2
        end do
     end if
  end do

  sumd=0.0d0
  do i=1,nsites
     ndens(i)=occ(i)+occ(i+nsites)
     nspin(i)=occ(i)-occ(i+nsites)
!     write(6,*)i,dens(i),spin(i)
     if(abs(ndens(i)-dens(i))/dens(i).lt.0.0001d0)then
        if(abs(nspin(i)-spin(i))/abs(spin(i)).lt.0.0001d0)then
           flag=0
        end if
     end if
     dens(i)=ndens(i)
     spin(i)=nspin(i)
!     write(6,*)i,dens(i),spin(i)
  end do
  count=count+1
  if(count.lt.cmin)goto 666
  if(count.lt.cmax)then
     if(flag.ne.0)goto 666
  end if

  write(6,*)count
!  gauss=1000.0d0

  fill=0.0d0
  do i=1,nsites
     fill=fill+dens(i)/dfloat(nsites)
  end do

  write(6,*)fill
  
  do i=1,nx
     do j=1,ny
        site=i+(j-1)*nx
        write(7,*)i,j,dens(site),spin(site)
     end do
     write(7,*)
  end do

999 continue
  stop
end program main

