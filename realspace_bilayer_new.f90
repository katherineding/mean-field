  program main
  implicit none
  integer, parameter :: nx=32
  integer, parameter :: ny=32
  integer, parameter :: nsites=nx*ny
  integer, parameter :: ndim=4*nsites
  integer, parameter :: cmax=10
  integer, parameter :: cmin=1
  integer, parameter :: q=16
  integer, parameter :: nfreq=100
  integer i,j,count,site,flag,k
  integer info
  double precision hubu,dens(2*nsites),spin(2*nsites),mu,mud,filld
  double precision occ(ndim),ndens(2*nsites),nspin(2*nsites),fill
  double precision t(ndim,ndim),w(ndim),work(3*ndim-1),pi,td,tz,var
  double precision total_e,dos1,dos2,w1,w2,gamma,omega,gauss
  character jobz,uplo

  open(unit=1111,file='bands.dat')

  pi=2.0d0*dasin(1.0d0)
  hubu=8.0d0
  jobz='V'
  uplo='U'
  fill=1.0d0
  var=0.0d0
  do i=1,nx
     do j=1,ny
        site=i+(j-1)*nx
!        dens(site)=fill*(1.0d0+var*dcos(4.0d0*pi*i/q))
!        spin(site)=var*dcos(2.0d0*pi*i/q)*(-1.0d0)**(i+j)
        dens(site)=fill
        spin(site)=(-1.0d0)**(i+j)
     end do
  end do
  filld=1.0d0
  do i=1,nx
     do j=1,ny
        site=nsites+i+(j-1)*nx
        dens(site)=filld+var*fill*dcos(4.0d0*pi*i/q)
        spin(site)=var*dcos(2.0d0*pi*i/q)*(-1.0d0)**(i+j)
     end do
  end do

  count=0
  
666 continue

  flag=1
  mu=0.5d0*hubu*fill
!  mud=-3.5d0
  mud=-3.5d0
  td=1.0d0
  tz=2.0d0

  t=0.0d0
  do i=1,nx
     do j=1,ny
        site=i+(j-1)*nx
        if(i.lt.nx)t(site,site+1)=-1.0d0
        if(i.gt.1)t(site,site-1)=-1.0d0
        if(j.lt.ny)t(site,site+nx)=-1.0d0
        if(j.gt.1)t(site,site-nx)=-1.0d0
        t(site,site)=0.5d0*hubu*(dens(site)-spin(site))-mu
        t(site,2*nsites+site)=-tz
     end do
  end do

  do i=1,nx
     do j=1,ny
        site=nsites+i+(j-1)*nx
        if(i.lt.nx)t(site,site+1)=-1.0d0
        if(i.gt.1)t(site,site-1)=-1.0d0
        if(j.lt.ny)t(site,site+nx)=-1.0d0
        if(j.gt.1)t(site,site-nx)=-1.0d0
        t(site,site)=0.5d0*hubu*(dens(site-nsites)+spin(site-nsites))-mu
        t(site,2*nsites+site)=-tz
     end do
  end do

  do i=1,nx
     do j=1,ny
        site=2*nsites+i+(j-1)*nx
        if(i.lt.nx)t(site,site+1)=-td
        if(i.gt.1)t(site,site-1)=-td
        if(j.lt.ny)t(site,site+nx)=-td
        if(j.gt.1)t(site,site-nx)=-td
        t(site,site)=-mud
        t(site,site-2*nsites)=-tz
     end do
  end do

  do i=1,nx
     do j=1,ny
        site=3*nsites+i+(j-1)*nx
        if(i.lt.nx)t(site,site+1)=-td
        if(i.gt.1)t(site,site-1)=-td
        if(j.lt.ny)t(site,site+nx)=-td
        if(j.gt.1)t(site,site-nx)=-td
        t(site,site)=-mud
        t(site,site-2*nsites)=-tz
     end do
  end do

!  do i=1,ndim
!     write(6,*)(t(i,j),j=1,ndim)
!  end do

!  goto 999
  
  call dsyev(jobz,uplo,ndim,t,ndim,w,work,3*ndim-1,info)
        
  occ=0.0d0
  total_e=0.0d0
  do i=1,ndim
     if(w(i).le.0.0d0)then
        total_e=total_e+w(i)/dfloat(ndim)
        do j=1,4*nsites
           occ(j)=occ(j)+t(j,i)**2
        end do
     end if
  end do

  do i=1,nsites
     ndens(i)=occ(i)+occ(i+nsites)
     nspin(i)=occ(i)-occ(i+nsites)
     ndens(i+nsites)=occ(i+2*nsites)+occ(i+3*nsites)
     nspin(i+nsites)=occ(i+2*nsites)-occ(i+3*nsites)
!     write(6,*)i,dens(i),spin(i)
     if(abs(ndens(i)-dens(i))/dens(i).lt.0.0001d0)then
        if(abs(nspin(i)-spin(i))/abs(spin(i)).lt.0.0001d0)then
           flag=0
        end if
     end if
     dens(i)=ndens(i)
     spin(i)=nspin(i)
     dens(i+nsites)=ndens(i+nsites)
     spin(i+nsites)=nspin(i+nsites)
!     write(6,*)i,dens(i),spin(i)
  end do
  fill=0.0d0
  filld=0.0d0
  do i=1,nsites
     fill=fill+dens(i)/dfloat(nsites)
     filld=filld+dens(i+nsites)/dfloat(nsites)
  end do

  write(6,*)fill,filld,total_e
  
  count=count+1
  if(count.lt.cmin)goto 666
  if(count.lt.cmax)then
     if(flag.ne.0)goto 666
  end if

  write(6,*)count

  do i=1,nx
     do j=1,ny
        site=i+(j-1)*nx
        write(8,*)i,j,dens(site),spin(site)
     end do
     write(8,*)
  end do

  do i=1,nx
     do j=1,ny
        site=nsites+i+(j-1)*nx
        write(9,*)i,j,dens(site),spin(site)
     end do
     write(9,*)
  end do

  gauss=10.0d0
  do j=-nfreq,nfreq
     omega=dfloat(10*j)/dfloat(nfreq)
     dos1=0.0d0
     dos2=0.0d0
     do i=1,ndim
        w1=0.0d0
        w2=0.0d0
        do k=1,2*nsites
           w1=w1+t(k,i)**2
           w2=w2+t(k+2*nsites,i)**2
        end do
        dos1=dos1+w1*dexp(-gauss*(omega-w(i))**2)/dfloat(ndim)
        dos2=dos2+w2*dexp(-gauss*(omega-w(i))**2)/dfloat(ndim)
     end do
     write(10,*)omega,dos1,dos2
  end do

999 continue
  stop
end program main

