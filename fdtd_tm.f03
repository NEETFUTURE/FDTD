program fdtd_tm

  implicit none

  integer :: i,j,n
  integer :: nx=60, ny=60, mstep=1800
  double precision,allocatable,dimension(:,:) :: ez
  double precision,allocatable,dimension(:,:) :: hx
  double precision,allocatable,dimension(:,:) :: hy
  
  double precision :: dx=0.01
  double precision :: dy=0.01
  double precision :: umu=1.257e-6, eps0=8.854e-12, c=2.998e8
  double precision :: sigma=0.0
  double precision :: pi=3.141592, freq=0.5e9
  double precision :: dt, ec1, ec2x, ec2y, hcx, hcy
  double precision :: t
  character fname*10

  allocate(ez(nx,ny))
  allocate(hx(nx,ny-1))
  allocate(hy(nx-1,ny))

  dt=0.1/(c*sqrt(1/(dx**2) + 1/(dy**2)))
  ec1=(1.0-sigma*dt/(2.0*eps0))/(1.0+sigma*dt/(2.0*eps0))
  ec2x=(dt/(eps0*dx)/(1.0+sigma*dt/(2.0*eps0)))
  ec2y=(dt/(eps0*dy)/(1.0+sigma*dt/(2.0*eps0)))
  hcx=-dt/(dx*umu)
  hcy=-dt/(dy*umu)

  t=0.0

  !配列初期化
  do i=1,nx
    do j=1,ny
      ez(i,j) = 0
      if(j<nx) hx(i,j) = 0
      if(i<ny) hy(i,j) = 0
    end do
  end do

  do n=1,mstep
    
    if(t < 0.5/freq) ez(nx/2,ny/2)=ez(10,10)+sin(2.*pi*freq*t)**4

    !電界を計算
    do i=2,nx-1
      do j=2,ny-1
        ez(i,j) = ec1*ez(i,j) &
          & + ec2x * (hy(i,j) - hy(i-1,j)) &
          & - ec2y * (hx(i,j) - hx(i,j-1))
      end do
    end do

    !結果をファイルに出力
!    if(n == 1 .or. mod(n,5) == 0) then
!      write(fname, "('e_',i4.4,'.txt')") n
!      open(77, file = "data_tm/"//fname, status='replace', action='write')
!
!      do i=1,nx
!        do j=1,ny
!          write(77, *)ez(i,j)
!        end do
!      end do
!      close(77)
!    endif
    if(n == 1 .or. mod(n,5) == 0) then
      write(fname, "('e_',i4.4,'.raw')") n
      open(77, file = "data_tm/"//fname, status='replace', action='write',form='unformatted', access='stream')
      write(77) ez
      close(77)
    end if

    t=t+dt/2.

    !磁界を計算
    do i=1,nx
      do j = 1,ny-1
        hx(i,j) = hx(i,j) + hcx * (ez(i,j+1) - ez(i,j))
      end do
    end do
    do i=1,nx-1
      do j = 1,ny
        hy(i,j) = hy(i,j) - hcy * (ez(i+1,j) - ez(i,j))
      end do
    end do

    t=t+dt/2.

  end do
end program fdtd_tm

