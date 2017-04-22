program fdtd1d
  integer :: nz=900, mstep=1800
  real :: e(900)=0, h(900)=0
  real :: dz=1.0e-2
  real :: umu=1.257e-6, eps0=8.854e-12, c=2.998e8
  real :: sigma=0.0
  real :: pi=3.141592, freq=0.5e9
  real :: dt, ec1, ec2, hc
  character fname*4

  t=0.0

  dt=dz/c
  ec1=(1.0-sigma*dt/(2.0*eps0))/(1.0+sigma*dt/(2.0*eps0))
  ec2=-(dt/(eps0*dz)/(1.0+sigma*dt/(2.0*eps0)))
  hc=-dt/(dz*umu)

  do n=1,mstep

    if(t < 0.5/freq) e(500)=e(500)+sin(2.*pi*freq*t)**4

    do k=2,nz-1
      e(k)=ec1*e(k)+ec2*(h(k)-h(k-1))
    end do

    if(n == 1 .or. mod(n,5) == 0) then
      write(fname, '(i4.4)') n
      open(77, file = "data_1d/"//"e_"//fname//".txt")
      do k=1,nz
        write(77, *)e(k)
      end do
      close(77)
    endif

    t=t+dt/2.

    do k=1,nz-1
      h(k)=h(k)+hc*(e(k+1)-e(k))
    end do

    t=t+dt/2.

  end do
end program fdtd1d
