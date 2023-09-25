!中性子星周りの四重極電場
!グローバル変数
module globals
    implicit none
    integer i,k
    real(8) r,theta,chi,phi,h,phi0,m0,E0,rho
    real(8) Ex,Ey,Ez,LE
    real(8) d1,d2,d3,d4,d5,d6
    real(8) k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z,k4x,k4y,k4z
    real(8) x,y,z,dx,dy,dz,x0,y0,z0
    real(8) lc,lcx,lcy,lcz

    !定数
    real(8) :: c=29979245800.0d0
    real(8) :: pi=acos(-1.0d0)

    !パラメータ
    integer :: imax=5*10**7
    integer :: iout=3*10*7
    real(8) :: alpha=acos(0.4d6/1.2d6)!acos(sqrt(1.2d6*1.2d6-1.0d5*1.0d5)/1.2d6)!0.78539816339!1.3089969389957!1.832595714594!
    real(8) :: dh=1.0d1
    real(8) :: B0=1.0d9
    real(8) :: t=0.0025
    real(8) :: r0=1.2d6
    real(8) :: omega=1256.637061d0
    real(8) :: anglec=1.0d0/6.0d0
    
contains
    subroutine xp(chi,x0)
            real(8)  chi,x0
            x0=1.5d0*cos(phi0*pi/5.0d0)!r0*sin(alpha)*cos(phi0*pi/5.0d0)*cos(chi)+r0*cos(alpha)*sin(chi)
    end subroutine xp

    subroutine yp(y0)
            real(8)  y0
            y0=1.5d0*sin(phi0*pi/5.0d0)!r0*sin(alpha)*sin(phi0*pi/5.0d0)
    end subroutine yp

    subroutine zp(chi,z0)
            real(8) chi,z0
            z0=1.5d0*sin(alpha)!-r0*sin(alpha)*cos(phi0*pi/5.0d0)*sin(chi)+r0*cos(alpha)*cos(chi)
    end subroutine zp
end module globals

!四重極電場(時間変化あり)
module quadrupole
    use globals
    implicit none
contains
    subroutine radious(x,y,z,r)
            real(8) :: x,y,z,r
            r=sqrt(x**2.0d0+y**2.0d0+z**2.0d0)
    end subroutine radious

    subroutine atack_angle(r,z,theta)
            real(8) :: r,z,theta
            theta=acos(z/r)
    end subroutine atack_angle

    subroutine xyradious(theta,rho)
        real(8) :: theta,rho
        rho=r*sin(theta)
    end subroutine xyradious

    subroutine fd1(d1)
        real(8) :: d1
        d1=1.0d0
    end subroutine fd1

    subroutine fd2(r,d2)
        real(8) :: r,d2
        d2=-r*omega/c
    end subroutine fd2

    subroutine efx(x,y,z,r,rho,d1,d2,Ex)
        real(8) :: x,y,z,r,rho,d1,d2,Ex
        Ex=E0*r0**4*x/(r**5)*cos(chi)*(1.0d0-5.0d0*z**2/r**2)&
            &-E0*r0**2*z/r**3*sin(chi)*((d1*x/rho-d2*y/rho)*(x*cos(omega*(t-r/c))/rho+y*sin(omega*(t-r/c))/rho)&
            &+(d2*x/rho+d1*y/rho)*(x*sin(omega*(t-r/c))/rho-y*cos(omega*(t-r/c))/rho))
    end subroutine efx

    subroutine efy(x,y,z,r,rho,d1,d2,Ey)
        real(8) :: x,y,z,r,rho,d1,d2,Ey
        Ey=E0*r0**4*y/(r**5)*cos(chi)*(1.0d0-5.0d0*z**2/r**2)&
            &-E0*r0**2*z/r**3*sin(chi)*((d1*y/rho+d2*x/rho)*(x*cos(omega*(t-r/c))/rho+y*sin(omega*(t-r/c))/rho)&
            &+(d2*y/rho-d1*x/rho)*(x*sin(omega*(t-r/c))/rho-y*cos(omega*(t-r/c))/rho))
    end subroutine efy

    subroutine efz(x,y,z,r,rho,d1,d2,Ez)
        real(8) :: x,y,z,r,rho,d1,d2,Ez
        Ez=E0*r0**4*z/(r**5)*cos(chi)*(3.0d0-5.0d0*z**2/r**2)&
            &+E0*r0**2/r**2*rho/r*sin(chi)*(d1*(x*cos(omega*(t-r/c))/rho+y*sin(omega*(t-r/c))/rho)&
            &+d2*(x*sin(omega*(t-r/c))/rho-y*cos(omega*(t-r/c))/rho))
    end subroutine efz

end module quadrupole

!四重極電場(時間変化あり)への代入
module sub_quadrupole
    use globals
    use quadrupole
contains
    subroutine assignment_dipole(x,y,z,r,theta,rho,d1,d2,Ex,Ey,Ez)
            real(8) :: x,y,z,r,theta,rho,d1,d2,Ex,Ey,Ez
            call radious(x,y,z,r)
            call atack_angle(r,z,theta)
            call xyradious(theta,rho)
            call fd1(d1)
            call fd2(r,d2)
            call efx(x,y,z,r,rho,d1,d2,Ex)
            call efy(x,y,z,r,rho,d1,d2,Ey)
            call efz(x,y,z,r,rho,d1,d2,Ez)
    end subroutine assignment_dipole
end module sub_quadrupole

!メインプログラム
program main
    use globals
    use quadrupole
    use sub_quadrupole
    implicit none

    open(11, file='2nelineE0.d')
    open(12, file='2nelineE1.d')
    open(13, file='2nelineE2.d')
    open(14, file='2nelineE3.d')
    open(15, file='2nelineE4.d')
    open(16, file='2nelineE5.d')
    open(17, file='2nelineE6.d')
    open(18, file='2nelineE7.d')
    open(19, file='2nelineE8.d')
    open(20, file='2nelineE9.d')
do k = 1,10
!k=10
    phi0=dble(k)-1.0d0
    chi=anglec*pi
    m0=r0**3*B0
    E0=omega*r0*B0
    lc=c/omega

    call xp(chi,x0)
    call yp(y0)
    call zp(chi,z0)

    !初期値
    x=x0!/lc
    y=y0!/lc
    z=z0!/lc

!    call assignment_dipole(x,y,z,r,theta,rho,d1,d2,Ex,Ey,Ez)
!    write(10+k,*) x,y,z,r,d2,Ex,Ey,Ez

    !runge-kutta法
    do i = 1,imax

            x=x*lc
            y=y*lc
            z=z*lc

            call assignment_dipole(x,y,z,r,theta,rho,d1,d2,Ex,Ey,Ez)

            LE=sqrt(Ex**2.0d0+Ey**2.0d0+Ez**2.0d0)

            k1x=dh*Ex/LE
            k1y=dh*Ey/LE
            k1z=dh*Ez/LE

            x=x+0.5d0*k1x
            y=y+0.5d0*k1y
            z=z+0.5d0*k1z

            call assignment_dipole(x,y,z,r,theta,rho,d1,d2,Ex,Ey,Ez)

            LE=sqrt(Ex**2.0d0+Ey**2.0d0+Ez**2.0d0)

            k2x=dh*Ex/LE
            k2y=dh*Ey/LE
            k2z=dh*Ez/LE

            x=x+0.5d0*k2x-0.5d0*k1x
            y=y+0.5d0*k2y-0.5d0*k1y
            z=z+0.5d0*k2z-0.5d0*k1z

            call assignment_dipole(x,y,z,r,theta,rho,d1,d2,Ex,Ey,Ez)

            LE=sqrt(Ex**2.0d0+Ey**2.0d0+Ez**2.0d0)

            k3x=dh*Ex/LE
            k3y=dh*Ey/LE
            k3z=dh*Ez/LE

            x=x+k3x-0.5d0*k2x
            y=y+k3y-0.5d0*k2y
            z=z+k3z-0.5d0*k2z

            call assignment_dipole(x,y,z,r,theta,rho,d1,d2,Ex,Ey,Ez)

            LE=sqrt(Ex**2.0d0+Ey**2.0d0+Ez**2.0d0)

            k4x=dh*Ex/LE
            k4y=dh*Ey/LE
            k4z=dh*Ez/LE

            dx=(k1x+2.0d0*k2x+2.0d0*k3x+k4x)/6.0d0/lc
            dy=(k1y+2.0d0*k2y+2.0d0*k3y+k4y)/6.0d0/lc
            dz=(k1z+2.0d0*k2z+2.0d0*k3z+k4z)/6.0d0/lc

            x=(x-k3x)/lc+dx
            y=(y-k3y)/lc+dy
            z=(z-k3z)/lc+dz

            if(mod(i-1,10**5) == 0)then
                    write(10+k,*) x,y,z,r
            end if

            if(r < 1.1d6)exit
            if(rho < 1.0d2)exit
    end do
end do 
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)

!光円柱の計算
    open(21, file='lc.d')
do i = 1,629
    phi=1.0d-2*dble(i)

    lc=c/omega
    lcx=lc*cos(phi)
    lcy=lc*sin(phi)
    lcz=0.0d0

    write(21,*)phi,lcx,lcy,lcz,lc
end do
close(21)
end program main
