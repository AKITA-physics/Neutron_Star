!グローバル変数を指定するモジュール
module globals
    implicit none
    integer i,k
    real(8) r,theta,chi,phi,h,phi0,rho,m0
    real(8) d1,d2,d3,d4,d5,d6
    real(8) Bx,By,Bz,LBx,LBy,LBz,LB
    real(8) k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z,k4x,k4y,k4z
    real(8) x,y,z,dx,dy,dz,x0,y0,z0
    real(8) lc,lcx,lcy,lcz

!定数
    !光速
    real(8) :: c=29979245800.0d0
    !円周率の計算
    real(8) :: pi=acos(-1.0d0)

!パラメータ
    !何回繰り返し計算するか
    integer :: imax=5*10**5
    !何回に一回出力するのか
    integer :: iout=10*3
    !初期位置が自転軸から何度傾いているか
    real(8) :: alpha=acos(1.15d6/1.2d6)
    !距離ステップの取り方
    real(8) :: dh=1.0d3
    !今回考える時間
    real(8) :: t=0.0d0
    !中性子星の半径
    real(8) :: r0=1.2d6
    !星表面の赤道磁場
    real(8) :: B0=1.0d9
    !星の自転に対する角周波数
    real(8) :: omega=1256.637061d0
    !自転軸と磁軸のなす角の係数
    real(8) :: anglec=1.0d0/6.0d0

!星がchiだけ傾いているときの初期位置
contains
    subroutine xp(chi,x0)
            real(8)  chi,x0
            x0=r0*sin(alpha)*cos(phi0*pi/5.0d0)*cos(chi)+r0*cos(alpha)*sin(chi)
    end subroutine xp

    subroutine yp(y0)
            real(8)  y0
            y0=r0*sin(alpha)*sin(phi0*pi/5.0d0)
    end subroutine yp

    subroutine zp(chi,z0)
            real(8) chi,z0
            z0=-r0*sin(alpha)*cos(phi0*pi/5.0d0)*sin(chi)+r0*cos(alpha)*cos(chi)
    end subroutine zp
end module globals

!双極子磁場
module dipole
    use globals
    implicit none
contains
    !星中心からの距離
    subroutine radious(x,y,z,r)
        real(8) :: x,y,z,r
        r=sqrt(x**2.0d0+y**2.0d0+z**2.0d0)
    end subroutine radious
    !自転軸に対する偏角
    subroutine atack_angle(r,z,theta)
        real(8) :: r,z,theta
        theta=acos(z/r)
    end subroutine atack_angle
    !自転軸に垂直な方向にrを射影する
    subroutine xyradious(theta,rho)
        real(8) :: theta,rho
        rho=r*sin(theta)
    end subroutine xyradious
    !双極子磁場の係数１
    subroutine fd1(d1)
        real(8) :: d1
        d1=1.0d0
    end subroutine fd1
    !双極子磁場の係数２
    subroutine fd2(r,d2)
        real(8) :: r,d2
        d2=-r*omega/c
    end subroutine fd2
    !双極子磁場の係数３
    subroutine fd3(r,d3)
        real(8) :: r,d3
        d3=1.0d0-r*r*omega*omega/(c*c)
    end subroutine fd3
    !双極子磁場の係数４
    subroutine fd4(r,d4)
        real(8) :: r,d4
        d4=-r*omega/c
    end subroutine fd4
    !双極子磁場の係数５
    subroutine fd5(r,d5)
        real(8) :: r,d5
        d5=-r*omega/c
    end subroutine fd5
    !双極子磁場の係数６
    subroutine fd6(r,d6)
        real(8) :: r,d6
        d6=-1.0d0+r*r*omega*omega/(c*c)
    end subroutine fd6
    !双極子磁場のｘ成分
    subroutine bfx(x,y,z,r,rho,d1,d2,d3,d4,d5,d6,Bx)
        real(8) :: x,y,z,r,rho,d1,d2,d3,d4,d5,d6,Bx
        Bx=m0/(r*r*r)*(3.0d0*x*z*cos(chi)/(r*r)+2.0d0*sin(chi)*x/r&
            &*(d1*x*cos(omega*(t-r/c))/r+d1*y*sin(omega*(t-r/c))/r+d2*x*sin(omega*(t-r/c))/r-d2*y*cos(omega*(t-r/c))/r)&
            &-sin(chi)*(cos(omega*(t-r/c))*x/rho+sin(omega*(t-r/c))*y/rho)*(d3*z*z*x/(r*r*rho)+d5*y/rho)&
            &-sin(chi)*(sin(omega*(t-r/c))*x/rho-cos(omega*(t-r/c))*y/rho)*(d4*z*z*x/(r*r*rho)+d6*y/rho))
    end subroutine bfx
    !双極子磁場のｙ成分
    subroutine bfy(x,y,z,r,rho,d1,d2,d3,d4,d5,d6,By)
        real(8) :: x,y,z,r,rho,d1,d2,d3,d4,d5,d6,By
        By=m0/(r*r*r)*(3.0d0*y*z*cos(chi)/(r*r)+2.0d0*sin(chi)*y/r&
            &*(d1*x*cos(omega*(t-r/c))/r+d1*y*sin(omega*(t-r/c))/r+d2*x*sin(omega*(t-r/c))/r-d2*y*cos(omega*(t-r/c))/r)&
            &-sin(chi)*(cos(omega*(t-r/c))*x/rho+sin(omega*(t-r/c))*y/rho)*(d3*z*z*y/(r*r*rho)-d5*x/rho)&
            &-sin(chi)*(sin(omega*(t-r/c))*x-cos(omega*(t-r/c))*y)/rho*(d4*z*z*y/(r*r*rho)-d6*x/rho))
    end subroutine bfy
    !双極子磁場のｚ成分
    subroutine bfz(x,y,z,r,d1,d2,d3,d4,Bz)
        real(8) :: x,y,z,r,d1,d2,d3,d4,Bz
        Bz=m0/(r*r*r)*((3.0d0*z*z/(r*r)-1.0d0)*cos(chi)+sin(chi)*z/r&
            &*((2.0d0*d1+d3)*(x*cos(omega*(t-r/c))/r+y*sin(omega*(t-r/c))/r)&
            &+(2.0d0*d2+d4)*(x*sin(omega*(t-r/c))/r-y*cos(omega*(t-r/c))/r)))
    end subroutine bfz

end module dipole

!双極子磁場への代入
module sub_dipole
    use globals
    use dipole
contains
    subroutine assignment_dipole(x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz)
            real(8) :: x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz
            call radious(x,y,z,r)
            call atack_angle(r,z,theta)
            call xyradious(theta,rho)
            call fd1(d1)
            call fd2(r,d2)
            call fd3(r,d3)
            call fd4(r,d4)
            call fd5(r,d5)
            call fd6(r,d6)
            call bfx(x,y,z,r,rho,d1,d2,d3,d4,d5,d6,Bx)
            call bfy(x,y,z,r,rho,d1,d2,d3,d4,d5,d6,By)
            call bfz(x,y,z,r,d1,d2,d3,d4,Bz)
    end subroutine assignment_dipole
end module sub_dipole

!メインプログラム-双極子磁場の磁力線を描くためのプログラム
program main
    !使うモジュールの指定
    use globals
    use dipole
    use sub_dipole
    implicit none
    !ファイルの用意
    open(11, file='lineB0.d')
    open(12, file='lineB1.d')
    open(13, file='lineB2.d')
    open(14, file='lineB3.d')
    open(15, file='lineB4.d')
    open(16, file='lineB5.d')
    open(17, file='lineB6.d')
    open(18, file='lineB7.d')
    open(19, file='lineB8.d')
    open(20, file='lineB9.d')
!位置を変化させるためのdoループ
do k = 1,10
    !初期位置の方位角を決めるパラメータ
    phi0=dble(k)-1.0d0
    !磁手軸と磁軸のなす角
    chi=anglec*pi
    !双極子モーメント
    m0=r0**3*B0
    !光円柱半径
    lc=c/omega

    !初期位置の呼び出し
    call xp(chi,x0)
    call yp(y0)
    call zp(chi,z0)
    !初期位置を光円柱半径で規格化
    x=x0/lc
    y=y0/lc
    z=z0/lc

    !4次のrunge-kutta法
    do i = 1,imax
            !規格化された位置を元のスケールに戻す
            x=x*lc
            y=y*lc
            z=z*lc
            !双極子磁場の呼び出し
            call assignment_dipole(x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz)
            !磁場の強度
            LB=sqrt(Bx**2.0d0+By**2.0d0+Bz**2.0d0)

            k1x=dh*Bx/LB
            k1y=dh*By/LB
            k1z=dh*Bz/LB

            x=x+0.5d0*k1x
            y=y+0.5d0*k1y
            z=z+0.5d0*k1z
            !双極子磁場の呼び出し
            call assignment_dipole(x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz)
            !磁場の強度
            LB=sqrt(Bx**2.0d0+By**2.0d0+Bz**2.0d0)

            k2x=dh*Bx/LB
            k2y=dh*By/LB
            k2z=dh*Bz/LB

            x=x+0.5d0*k2x-0.5d0*k1x
            y=y+0.5d0*k2y-0.5d0*k1y
            z=z+0.5d0*k2z-0.5d0*k1z
            !双極子磁場の呼び出し
            call assignment_dipole(x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz)
            !磁場の強度
            LB=sqrt(Bx**2.0d0+By**2.0d0+Bz**2.0d0)

            k3x=dh*Bx/LB
            k3y=dh*By/LB
            k3z=dh*Bz/LB

            x=x+k3x-0.5d0*k2x
            y=y+k3y-0.5d0*k2y
            z=z+k3z-0.5d0*k2z
            !双極子磁場の呼び出し
            call assignment_dipole(x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz)
            !磁場の強度
            LB=sqrt(Bx**2.0d0+By**2.0d0+Bz**2.0d0)

            k4x=dh*Bx/LB
            k4y=dh*By/LB
            k4z=dh*Bz/LB
            !微小距離の導出
            dx=(k1x+2.0d0*k2x+2.0d0*k3x+k4x)/6.0d0
            dy=(k1y+2.0d0*k2y+2.0d0*k3y+k4y)/6.0d0
            dz=(k1z+2.0d0*k2z+2.0d0*k3z+k4z)/6.0d0
            !規格化された距離の計算後
            x=(x-k3x)/lc+dx/lc
            y=(y-k3y)/lc+dy/lc
            z=(z-k3z)/lc+dz/lc
            !結果の書き出し
            if(mod(i-1,iout) == 0)then
                write(10+k,*) x,y,z,r,Bx,By,Bz
            end if
            !星中心からの距離が一定以上より小さいと繰り返し計算を出る。
            if(r < 1.1d6)exit
    end do
end do
!ファイルを閉じる
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