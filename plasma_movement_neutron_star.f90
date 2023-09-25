!中性子星周りでのプラズマ粒子の運動の解法
!グローバル変数
module globals
        implicit none
        real(8) r,omega
        real(8) theta
        real(8) d1,d2,d3,d4,d5,d6
        real(8) Bx,By,Bz
        real(8) E0,rho
        real(8) Ex,Ey,Ez
        real(8) chi,m0,r0
        real(8) x,y,z,x0,y0,z0,t

        !定数
        real(8) :: c=2.99792458d0*1.0d10
        real(8) :: e=-4.803204673*1.0d-10
        real(8) :: me=9.1093837*1.0d-28
        real(8) :: mp=1.67262192369*1.0d-24
        real(8) :: pi=acos(-1.0d0)     

        !パラメータ
        integer :: imax=5*10**8
        integer :: iout=5*10**6
        real(8) :: gamma0=1.0d3
        real(8) :: dt=1.0d-8
        real(8) :: P=0.005d0
        real(8) :: Rns=1.2d6
        real(8) :: Bs=1.0d9
        real(8) :: phi0=acos(sqrt(0.99d0))

contains
        subroutine radious0(r0)
                real(8) r0
                r0=Rns
        end subroutine radious0

        subroutine xp(r0,chi,omega,x0)
!                real(8)  chi,omega,x0
!                x0=c*sin(chi+phi0)/omega
                real(8)  r0,chi,omega,x0
                x0=sqrt(r0*c/omega)*sin(chi+phi0)
!                real(8)  r0,chi,x0
!                x0=r0*sin(chi+phi0)
       end subroutine xp

        subroutine yp(y0)
                real(8)  y0
                y0=0.0d0
        end subroutine yp

        subroutine zp(r0,chi,omega,z0)
!                real(8) chi,omega,z0
!                z0=c*cos(chi+phi0)/omega
                real(8) r0,chi,omega,z0
                z0=sqrt(r0*c/omega)*cos(chi+phi0)
!                real(8) r0,chi,z0
!                z0=r0*cos(chi+phi0)
        end subroutine zp
end module globals

!電磁場の導出
module dipole
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

    subroutine fd3(r,d3)
        real(8) :: r,d3
        d3=1.0d0-r*r*omega*omega/(c*c)
    end subroutine fd3

    subroutine fd4(r,d4)
        real(8) :: r,d4
        d4=-r*omega/c
    end subroutine fd4

    subroutine fd5(r,d5)
        real(8) :: r,d5
        d5=-r*omega/c
    end subroutine fd5

    subroutine fd6(r,d6)
        real(8) :: r,d6
        d6=-1.0d0+r*r*omega*omega/(c*c)
    end subroutine fd6

    subroutine bfx(t,x,y,z,r,rho,d1,d2,d3,d4,d5,d6,Bx)
        real(8) :: t,x,y,z,r,rho,d1,d2,d3,d4,d5,d6,Bx
        Bx=m0/(r*r*r)*(3.0d0*x*z*cos(chi)/(r*r)+2.0d0*sin(chi)*x/r&
            &*(d1*x*cos(omega*(t-r/c))/r+d1*y*sin(omega*(t-r/c))/r+d2*x*sin(omega*(t-r/c))/r-d2*y*cos(omega*(t-r/c))/r)&
            &-sin(chi)*(cos(omega*(t-r/c))*x/rho+sin(omega*(t-r/c))*y/rho)*(d3*z*z*x/(r*r*rho)+d5*y/rho)&
            &-sin(chi)*(sin(omega*(t-r/c))*x/rho-cos(omega*(t-r/c))*y/rho)*(d4*z*z*x/(r*r*rho)+d6*y/rho))
    end subroutine bfx

    subroutine bfy(t,x,y,z,r,rho,d1,d2,d3,d4,d5,d6,By)
        real(8) :: t,x,y,z,r,rho,d1,d2,d3,d4,d5,d6,By
        By=m0/(r*r*r)*(3.0d0*y*z*cos(chi)/(r*r)+2.0d0*sin(chi)*y/r&
            &*(d1*x*cos(omega*(t-r/c))/r+d1*y*sin(omega*(t-r/c))/r+d2*x*sin(omega*(t-r/c))/r-d2*y*cos(omega*(t-r/c))/r)&
            &-sin(chi)*(cos(omega*(t-r/c))*x/rho+sin(omega*(t-r/c))*y/rho)*(d3*z*z*y/(r*r*rho)-d5*x/rho)&
            &-sin(chi)*(sin(omega*(t-r/c))*x-cos(omega*(t-r/c))*y)/rho*(d4*z*z*y/(r*r*rho)-d6*x/rho))
    end subroutine bfy

    subroutine bfz(t,x,y,z,r,d1,d2,d3,d4,Bz)
        real(8) :: t,x,y,z,r,d1,d2,d3,d4,Bz
        Bz=m0/(r*r*r)*((3.0d0*z*z/(r*r)-1.0d0)*cos(chi)+sin(chi)*z/r&
            &*((2.0d0*d1+d3)*(x*cos(omega*(t-r/c))/r+y*sin(omega*(t-r/c))/r)&
            &+(2.0d0*d2+d4)*(x*sin(omega*(t-r/c))/r-y*cos(omega*(t-r/c))/r)))
    end subroutine bfz

    subroutine efx(t,x,y,z,r,rho,d1,d2,Ex)
        real(8) :: t,x,y,z,r,rho,d1,d2,Ex
        Ex=E0*r0**4*x/(r**5)*cos(chi)*(1.0d0-5.0d0*z**2/r**2)&
            &-E0*r0**2*z/r**3*sin(chi)*((d1*x/rho-d2*y/rho)*(x*cos(omega*(t-r/c))/rho+y*sin(omega*(t-r/c))/rho)&
            &+(d2*x/rho+d1*y/rho)*(x*sin(omega*(t-r/c))/rho-y*cos(omega*(t-r/c))/rho))
    end subroutine efx

    subroutine efy(t,x,y,z,r,rho,d1,d2,Ey)
        real(8) :: t,x,y,z,r,rho,d1,d2,Ey
        Ey=E0*r0**4*y/(r**5)*cos(chi)*(1.0d0-5.0d0*z**2/r**2)&
            &-E0*r0**2*z/r**3*sin(chi)*((d1*y/rho+d2*x/rho)*(x*cos(omega*(t-r/c))/rho+y*sin(omega*(t-r/c))/rho)&
            &+(d2*y/rho-d1*x/rho)*(x*sin(omega*(t-r/c))/rho-y*cos(omega*(t-r/c))/rho))
    end subroutine efy

    subroutine efz(t,x,y,z,r,rho,d1,d2,Ez)
        real(8) :: t,x,y,z,r,rho,d1,d2,Ez
        Ez=E0*r0**4*z/(r**5)*cos(chi)*(3.0d0-5.0d0*z**2/r**2)&
            &+E0*r0**2/r**2*rho/r*sin(chi)*(d1*(x*cos(omega*(t-r/c))/rho+y*sin(omega*(t-r/c))/rho)&
            &+d2*(x*sin(omega*(t-r/c))/rho-y*cos(omega*(t-r/c))/rho))
    end subroutine efz

end module dipole

!電磁場への代入
module sub_dipole
    use globals
    use dipole
contains
        subroutine assignment_dipole(t,x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz,Ex,Ey,Ez)
                real(8) :: t,x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz,Ex,Ey,Ez
                call radious(x,y,z,r)
                call atack_angle(r,z,theta)
                call xyradious(theta,rho)
                call fd1(d1)
                call fd2(r,d2)
                call fd3(r,d3)
                call fd4(r,d4)
                call fd5(r,d5)
                call fd6(r,d6)
                call bfx(t,x,y,z,r,rho,d1,d2,d3,d4,d5,d6,Bx)
                call bfy(t,x,y,z,r,rho,d1,d2,d3,d4,d5,d6,By)
                call bfz(t,x,y,z,r,d1,d2,d3,d4,Bz)
                call efx(t,x,y,z,r,rho,d1,d2,Ex)
                call efy(t,x,y,z,r,rho,d1,d2,Ey)
                call efz(t,x,y,z,r,rho,d1,d2,Ez)
        end subroutine assignment_dipole
end module sub_dipole

program main
        use globals
        use dipole
        use sub_dipole
        implicit none
        integer i,j
        real(8) rr,s
        real(8) LB0
        real(8) b
        real(8) NE(1:3),NB(1:3),LB(1:3),LE(1:3),u(1:4),xmat(1:4)
        real(8) bt(1:3),um(1:4),ut(1:4),utt(1:4),up(1:4)
        real(8) VB(1:3),K,I1,I2
        real(8) V(1:4),F(1:3)

        open(10,file='pnR0.d')
        open(11,file='pnR60.d')
        open(12,file='pnR120.d')
        open(13,file='pnR180.d')
        open(14,file='pnR240.d')
        open(15,file='pnR300.d')
        open(20,file='pnRr0.d')
        open(21,file='pnRr60.d')
        open(22,file='pnRr120.d')
        open(23,file='pnRr180.d')
        open(24,file='pnRr240.d')
        open(25,file='pnRr300.d')
        open(30,file='pnrL0.d')
        open(31,file='pnrL60.d')
        open(32,file='pnrL120.d')
        open(33,file='nrL180.d')
        open(34,file='penrL240.d')
        open(35,file='pnrL300.d')

!do j = 1,2
j=3
        !星の傾き
        chi=dble(j-1)*pi/3.0d0

        !規格化するためのシンクロトロン周波数
        omega=2.0d0*pi/P
        
        !初期位置の呼び出し
        call radious0(r0)
!        call xp(r0,chi,x0)
!        call yp(y0)
!        call zp(r0,chi,z0)
        call xp(r0,chi,omega,x0)
        call yp(y0)
        call zp(r0,chi,omega,z0)
!        call xp(chi,omega,x0)
!        call yp(y0)
!        call zp(chi,omega,z0)

        !磁気モーメントの大きさ
        m0=Rns**3*Bs

        !表面電場の大きさ
        E0=omega*Rns*Bs/c

        !初期位置
        t=0.0d0
        x=x0
        y=y0
        z=z0

        !電磁場の呼び出し
        call assignment_dipole(t,x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz,Ex,Ey,Ez)

        !電磁場
        LB(1)=Bx
        LB(2)=By
        LB(3)=Bz
        LE(1)=Ex
        LE(2)=Ey
        LE(3)=Ez

        !位置を行列形式へ
        xmat(1)=t
        xmat(2)=x0*omega/c
        xmat(3)=y0*omega/c
        xmat(4)=z0*omega/c

        !初速度
        u(1)=gamma0
        u(2)=0.0d0
        u(3)=sqrt(u(1)*u(1)-1.0d0)
        u(4)=0.0d0

        !放射減衰項の係数
        b=omega*2.0d0*e**2/(3.0d0*mp*c**3)

        !規格化するための磁場
        LB0=abs(omega*mp*c/e)

!                write(9+k,*) xmat(1:4),x,y,z,r,sqrt(dot_product(LB(1:3),LB(1:3)))

        do i = 1,imax
                !規格化
                NE(1:3)=LE(1:3)/LB0
                NB(1:3)=LB(1:3)/LB0
                !四元速度から速度への変換
                V(1:4)=u(1:4)/u(1)
                !V×B
                VB(1)=V(3)*LB(3)-V(4)*LB(2)
                VB(2)=V(4)*LB(1)-V(2)*LB(3)
                VB(3)=V(2)*LB(2)-V(3)*LB(1)
                !ローレンツ力
                F(1)=LE(1)+VB(1)
                F(2)=LE(2)+VB(2)
                F(3)=LE(3)+VB(3)
                !電磁場テンソルの固有値
                I1=dot_product(LE(1:3),LE(1:3))-dot_product(LB(1:3),LB(1:3))
                I2=dot_product(LE(1:3),LB(1:3))
                K=-sqrt((I1+sqrt(I1*I1+4.0d0*I2*I2))/2.0d0)

                !放射減衰の第二項目の係数の計算
                rr=u(1)*(-dot_product(NE(1:3),NE(1:3))*u(1)+(NE(2)*NB(3)-NE(3)*NB(2))*u(2)+(NE(3)*NB(1)-NE(1)*NB(3))*u(3)&
                &+(NE(1)*NB(2)-NE(2)*NB(1))*u(4))&
                &+u(2)*((NE(2)*NB(3)-NE(3)*NB(2))*u(1)+(NE(1)**2.0d0-NB(3)**2.0d0-NB(2)**2.0d0)*u(2)&
                &+(NE(1)*NE(2)+NB(1)*NB(2))*u(3)+(NE(1)*NE(3)+NB(1)*NB(3))*u(4))&
                &+u(3)*((NE(3)*NB(1)-NE(1)*NB(3))*u(1)+(NE(1)*NE(2)+NB(1)*NB(2))*u(2)&
                &+(NE(2)**2.0d0-NB(3)**2.0d0-NB(1)**2.0d0)*u(3)+(NE(2)*NE(3)+NB(2)*NB(3))*u(4))&
                &+u(4)*((NE(1)*NB(2)-NE(2)*NB(1))*u(1)+(NE(1)*NE(3)+NB(1)*NB(3))*u(2)&
                &+(NE(2)*NE(3)+NB(2)*NB(3))*u(3)+(NE(3)**2.0d0-NB(1)**2.0d0-NB(2)**2.0d0)*u(4))
                !空間成分の半ステップの電場の加速と放射減衰
                um(2)=u(2)+0.5d0*NE(1)*dt&
                        &+0.5d0*b*((NE(2)*NB(3)-NE(3)*NB(2))*u(1)+(NE(1)*NE(1)-NB(2)*NB(2)-NB(3)*NB(3))*u(2)&
                        &+(NE(1)*NE(2)+NB(1)*NB(2))*u(3)+(NE(1)*NE(3)+NB(1)*NB(3))*u(4)&
                        &+rr*u(2))*dt/u(1)
                um(3)=u(3)+0.5d0*NE(2)*dt&
                        &+0.5d0*b*((NE(3)*NB(1)-NE(1)*NB(3))*u(1)+(NE(1)*NE(2)+NB(1)*NB(2))*u(2)&
                        &+(NE(2)*NE(2)-NB(1)*NB(1)-NB(3)*NB(3))*u(3)+(NE(2)*NE(3)+NB(2)*NB(3))*u(4)&
                        &+rr*u(3))*dt/u(1)
                um(4)=u(4)+0.5d0*NE(3)*dt&
                        &+0.5d0*b*((NE(1)*NB(2)-NE(2)*NB(1))*u(1)+(NE(1)*NE(3)+NB(1)*NB(3))*u(2)&
                        &+(NE(2)*NE(3)+NB(2)*NB(3))*u(3)+(NE(3)*NE(3)-NB(1)*NB(1)-NB(2)*NB(2))*u(4)&
                        &+rr*u(4))*dt/u(1)
                !ローレンツ因子の書き換え
                um(1)=sqrt(1.0d0+um(2)*um(2)+um(3)*um(3)+um(4)*um(4))
                !磁場による回転
                bt(1:3)=0.5d0*NB(1:3)*dt/um(1)
                s=2.0d0/(1.0d0+dot_product(bt(1:3),bt(1:3)))
                ut(2)=um(3)*bt(3)-um(4)*bt(2)
                ut(3)=um(4)*bt(1)-um(2)*bt(3)
                ut(4)=um(2)*bt(2)-um(3)*bt(1)
                utt(2)=ut(3)*bt(3)-ut(4)*bt(2)
                utt(3)=ut(4)*bt(1)-ut(2)*bt(3)
                utt(4)=ut(2)*bt(2)-ut(3)*bt(1)
                up(1:4)=um(1:4)+(ut(1:4)+utt(1:4))*s
                !放射減衰の第二項目の係数の計算
                rr=up(1)*(-dot_product(NE(1:3),NE(1:3))*up(1)+(NE(2)*NB(3)-NE(3)*NB(2))*up(2)+(NE(3)*NB(1)-NE(1)*NB(3))*up(3)&
                        &+(NE(1)*NB(2)-NE(2)*NB(1))*up(4))&
                        &+up(2)*((NE(2)*NB(3)-NE(3)*NB(2))*up(1)+(NE(1)**2.0d0-NB(3)**2.0d0-NB(2)**2.0d0)*up(2)&
                        &+(NE(1)*NE(2)+NB(1)*NB(2))*up(3)+(NE(1)*NE(3)+NB(1)*NB(3))*up(4))&
                        &+up(3)*((NE(3)*NB(1)-NE(1)*NB(3))*up(1)+(NE(1)*NE(2)+NB(1)*NB(2))*up(2)&
                        &+(NE(2)**2.0d0-NB(3)**2.0d0-NB(1)**2.0d0)*up(3)+(NE(2)*NE(3)+NB(2)*NB(3))*up(4))&
                        &+up(4)*((NE(1)*NB(2)-NE(2)*NB(1))*up(1)+(NE(1)*NE(3)+NB(1)*NB(3))*up(2)&
                        &+(NE(2)*NE(3)+NB(2)*NB(3))*up(3)+(NE(3)**2.0d0-NB(1)**2.0d0-NB(2)**2.0d0)*up(4))
                !空間成分の半ステップの電場の加速と放射減衰
                u(2)=up(2)+0.5d0*NE(1)*dt&
                        &+0.5d0*b*((NE(2)*NB(3)-NE(3)*NB(2))*up(1)+(NE(1)*NE(1)-NB(2)*NB(2)-NB(3)*NB(3))*up(2)&
                        &+(NE(1)*NE(2)+NB(1)*NB(2))*up(3)+(NE(1)*NE(3)+NB(1)*NB(3))*up(4)&
                        &+rr*up(2))*dt/up(1)
                u(3)=up(3)+0.5d0*NE(2)*dt&
                        &+0.5d0*b*((NE(3)*NB(1)-NE(1)*NB(3))*up(1)+(NE(1)*NE(2)+NB(1)*NB(2))*up(2)&
                        &+(NE(2)*NE(2)-NB(1)*NB(1)-NB(3)*NB(3))*up(3)+(NE(2)*NE(3)+NB(2)*NB(3))*up(4)&
                        &+rr*up(3))*dt/up(1)
                u(4)=up(4)+0.5d0*NE(3)*dt&
                        &+0.5d0*b*((NE(1)*NB(2)-NE(2)*NB(1))*up(1)+(NE(1)*NE(3)+NB(1)*NB(3))*up(2)&
                        &+(NE(2)*NE(3)+NB(2)*NB(3))*up(3)+(NE(3)*NE(3)-NB(1)*NB(1)-NB(2)*NB(2))*up(4)&
                        &+rr*up(4))*dt/up(1)
                !ローレンツ因子の書き換え
                u(1)=sqrt(1.0d0+u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
                !位置の導出
                xmat(1)=xmat(1)+dt/omega!u(1)*dt/u(1)**2.0d0
                xmat(2)=xmat(2)+u(2)*dt/u(1)
                xmat(3)=xmat(3)+u(3)*dt/u(1)
                xmat(4)=xmat(4)+u(4)*dt/u(1)
                !距離の変換
                t=xmat(1)
                x=xmat(2)*c/omega
                y=xmat(3)*c/omega
                z=xmat(4)*c/omega
                !電磁場の呼び出し
                call assignment_dipole(t,x,y,z,r,theta,rho,d1,d2,d3,d4,d5,d6,Bx,By,Bz,Ex,Ey,Ez)
                !電磁場の更新
                LB(1)=Bx
                LB(2)=By
                LB(3)=Bz
                LE(1)=Ex
                LE(2)=Ey
                LE(3)=Ez
                !プロット
                if(mod(i-1,iout) == 0)then
                        write(19+j,*) xmat(1:4),r*omega/c,sqrt(dot_product(F(1:3),F(1:3)))&
                                &,abs(K)*sqrt(V(2)*V(2)+V(3)*V(3)+V(4)*V(4)),sqrt(V(2)*V(2)+V(3)*V(3)+V(4)*V(4))&
                                &,sqrt(dot_product(VB(1:3),VB(1:3)))/sqrt(dot_product(LE(1:3),LE(1:3)))&
                                &,dot_product(LE(1:3),LB(1:3))/dot_product(LE(1:3),LE(1:3))&
                                &,sqrt(dot_product(LE(1:3),LE(1:3)))&
                                &,sqrt(dot_product(LB(1:3),LB(1:3)))/sqrt(dot_product(LE(1:3),LE(1:3)))
                end if
                if(r<1.19d6)exit
        end do
!end do
        close(10)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)
        close(20)
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)
        close(26)
        close(30)
        close(31)
        close(32)
        close(33)
        close(34)
        close(35)
        close(36)
end program main
