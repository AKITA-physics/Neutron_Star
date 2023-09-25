!新バージョン(成分分解)
implicit none
integer i,k,l,imax,iout
real(8) t,gamma0
real(8) rr,s
real(8) b,dt,c,e,me
real(8) NE(1:3),NB(1:3),LB(1:3),LE(1:3),u(1:4),x(1:4)
real(8) bt(1:3),um(1:4),ut(1:4),utt(1:4),up(1:4)
real(8) plot(1000000000,2)
data    x/0.0d0,0.0d0,0.0d0,0.0d0/

open(20,file='es00.d')
open(21,file='es05.d')
open(22,file='es010.d')
open(23,file='es015.d')
open(30,file='es40.d')
open(31,file='msp45.d')
open(32,file='es410.d')
open(33,file='es415.d')
open(40,file='es80.d')
open(41,file='es85.d')
open(42,file='es810.d')
open(43,file='es815.d')

c=2.99792458d0*1.0d8
e=1.60217664*1.0d-19
me=9.1093837*1.0d-31

!パラメーター
k=2
gamma0=10.0d0**(dble(k-1)*4.0d0)
l=2
b=10.0d0**(-5.0d0*(dble(l)-1.0d0))
dt=5.0d-4
imax=10**7
iout=5*10**4

!速度、電磁場の更新
u(1)=gamma0
u(2)=0.0d0
u(3)=sqrt(u(1)*u(1)-1.0d0)
u(4)=0.0d0
LB(1)=0.0d0
LB(2)=0.0d0
LB(3)=1.5d0*(me*c**2.0d0)**2.0d0*b*1.0d4/e**3.0d0
LE(1)=0.0d0
LE(2)=0.0d0
LE(3)=0.0d0

NE(1:3)=LE(1:3)/sqrt(dot_product(LB(1:3),LB(1:3)))
NB(1:3)=LB(1:3)/sqrt(dot_product(LB(1:3),LB(1:3)))


t=0.0d0
do i = 1,imax
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
        u(2)=up(2)+0.5*NE(1)*dt&
                &+0.5*b*((NE(2)*NB(3)-NE(3)*NB(2))*up(1)+(NE(1)*NE(1)-NB(2)*NB(2)-NB(3)*NB(3))*up(2)&
                &+(NE(1)*NE(2)+NB(1)*NB(2))*up(3)+(NE(1)*NE(3)+NB(1)*NB(3))*up(4)&
                &+rr*up(2))*dt/up(1)
        u(3)=up(3)+0.5*NE(2)*dt&
                &+0.5*b*((NE(3)*NB(1)-NE(1)*NB(3))*up(1)+(NE(1)*NE(2)+NB(1)*NB(2))*up(2)&
                &+(NE(2)*NE(2)-NB(1)*NB(1)-NB(3)*NB(3))*up(3)+(NE(2)*NE(3)+NB(2)*NB(3))*up(4)&
                &+rr*up(3))*dt/up(1)
        u(4)=up(4)+0.5*NE(3)*dt&
                &+0.5*b*((NE(1)*NB(2)-NE(2)*NB(1))*up(1)+(NE(1)*NE(3)+NB(1)*NB(3))*up(2)&
                &+(NE(2)*NE(3)+NB(2)*NB(3))*up(3)+(NE(3)*NE(3)-NB(1)*NB(1)-NB(2)*NB(2))*up(4)&
                &+rr*up(4))*dt/up(1)
        !ローレンツ因子の書き換え
        u(1)=sqrt(1.0d0+u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
        !位置の導出
        x(1)=x(1)+u(1)*dt/u(1)**2.0d0
                x(2)=x(2)+u(2)*dt/(u(1)*gamma0)
                x(3)=x(3)+u(3)*dt/(u(1)*gamma0)
                x(4)=x(4)+u(4)*dt/(u(1)*gamma0)
        !時間の更新
        t=t+dt
        !プロット
        if(mod(i-1,iout) == 0)then
                plot(i,1)=x(1)
                plot(i,2)=u(1)
                write(31,*) x(1:4)!(plot(i,j),j=1,2)
        end if
end do
close(20)
close(21)
close(22)
close(23)
close(30)
close(31)
close(32)
close(33)
close(40)
close(41)
close(42)
close(43)
end

