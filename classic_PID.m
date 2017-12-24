function [time,yout]=classic_PID()

ts=0.01;
sys=tf(400,[1,50,0]);
dsys=c2d(sys,ts,'z');
[num,den]=tfdata(dsys,'v');

u_1=0; u_2=0; u_3=0; u_4=0; u_5=0;
y_1=0; y_2=0; y_3=0;
error_1=0; error_2=0;
x=[0,0,0]';

for k = 1:1:200
    time(k)=k*ts;
    rin(k)=1.0;
%无超调参数
%     kp=8;
%     ki=0.10;
%     kd=10;
%有超调与震荡参数
    kp=4.0795;
    ki=0.8703;
    kd=12.0019;
   
    du(k)=kp*x(1)+kd*x(2)+ki*x(3);
    u(k)=u_1+du(k);
    if(u(k)>10)
        u(k)=10;
    end
    if(u(k)<-10)
        u(k)=-10;
    end
    %yout(k)=-den(2)*y_1+num(2)*u_5;
    yout(k)=-den(2)*y_1-den(3)*y_2+num(2)*u_1+num(3)*u_2;
    error(k)=rin(k)-yout(k);
    u_5=u_4; u_4=u_3; u_3=u_2; u_2=u_1; u_1=u(k);
    y_3=y_2; y_2=y_1; y_1=yout(k);
    x(1)=error(k)-error_1;
    x(2)=error(k)-2*error_1+error_2;
    x(3)=error(k);
    error_2=error_1;
    error_1=error(k);
end
% plot(time,rin,'r',time,yout,'b');
% xlabel('time(s)');
% ylabel('rin,yout');