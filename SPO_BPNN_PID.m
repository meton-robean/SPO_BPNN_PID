
clear all;
close all;

xite=0.50;
alfa=0.05; 
IN=4; H=5; Out=3;

ts=0.01;

% wi=0.50*rands(H,IN);
% wo=0.50*rands(Out,H);

%结合粒子群数值优化算法的神经网络权重初始化矩阵
wi=[
   -4.7730    5.0000    4.8238   -4.7085
    4.9337    4.8659   -4.9651   -4.8618
    4.9470    4.9485   -4.6965    4.8327
    5.0000    4.9508    4.5696    4.9154
    4.6047   -4.7804    5.0000   -4.7332
   ];
wo=[
    4.9438   -4.6610    4.9293    4.7416   -4.9323
    5.0000    4.9414    4.0693   -4.7460    5.0000
   -4.9147    4.8119    4.9292    4.9434    4.9398
   ];

%[wi , wo]=SPO_InitW();
wi_init_save=wi;   wo_init_save=wo;   
wo_1=wo;  wo_2=wo; 
wi_1=wi;  wi_2=wi;
 
%M=[10,1,10];
%激活函数系数
M=[9.9,9.8,9.4];

x=[0,0,0];
du_1=0;
u_1=0; u_2=0; u_3=0; u_4=0; u_5=0;u_6=0;u_7=0;
y_1=0; y_2=0; y_3=0; 
error_1=0; error_2=0;
  
Oh=zeros(H,1);
I=Oh;

  
sys=tf(400,[1,50,0]);
dsys=c2d(sys,ts,'z');
[num,den]=tfdata(dsys,'v');

for k=1:1:200 
    time(k)=k*ts;
    rin(k)=1.0; 
    yout(k)=-den(2)*y_1-den(3)*y_2+num(2)*u_1+num(3)*u_2;
    error(k)=rin(k)-yout(k);
    X(1)=error(k)-error_1;
    X(2)=error(k);
    X(3)=error(k)-2*error_1+error_2;
    xii=[X(1),X(2),X(3),1];
    xi=xii/norm(xii);
    epid=[X(1);X(2);X(3)];

    %%%前向传播----------------------------------------
    
    net2=xi*(wi');
     
    for j=1:1:H
        Oh(j)=( exp( net2(j)-exp(-net2(j)) ) )/(exp( net2(j)+exp(-net2(j)) ));
    end
   
    net3=wo*Oh;
      
    for l=1:1:Out
        K(l)=exp(net3(l))/(exp(net3(l))+exp(-net3(l)));
        %K(l)=M*net3(l);
    end
    kp(k)=M(1)*K(1); ki(k)=M(2)*K(2); kd(k)=M(3)*K(3);
    Kpid=[kp(k),ki(k),kd(k)];
    du(k)=Kpid*epid;
    u(k)=u_1+du(k); 
    if u(k)>10
        u(k)=10;
    end
    if u(k)<-10
        u(k)=-10;
    end
    
    %%%后向传播------------------------------------------------
    dyu(k)=sign((yout(k)-y_1)/(du(k)-du_1+0.0001));
    for j=1:1:Out
            dK(j)=1/(exp(net3(j))+exp(-net3(j)));
           %dK(j)=M;
    end
      
    for l=1:1:Out
        delta3(l)=error(k)*dyu(k)*epid(l)*dK(l);
    end

    for l=1:1:Out
        for i=1:1:H
            d_wo=xite*delta3(l)*Oh(i)+alfa*(wo_1-wo_2);
        end
    end
      
    wo=wo_1+d_wo+alfa*(wo_1-wo_2);
    %更新M系数
%     for h = 1:Out
%         dM(h)=error(k)*K(h);
%         M(h)=M(h)+rite*dM(h);
%     end
    
    for i=1:1:H
        dO(i)=4/(exp(net2(i))+exp(-net2(i)))^2;
    end
    segma=delta3*wo;
    for i=1:1:H
        delta2(i)=dO(i)*segma(i);
    end
    d_wi=xite*delta2'*xi;
      
    wi=wi_1+d_wi+alfa*(wi_1-wi_2);
      
    wo_2=wo_1; wo_1=wo;
    wi_2=wi_1; wi_1=wi;
    du_1=du(k);
    
    u_7=u_6;u_6=u_5;u_5=u_4; u_4=u_3;u_3=u_2;u_2=u_1;u_1=u(k);   
    y_2=y_1; y_1=yout(k); 
    error_2=error_1; error_1=error(k);

end
figure(2);
[t,y]=BPNN_PID();
plot(t,y,'g','Linewidth',2);
hold on;

plot(time,rin,'r','Linewidth',2);
xlabel('t/s');  ylabel('rin,yout');
hold on ;

plot(time,yout,'c','Linewidth',2);
hold on;

[x,y]=classic_PID();
plot(x,y,'b','Linewidth',2);
%-----------------

figure(3);
plot(time,error,'r','Linewidth',2);
xlabel('t/s');  ylabel('error');

figure(4);
plot(time,u,'r','Linewidth',2);
xlabel('t/s');  ylabel('u');

figure(5);
subplot(311);
plot(time,kp,'r','Linewidth',2);
xlabel('t/s');  ylabel('kp');
subplot(312);
plot(time,ki,'g','Linewidth',2);
xlabel('t/s');  ylabel('ki');
subplot(313);
plot(time,kd,'b','Linewidth',2);
xlabel('t/s'); ylabel('kd');













