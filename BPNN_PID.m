function [time,yout]=BPNN_PID()
%%惯性系数
xite=0.50;
alfa=0.05; 
  %神经网络结构参数
IN=4; H=5; Out=3;
  %采样时间
ts=0.01;
  %初始化权值参数
% wi=0.50*rands(H,IN);
% wo=0.50*rands(Out,H);
 wi=[
    -0.3234    0.3269   -0.3827   -0.0136
    0.4135   -0.2866   -0.3405   -0.4651
    0.2692    0.2864   -0.4887    0.1367
   -0.1759    0.1077   -0.1479    0.2126
    0.3873   -0.0634    0.3556    0.1271
     ];
 
wo=[
    0.0208    0.3921    0.0192   -0.4522   -0.0991
   -0.4179    0.0088   -0.1617   -0.4256   -0.3543
   -0.3415    0.1173    0.1547    0.4085    0.3075
           ];
wi_init_save=wi;   wo_init_save=wo;   
wo_1=wo;  wo_2=wo; 
wi_1=wi;  wi_2=wi;

M=[10,1,10];

x=[0,0,0];
du_1=0;
u_1=0; u_2=0; u_3=0; u_4=0; u_5=0;u_6=0;u_7=0;
y_1=0; y_2=0; y_3=0; 
error_1=0; error_2=0;
Oh=zeros(H,1);
I=Oh;

  %被控对象传递函数的z变换
sys=tf(400,[1,50,0]);
dsys=c2d(sys,ts,'z');
[num,den]=tfdata(dsys,'v');
%%开始采样
for k=1:1:200 
    time(k)=k*ts;
    rin(k)=1.0;   %阶跃信号
    yout(k)=-den(2)*y_1-den(3)*y_2+num(2)*u_1+num(3)*u_2;
    error(k)=rin(k)-yout(k); %误差
    X(1)=error(k)-error_1;
    X(2)=error(k);
    X(3)=error(k)-2*error_1+error_2;
    xii=[X(1),X(2),X(3),1];
    xi=xii/norm(xii);
    epid=[X(1);X(2);X(3)];

    %%%前向传播----------------------------------------
     %输入层到隐层线性变换
    net2=xi*(wi');
      %tanh激活函数
    for j=1:1:H
        Oh(j)=( exp( net2(j)-exp(-net2(j)) ) )/(exp( net2(j)+exp(-net2(j)) ));
    end
    %隐层到输出层线性变换
    net3=wo*Oh;
      %sigmoid激活函数，输出
    for l=1:1:Out
        K(l)=exp(net3(l))/(exp(net3(l))+exp(-net3(l)));
    end
    kp(k)=M(1)*K(1); ki(k)=M(2)*K(2); kd(k)=M(3)*K(3);
    Kpid=[kp(k),ki(k),kd(k)];
    du(k)=Kpid*epid;
    u(k)=u_1+du(k); %更新控制器输出量
    if u(k)>10
        u(k)=10;
    end
    if u(k)<-10
        u(k)=-10;
    end
    %%%反向传播更新权值------------------------------------------------
    dyu(k)=sign((yout(k)-y_1)/(du(k)-du_1+0.0001));
    for j=1:1:Out
        dK(j)=1/(exp(net3(j))+exp(-net3(j)));
    end
      %求局部梯度因子delta3
    for l=1:1:Out
        delta3(l)=error(k)*dyu(k)*epid(l)*dK(l);
    end

    for l=1:1:Out
        for i=1:1:H
            d_wo=xite*delta3(l)*Oh(i)+alfa*(wo_1-wo_2);
        end
    end
      %更新wo
    wo=wo_1+d_wo+alfa*(wo_1-wo_2);
    for i=1:1:H
        dO(i)=4/(exp(net2(i))+exp(-net2(i)))^2;
    end
    segma=delta3*wo;
    for i=1:1:H
        delta2(i)=dO(i)*segma(i);
    end
    d_wi=xite*delta2'*xi;
      %更新wi
    wi=wi_1+d_wi+alfa*(wi_1-wi_2);
    
    wo_2=wo_1; wo_1=wo;
    wi_2=wi_1; wi_1=wi;
    du_1=du(k);   
    u_7=u_6;u_6=u_5;u_5=u_4; u_4=u_3;u_3=u_2;u_2=u_1;u_1=u(k);   
    y_2=y_1; y_1=yout(k); 
    error_2=error_1; error_1=error(k);

end













