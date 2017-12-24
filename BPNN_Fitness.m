function tol_fitness=BPNN_Fitness(x)
ts=0.01; 

IN=4; H=5; Out=3;
%离线训练样本数（离线迭代训练次数）
train_num=100; 
%从粒子向量中恢复出神经网络的权值参数矩阵  
for t=1:H
       wi(t,:)=x(1,(t-1)*IN+1:t*IN);
end
for r=1:Out
       wo(r,:)=x(1, ( (IN*H+1)+(r-1)*H ): ( (IN*H+1)+r*H-1 ) );
end

%M=[10,1,10];
M=[9.9,9.8,9.4];
%M=[1.8300,0.6629,0.6088];

u_1=0; u_2=0; u_3=0; u_4=0; u_5=0;u_6=0;u_7=0;
y_1=0; y_2=0; y_3=0; 
error_1=0; error_2=0;
  
Oh=zeros(H,1);

sys=tf(400,[1,50,0]);
dsys=c2d(sys,ts,'z');
[num,den]=tfdata(dsys,'v');

for k=1:1:train_num
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

    %%%只进行前向传播---------------------------------------
    net2=xi*(wi');
    for j=1:1:H
        Oh(j)=( exp( net2(j)-exp(-net2(j)) ) )/(exp( net2(j)+exp(-net2(j)) ));
    end
    net3=wo*Oh;
    for l=1:1:Out
        K(l)=exp(net3(l))/(exp(net3(l))+exp(-net3(l)));
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
    u_7=u_6;u_6=u_5;u_5=u_4; u_4=u_3;u_3=u_2;u_2=u_1;u_1=u(k);   
    y_2=y_1; y_1=yout(k); 
    error_2=error_1; error_1=error(k);
end

%统计粒子总的fitness值
tol_fitness=0;
for k=1:train_num
        tol_fitness=tol_fitness+abs(error(k)^2);
end
tol_fitness=tol_fitness/train_num;
end













