clear variables
clc

q_e=1.602e-19; %电子电量，单位：C
m_e=9.1e-31;  %电子质量，单位：kg
EPS0 = 8.854e-12; %真空中的介电常数，单位：F/m
MU0=4*pi()*1e-7; %真空中的磁导率，单位：H/m
k=1.38e-23; %玻尔兹曼常数，单位：J/K
c=3e8; %光速，单位：m/s
Te_average=10; %平均电子温度，单位：eV
n_e_average=5e18; %平均电子密度，单位：m^-3
lamda=sqrt(EPS0*Te_average/(n_e_average*q_e)); %德拜长度，单位：m
z_min=0; %z方向的最小值，单位：m
L_real=0.12; %离子源z方向的最大值，单位：m
dz=lamda; %设置空间间隔
dt=dz/(2*c); %设置时间间隔，为保证模拟精确度，保证dt<=dz/c，一般设为dz/(2*c)单位：s
nz=round(L_real/dz)+1; %格点数
z_max=dz*(nz-1); %模拟区域z方向的最大值
E0=1.5; %微波电场强度的幅值，单位：V/m（不确定）
B0=1.0; %微波磁感应强度的幅值，单位:T（不确定）
w=2.45; %微波的频率,单位:GHz
B_mic=zeros(nz+1,3); %微波的磁场
E_mic=zeros(nz+2,3); %微波的电场

num_step=30000; %运行步数
E_bound=zeros(num_step,3);

time=dt*[0:num_step-1]'; %运行时间
E_x_start=E0*cos(w*1e9*2*pi()*time); %Ex初始值
E_y_start=E0*sin(w*1e9*2*pi()*time); %Ey初始值
B_x_start=B0*sin(w*1e9*2*pi()*time); %Bx初始值
B_y_start=B0*cos(w*1e9*2*pi()*time); %By初始值

for ts=1:num_step
    E_mic(1,1:2)=[E_x_start(ts),E_y_start(ts)];
   % B_mic(1,1:2)=[B_x_start(ts),B_y_start(ts)];
    B_mic(1:nz+1,1)=B_mic(1:nz+1,1)+dt/dz*(E_mic(2:nz+2,2)-E_mic(1:nz+1,2));
    B_mic(1:nz+1,2)=B_mic(1:nz+1,2)-dt/dz*(E_mic(2:nz+2,1)-E_mic(1:nz+1,1));
    E_mic(2:nz+1,1)=E_mic(2:nz+1,1)-c^2*dt*((B_mic(2:nz+1,2)-B_mic(1:nz,2))/dz);
    E_mic(2:nz+1,2)=E_mic(2:nz+1,2)+c^2*dt*((B_mic(2:nz+1,1)-B_mic(1:nz,1))/dz);
    %E_mic(nz+2,1)=E_mic(nz+2,1)+(c*dt-dz)/(c*dt+dz)*(E_mic(nz+1,1)-E_mic(nz+2,1));
    %E_mic(nz+2,2)=E_mic(nz+2,2)+(c*dt-dz)/(c*dt+dz)*(E_mic(nz+1,2)-E_mic(nz+2,2));
    
    E_bound(ts,:)=E_mic(nz+1,:);
    
    if ts>2
        E_mic(nz+2,:)=E_bound(ts-2,:);
    end

    %E_mic(nz+2,:)=E_mic(nz+1,:);
    %B_mic(nz+1

%输出微波电磁场分布的图像
figure(1);
subplot(2,3,1);
plot(1:nz,E_mic(1:nz,1));
xlabel('grid');
ylabel('E_x');
subplot(2,3,2);
plot(1:nz,E_mic(1:nz,2));
xlabel('grid');
ylabel('E_y');
subplot(2,3,3);
plot(1:nz,E_mic(1:nz,3));
xlabel('grid');
ylabel('E_z');
subplot(2,3,4);
plot(1:nz,B_mic(1:nz,1));
xlabel('grid');
ylabel('B_x');
subplot(2,3,5);
plot(1:nz,B_mic(1:nz,2));
xlabel('grid');
ylabel('B_y');
subplot(2,3,6);
plot(1:nz,B_mic(1:nz,3));
xlabel('grid');
ylabel('B_z');
pause(0.0000000001);
end
