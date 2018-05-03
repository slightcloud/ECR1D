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
L_real=0.05; %离子源z方向的最大值，单位：m
dz=lamda; %设置空间间隔
dt=dz/(2*c); %设置时间间隔，为保证模拟精确度，保证dt<=dz/c，一般设为dz/(2*c)单位：s
nz=round(L_real/dz)+1; %格点数
z_max=dz*(nz-1); %模拟区域z方向的最大值
E0=1.5; %微波电场强度的幅值，单位：V/m（不确定）
B0=0.0001; %微波磁感应强度的幅值，单位:T（不确定）
w=2.45; %微波的频率,单位:GHz
step_num=10000; %总共运行的时间步数
B_mic=zeros(nz,3); %微波的磁场
E_mic=zeros(nz,3); %微波的电场
E_boun=zeros(step_num,3); %用于设置吸收边界
B_boun=zeros(step_num,3); %用于设置吸收边界

%计算微波电磁场的各个分量
for ts=1:step_num
    E_mic_xold=E0*cos(w*1e9*2*pi()*(ts-2)*dt);
    E_mic_yold=E0*sin(w*1e9*2*pi()*(ts-2)*dt);
    E_mic(1,1)=E0*cos(w*1e9*2*pi()*(ts-1)*dt); %入口处x方向的电场强度E_x=E0cos(wt)
    % B_mic_xin_full=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt); %入口处x方向磁感应强度，时间为整
    B_mic_xin=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt); %入口处x方向的磁感应强度为B_x=B0sin(wt)，且时间上向后移动Δt/2
    E_mic(1,2)=E0*sin(w*1e9*2*pi()*(ts-1)*dt); %入口处y方向的电场强度E_y=E0sin(wt)
    %B_mic_yin_full=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt); %入口处y方向磁感应强度，时间为整
    B_mic_yin=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt); %入口处y方向的磁感应强度B_y=B0cos(wt)，且时间上向后移动Δt/2
    %B_mic(1,1)=B_mic_xin+dt/(2*dz)*(E_mic(2,2)-E_mic(1,2)); %将Bx在空间和时间上都向后移动Δz/2和Δt/2
    B_mic(1,1)=B_mic_xin+dz/2*(1/(c^2)*(E_mic(1,2)-E_mic_yold)/dt); %入口处的Bx，将Bx在空间上向后移动Δz/2，向后差分
    B_mic(1,2)=B_mic_yin-dz/2*(1/(c^2)*(E_mic(1,1)-E_mic_xold)/dt); %入口处的By，将By在空间上向后移动Δz/2，向后差分
    %B_mic(1,2)=B_mic_yin-dt/(2*dz)*(E_mic(2,1)-E_mic(1,1)); %将By在空间和时间上都向后移动Δz/2和Δt/2
    
    % E_mic(nz,1)=0; %出口处x方向的电场强度0
    % E_mic(nz,2)=0; %出口处y方向的电场强度0
    
    B_mic_xout=B_mic(nz,1);
    B_mic_yout=B_mic(nz,2);
    
    %B_mic_xout=0; %出口处x方向的磁感应强度0
    % B_mic_yout=0; %出口处y方向的磁感应强度0
    % B_mic(nz-1,1)=B_mic_xout-dz/2*MU0*J(nz,2); %出口处的Bx，向前差分
    % B_mic(nz-1,2)=B_mic_yout+dz/2*MU0*J(nz,1); %出口处的By，向前差分
    
    B_mic(2:nz-1,1)=B_mic(2:nz-1,1)+dt/dz*(E_mic(3:nz,2)-E_mic(2:nz-1,2)); %Bx
    B_mic(2:nz-1,2)=B_mic(2:nz-1,2)-dt/dz*(E_mic(3:nz,1)-E_mic(2:nz-1,1)); %By
    E_mic(2:nz-1,1)=E_mic(2:nz-1,1)-c^2*dt*((B_mic(2:nz-1,2)-B_mic(1:nz-2,2))/dz); %Ex
    E_mic(2:nz-1,2)=E_mic(2:nz-1,2)+c^2*dt*((B_mic(2:nz-1,1)-B_mic(1:nz-2,1))/dz); %Ey
    %E_mic(1:nz,3)=E_mic(1:nz,3)-c^2*MU0*dt*J(1:nz,3); %Ez
    
    E_boun(ts,:)=E_mic(nz-1,:);
    B_boun(ts,:)=B_mic(nz-1,:);
    if ts>2
        E_mic(nz,:)=E_boun(ts-2,:);
        B_mic(nz,:)=B_boun(ts-2,:);
    end
    
end
%输出微波电磁场分布的图像
figure(1);
subplot(2,3,1);
plot(1:nz,E_mic(:,1));
xlabel('grid');
ylabel('E_x');
subplot(2,3,2);
plot(1:nz,E_mic(:,2));
xlabel('grid');
ylabel('E_y');
subplot(2,3,3);
plot(1:nz,E_mic(:,3));
xlabel('grid');
ylabel('E_z');
subplot(2,3,4);
plot(1:nz,B_mic(:,1));
xlabel('grid');
ylabel('B_x');
subplot(2,3,5);
plot(1:nz,B_mic(:,2));
xlabel('grid');
ylabel('B_y');
subplot(2,3,6);
plot(1:nz,B_mic(:,3));
xlabel('grid');
ylabel('B_z');
%pause(0.0000000000000001);

