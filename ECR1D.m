%1 2/2D ECR plasma simulation using PIC-MCC
clear variables
clc

global q_e m_e

%定义基本参量
q_e=1.602e-19; %电子电量，单位：C
AMU = 1.661e-27; %质子质量，单位：kg
m_e=9.1e-31;  %电子质量，单位：kg
m_Li=7*AMU; %Li-7质量，单位：kg
EPS0 = 8.854e-12; %真空中的介电常数，单位：F/m
MU0=4*pi()*1e-7; %真空中的磁导率，单位：H/m
k=1.38e-23; %玻尔兹曼常数，单位：J/K
c=3e8; %光速，单位：m/s
%c=3e2; %将光速改小，以观察粒子对电磁波的影响（结果显示这样改有问题）
Te_average=10; %平均电子温度，单位：eV
n_e_average=5e17; %平均电子密度，单位：m^-3
lamda=sqrt(EPS0*Te_average/(n_e_average*q_e)); %德拜长度，单位：m
z_min=0; %z方向的最小值，单位：m
L_real=0.12; %离子源z方向的最大值，单位：m
dz=lamda; %设置空间间隔
dt=dz/(2*c); %设置时间间隔，为保证模拟精确度，保证dt<=dz/c，一般设为dz/(2*c)单位：s
nz=round(L_real/dz)+1; %格点数
z_max=dz*(nz-1); %模拟区域z方向的最大值
spwt=1e10; %每个宏粒子包含50个真实粒子
vth_e=sqrt(2*q_e*0.5/m_e); %电子的热速度
N_e=1000; %设置500个电子用于程序测试
tol=0.1; %电势求解的精度值
max_part=1000; %粒子容器的最大值
B_ex_z=0.0875; %外部永磁铁产生的磁感应强度，单位：T
E0=1.5e5; %微波电场强度的幅值，单位：V/m（不确定）
B0=0.0001; %微波磁感应强度的幅值，单位:T（不确定）
w=2.45; %微波的频率,单位:GHz
step_num=10000; %总共运行的时间步数

%预分配空间
phi=zeros(nz,1); %电势
phi_tem=zeros(nz,1); %电势求解的中间变量
B_mic=zeros(nz,3); %微波的磁场
E_mic=zeros(nz,3); %微波的电场
E_sta=zeros(nz,1); %电子与离子自身产生的静电场
E_e=zeros(max_part,3); %粒子位置处的电场强度
B_e=zeros(max_part,3); %粒子位置处的磁感应强度
E_boun=zeros(step_num,3); %用于设置吸收边界
B_boun=zeros(step_num,3); %用于设置吸收边界
%den=zeros(nz,1); %每个格点上的粒子数密度
%chg=zeros(nz,1); %每个格点上的电荷比例
J=zeros(nz,3); %每个格点上的电流密度
pos=zeros(max_part,1); %粒子的位置信息
vel=zeros(max_part,3); %粒子的速度信息



%给定初始电子的位置和速度
pos(1:N_e)=rand(N_e,1)*z_max; %初始电子的位置在0到z_max之间

%vel(1:N_e/2,:)=vth_e*rand(N_e/2,3);
%vel(N_e/2+1:N_e,:)=-vth_e*rand(N_e/2,3);

%vel(1:N_e,:)=1e8;
vel(1:N_e,:)=vth_e*2*(rand(N_e,3)+rand(N_e,3)+rand(N_e,3)-1.5); %设置初始电子的速度在热速度的-3倍到3倍之间
vel(1:N_e,:)=UpdateVelocity(E_e(1:N_e,:),B_e(1:N_e,:),vel(1:N_e,:),-0.5*dt); %将电子的速度向前移动半个时间步长

%开始主循环，即对时间的循环
for ts=1:step_num %运行的时间步数
    
    chg=zeros(nz,1);
    den=zeros(nz,1); 
    den_vel=zeros(nz,3); %用于存放每个格点上的电流密度分量
    B_mic_old=B_mic; %B_mic_old用于和B_mic求平均，从而用于求解速度
    pos_half=zeros(N_e,1);
    
    for p_half=1:N_e
        pos_half(p_half,1)=pos(p_half,1)-vel(p_half,3)*0.5*dt;
        if pos_half(p_half,1)<0
            pos_half(p_half,1)=pos_half(p_half,1)+z_max;
        end
        if pos_half(p_half,1)>z_max
            pos_half(p_half,1)=pos_half(p_half,1)-z_max;
        end
    end
   % pos_half=pos-vel(:,3)*0.5*dt; %将电子的位置向前移动半个时间步长，以求得半时刻处的电流密度
    
    %首先将所有电子的电荷按照最近格点原则分配
    for p=1:N_e
        fi=1+pos(p)/dz; %实际格点位置，为浮点数
        i=floor(fi);  %对应整数格点位置
        hz=fi-i;      %粒子与第i个格点之间的距离占距离步长的比例
        
        %将电子电荷按照比例分配到临近的两个格点上
        chg(i)=chg(i)+(1-hz); %分配到第i个格点上的电荷比例
        chg(i+1)=chg(i+1)+hz; %分配到第i+1个格点上的电荷比例
        
        fi_half=1+pos_half(p)/dz; %半时刻处的格点位置
        i_half=floor(fi_half); %半时刻处格点位置对应的整数格点位置
        hz_half=fi_half-i_half; %半时刻处粒子与第i_half个格点之间的距离占距离步长的比例
        
        
        
        %将速度按比例分配到临近的两个格点上
        den_vel(i_half,:)=den_vel(i_half,:)+spwt*q_e*vel(p,:)*(1-hz_half); %分配到第i个格点上的速度
        den_vel(i_half+1,:)=den_vel(i_half+1,:)+spwt*q_e*vel(p,:)*hz_half; %分配到第i+1个格点上的速度
    end
    
    den=spwt*q_e*chg/dz; %每个格点上的电荷密度
    %J=spwt*q_e*(-1)*den_vel; %每个格点上的电流密度
    %J=den_vel.*den;
    J=den_vel;
    
    %应用边界条件，在边界上仅有一半的长度有贡献，因此边界上的密度应当×2
    den(1)=2*den(1); %左边界
    den(nz)=2*den(nz); %右边界
    
    %开始求解泊松方程，得到格点上的电势，所用方法为高斯-赛德尔迭代法（G-S）
    % phi_tem=phi;
    for j=1:2000 %迭代步数
        
        for i=2:nz-1 %中间格点的电势
            phi_tem(i)=0.5*(den(i)/EPS0*dz*dz+phi(i-1)+phi(i+1)); %根据有限差分方法得到
        end
        
        phi_tem(1)=phi(2); %入口处为第一类边界条件
        phi_tem(nz)=phi(nz-1); %出口处也为第一类边界条件
        
        if mod(j,10)==0
            R=norm(phi_tem-phi); %求出本次迭代的残差
            if (R<=tol)   %收敛条件
                phi=phi_tem;
                %fprintf('Congratulations! The solver of the electric potential is converged at %d step\n',j);
                break;
            end
        end
        phi=phi_tem;
    end
    if j>=2000 %如果电势求解不收敛，则终止程序的运行
        fprintf('The result of the electric potential is not convergent, the program is stopped!\n');
        break;
    end
    
    %求解每个格点上的静电场，由于是一维，因此只能得到z方向的电场强度
    E_sta(2:nz-1)=(phi(1:nz-2)-phi(3:nz))/2/dz; %中间格点的静电场
    E_sta(1)=(phi(1)-phi(2))/dz;  %入口处静电场
    E_sta(nz)=(phi(nz-1)-phi(nz))/dz; %出口处静电场
    
    %计算微波电磁场的各个分量
    E_mic_xold=E0*cos(w*1e9*2*pi()*(ts-2)*dt);
    E_mic_yold=E0*sin(w*1e9*2*pi()*(ts-2)*dt);
    E_mic(1,1)=E0*cos(w*1e9*2*pi()*(ts-1)*dt); %入口处x方向的电场强度E_x=E0cos(wt)
   % B_mic_xin_full=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt); %入口处x方向磁感应强度，时间为整
    B_mic_xin=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt); %入口处x方向的磁感应强度为B_x=B0sin(wt)，且时间上向后移动Δt/2
    E_mic(1,2)=E0*sin(w*1e9*2*pi()*(ts-1)*dt); %入口处y方向的电场强度E_y=E0sin(wt)
    %B_mic_yin_full=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt); %入口处y方向磁感应强度，时间为整
    B_mic_yin=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt); %入口处y方向的磁感应强度B_y=B0cos(wt)，且时间上向后移动Δt/2
    %B_mic(1,1)=B_mic_xin+dt/(2*dz)*(E_mic(2,2)-E_mic(1,2)); %将Bx在空间和时间上都向后移动Δz/2和Δt/2
    B_mic(1,1)=B_mic_xin+dz/2*(1/(c^2)*(E_mic(1,2)-E_mic_yold)/dt+MU0*J(1,2)); %入口处的Bx，将Bx在空间上向后移动Δz/2，向后差分
    B_mic(1,2)=B_mic_yin-dz/2*(1/(c^2)*(E_mic(1,1)-E_mic_xold)/dt+MU0*J(1,1)); %入口处的By，将By在空间上向后移动Δz/2，向后差分
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
    E_mic(2:nz-1,1)=E_mic(2:nz-1,1)-c^2*dt*(MU0*J(2:nz-1,1)+(B_mic(2:nz-1,2)-B_mic(1:nz-2,2))/dz); %Ex
    E_mic(2:nz-1,2)=E_mic(2:nz-1,2)+c^2*dt*(-1*MU0*J(2:nz-1,2)+(B_mic(2:nz-1,1)-B_mic(1:nz-2,1))/dz); %Ey
    %E_mic(1:nz,3)=E_mic(1:nz,3)-c^2*MU0*dt*J(1:nz,3); %Ez
    
    E_boun(ts,:)=E_mic(nz-1,:);
    B_boun(ts,:)=B_mic(nz-1,:);
    if ts>2
        E_mic(nz,:)=E_boun(ts-2,:);
        B_mic(nz,:)=B_boun(ts-2,:);
    end
    
    %将格点上的静电场和电磁场全部分配到每个粒子上，另外还要加上外部永磁体产生的0.0875T的磁场
    E_e=zeros(max_part,3); %每个电子所在位置处的电场强度
    B_e=zeros(max_part,3); %每个电子所在位置处的磁感应强度
    B_e(:,3)=B_e(:,3)+B_ex_z; %z方向的磁感应强度为0.0875T
    B_mic_ave=0.5*(B_mic_old+B_mic); %利用平均，将两个半时间的值平均到到整时间上
    p=1; %用于遍历所有粒子
    while p<=N_e
        fi=1+pos(p)/dz; %实际格点位置，为浮点数
        i=floor(fi);  %对应整数格点位置
        hz=fi-i;      %粒子与第i个格点之间的距离占距离步长的比例
        E_e(p,3)=E_e(p,3)+E_sta(i)*(1-hz)+E_sta(i+1)*hz; %分配粒子之间的静电场
        E_e(p,:)=E_e(p,:)+E_mic(i,:)*(1-hz)+E_mic(i+1,:)*hz; %分配微波产生的电场
        if fi<=1.5    %四个if用来确定粒子位置，从而为之分配磁感应强度
            B_e(p,:)=B_e(p,:)+2*hz*B_mic_ave(1,:)+(1-2*hz)*[B_mic_xin,B_mic_yin,0];
        end
        if fi>=nz-0.5
            B_e(p,:)=B_e(p,:)+2*(1-hz)*B_mic_ave(nz-1,:)+(2*hz-1)*[B_mic_xout,B_mic_yout,0];
        end
        if fi>1.5 && fi<nz-0.5 && hz<=0.5
            B_e(p,:)=B_e(p,:)+(hz+0.5)*B_mic_ave(i,:)+(0.5-hz)*B_mic_ave(i+1,:);
        end
        if fi>1.5 && fi<nz-0.5 && hz>0.5
            B_e(p,:)=B_e(p,:)+(1.5-hz)*B_mic_ave(i,:)+(hz-0.5)*B_mic_ave(i+1,:);
        end
        
        %利用BORIS方法求出下一时刻的速度和位置
        vel(p,:)=UpdateVelocity(E_e(p,:),B_e(p,:),vel(p,:),dt); %更新速度
        pos(p,1)=pos(p,1)+vel(p,3)*dt;   %更新位置
        
        %设置周期性边界条件
        if pos(p,1)<0
            pos(p,1)=pos(p,1)+z_max;
        end
        if pos(p,1)>z_max
            pos(p,1)=pos(p,1)-z_max;
        end
        
        p=p+1;
    end
   % plot(1:nz,E_mic(:,1));
   % pause(0.1);
end

fprintf('finished!\n');

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

%输出v-z图像
figure(2);
subplot(3,1,1);
scatter(pos(1:N_e),vel(1:N_e,1));
xlabel('z(m)');
ylabel('v_x(m/s)');
subplot(3,1,2);
scatter(pos(1:N_e),vel(1:N_e,2));
xlabel('z(m)');
ylabel('v_y(m/s)');
subplot(3,1,3);
scatter(pos(1:N_e),vel(1:N_e,3));
xlabel('z(m)');
ylabel('v_z(m/s)');


