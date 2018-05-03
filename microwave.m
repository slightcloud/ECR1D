clear variables
clc

q_e=1.602e-19; %���ӵ�������λ��C
m_e=9.1e-31;  %������������λ��kg
EPS0 = 8.854e-12; %����еĽ�糣������λ��F/m
MU0=4*pi()*1e-7; %����еĴŵ��ʣ���λ��H/m
k=1.38e-23; %����������������λ��J/K
c=3e8; %���٣���λ��m/s
Te_average=10; %ƽ�������¶ȣ���λ��eV
n_e_average=5e18; %ƽ�������ܶȣ���λ��m^-3
lamda=sqrt(EPS0*Te_average/(n_e_average*q_e)); %�°ݳ��ȣ���λ��m
z_min=0; %z�������Сֵ����λ��m
L_real=0.05; %����Դz��������ֵ����λ��m
dz=lamda; %���ÿռ���
dt=dz/(2*c); %����ʱ������Ϊ��֤ģ�⾫ȷ�ȣ���֤dt<=dz/c��һ����Ϊdz/(2*c)��λ��s
nz=round(L_real/dz)+1; %�����
z_max=dz*(nz-1); %ģ������z��������ֵ
E0=1.5; %΢���糡ǿ�ȵķ�ֵ����λ��V/m����ȷ����
B0=0.0001; %΢���Ÿ�Ӧǿ�ȵķ�ֵ����λ:T����ȷ����
w=2.45; %΢����Ƶ��,��λ:GHz
step_num=10000; %�ܹ����е�ʱ�䲽��
B_mic=zeros(nz,3); %΢���Ĵų�
E_mic=zeros(nz,3); %΢���ĵ糡
E_boun=zeros(step_num,3); %�����������ձ߽�
B_boun=zeros(step_num,3); %�����������ձ߽�

%����΢����ų��ĸ�������
for ts=1:step_num
    E_mic_xold=E0*cos(w*1e9*2*pi()*(ts-2)*dt);
    E_mic_yold=E0*sin(w*1e9*2*pi()*(ts-2)*dt);
    E_mic(1,1)=E0*cos(w*1e9*2*pi()*(ts-1)*dt); %��ڴ�x����ĵ糡ǿ��E_x=E0cos(wt)
    % B_mic_xin_full=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt); %��ڴ�x����Ÿ�Ӧǿ�ȣ�ʱ��Ϊ��
    B_mic_xin=B0*sin(w*1e9*2*pi()*(ts-0.5)*dt); %��ڴ�x����ĴŸ�Ӧǿ��ΪB_x=B0sin(wt)����ʱ��������ƶ���t/2
    E_mic(1,2)=E0*sin(w*1e9*2*pi()*(ts-1)*dt); %��ڴ�y����ĵ糡ǿ��E_y=E0sin(wt)
    %B_mic_yin_full=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt); %��ڴ�y����Ÿ�Ӧǿ�ȣ�ʱ��Ϊ��
    B_mic_yin=B0*cos(w*1e9*2*pi()*(ts-0.5)*dt); %��ڴ�y����ĴŸ�Ӧǿ��B_y=B0cos(wt)����ʱ��������ƶ���t/2
    %B_mic(1,1)=B_mic_xin+dt/(2*dz)*(E_mic(2,2)-E_mic(1,2)); %��Bx�ڿռ��ʱ���϶�����ƶ���z/2�ͦ�t/2
    B_mic(1,1)=B_mic_xin+dz/2*(1/(c^2)*(E_mic(1,2)-E_mic_yold)/dt); %��ڴ���Bx����Bx�ڿռ�������ƶ���z/2�������
    B_mic(1,2)=B_mic_yin-dz/2*(1/(c^2)*(E_mic(1,1)-E_mic_xold)/dt); %��ڴ���By����By�ڿռ�������ƶ���z/2�������
    %B_mic(1,2)=B_mic_yin-dt/(2*dz)*(E_mic(2,1)-E_mic(1,1)); %��By�ڿռ��ʱ���϶�����ƶ���z/2�ͦ�t/2
    
    % E_mic(nz,1)=0; %���ڴ�x����ĵ糡ǿ��0
    % E_mic(nz,2)=0; %���ڴ�y����ĵ糡ǿ��0
    
    B_mic_xout=B_mic(nz,1);
    B_mic_yout=B_mic(nz,2);
    
    %B_mic_xout=0; %���ڴ�x����ĴŸ�Ӧǿ��0
    % B_mic_yout=0; %���ڴ�y����ĴŸ�Ӧǿ��0
    % B_mic(nz-1,1)=B_mic_xout-dz/2*MU0*J(nz,2); %���ڴ���Bx����ǰ���
    % B_mic(nz-1,2)=B_mic_yout+dz/2*MU0*J(nz,1); %���ڴ���By����ǰ���
    
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
%���΢����ų��ֲ���ͼ��
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

