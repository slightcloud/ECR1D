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
L_real=0.12; %����Դz��������ֵ����λ��m
dz=lamda; %���ÿռ���
dt=dz/(2*c); %����ʱ������Ϊ��֤ģ�⾫ȷ�ȣ���֤dt<=dz/c��һ����Ϊdz/(2*c)��λ��s
nz=round(L_real/dz)+1; %�����
z_max=dz*(nz-1); %ģ������z��������ֵ
E0=1.5; %΢���糡ǿ�ȵķ�ֵ����λ��V/m����ȷ����
B0=1.0; %΢���Ÿ�Ӧǿ�ȵķ�ֵ����λ:T����ȷ����
w=2.45; %΢����Ƶ��,��λ:GHz
B_mic=zeros(nz+1,3); %΢���Ĵų�
E_mic=zeros(nz+2,3); %΢���ĵ糡

num_step=30000; %���в���
E_bound=zeros(num_step,3);

time=dt*[0:num_step-1]'; %����ʱ��
E_x_start=E0*cos(w*1e9*2*pi()*time); %Ex��ʼֵ
E_y_start=E0*sin(w*1e9*2*pi()*time); %Ey��ʼֵ
B_x_start=B0*sin(w*1e9*2*pi()*time); %Bx��ʼֵ
B_y_start=B0*cos(w*1e9*2*pi()*time); %By��ʼֵ

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

%���΢����ų��ֲ���ͼ��
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
