function v_new=UpdateVelocity(E,B,v_old,dt)
global m_e q_e

%np=size(B,1);

t=q_e/m_e*B*0.5*dt; %t vector
t_mag2=sum(t'.^2); %t*t
s=2*t./(1+t_mag2'); %s vector
v_minus=v_old+q_e/m_e*E*0.5*dt; %v_minus
v_prime=v_minus+cross(v_minus,t); %v_prime
v_plus=v_minus+cross(v_prime,s); %v_plus
v_new=v_plus+q_e/m_e*E*0.5*dt; %v_n+1/2