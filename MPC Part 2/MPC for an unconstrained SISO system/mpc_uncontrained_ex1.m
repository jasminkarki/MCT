clear
Am=[1 1;0 1];
Bm=[0.5;1];
Cm=[1 0];
Dm=0;
sys_p = ss(Am,Bm,Cm,Dm,1);
%step(sys_p); % step response of the system

% Prediction horizon
Np=20;

% control horizon
Nc=4;

% Calculate the matrices 
[Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e]= mpcgain(Am,Bm,Cm,Nc,Np);
%[Phi_Phi,Phi_F,Phi_R,Phi, F,BarRs,A_e, B_e,C_e] = mpcgain_MIMO(Am,Bm,Cm,Nc,Np);
[n,n_in]=size(B_e);
xm=[0;0];
Xf=zeros(n,1);
N_sim=100;
r=ones(N_sim,1);
u=0; % u(k-1) =0
y=0;
R = 5; % Input weight
for kk=1:N_sim;
DeltaU=inv(Phi_Phi+R*eye(Nc,Nc))*(Phi_R*r(kk)-Phi_F*Xf);
deltau=DeltaU(1,1);
u=u+deltau;
u1(kk)=u;
y1(kk)=y;
xm_old=xm;
xm=Am*xm+Bm*u;
y=Cm*xm;
Xf=[xm-xm_old;y]; % Full state feedback in delta
end
k=0:(N_sim-1);
figure
subplot(211)
plot(k,y1)
xlabel('Sampling Instant')
legend('Output')
subplot(212)
plot(k,u1)
xlabel('Sampling Instant')
legend('Control')
title(num2str(R))