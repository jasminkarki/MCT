% Assignment 5

%{
A = [0.9851 -0.9934 0
     0.0099  0.9950 0
     0.9934  99.5021 1];
B = [0.0099
     0.0000 
     0.0050];
C = [0 0 1];
%}
%Am = [0.9851   -0.9934; 0.009934 0.995];
%Bm = [0.009934; 4.979e-05];
%Cm = [0 100];
%Dm = [0];

%Np = 200;            % Prediction horizon
%Nc = 3;             % Control horizon
%R = 0.5;              % Input weight/Importance to control moves
%Ts = 0.01;

Am = [2 0; 6 1];
Bm = [-1; -3];
Cm = [0 -1];
Dm = 0;
set_point = 1.5;
%sys_p = ss(Am, Bm, Cm, Dm, Ts);
%step(sys_p);        % Step response of the system

Np = 3;            % Prediction horizon
Nc = 1;             % Control horizon
R = 5;              % Input weight/Importance to control moves

% Calculate the matrices
[Phi_Phi, Phi_F, Phi_R,A_e,B_e,C_e] = mpcgain(Am, Bm, Cm, Nc, Np);

[n, n_in] = size(B_e);
xm = [0;0];
%Xf = [0.1 0.2 0.3]';
%zeros(n,1);
N_sim = 100;
r = set_point* ones(N_sim, 1);
u = 0;
y = 0;
Xf = 0;
for kk=1:N_sim
    % delta_U = inv(Phi_Phi+R*eye(Nc,Nc))*(Phi_R*r(kk)-Phi_F*Xf);
    inverse_val = inv(Phi_Phi+R*eye(Nc,Nc));
    Ky = [1 0 0]*inverse_val*(Phi_R);
    Kmpc = inverse_val*(Phi_F);
    DeltaU = Ky*r(kk)-Kmpc*Xf;
    % DeltaU = inv(Phi_Phi+R*eye(Nc,Nc))*(Phi_R*r(kk))-inv(Phi_Phi+R*eye(Nc,Nc))*(Phi_F*Xf);
    disp(kk)
    deltau = DeltaU(1,1);
    disp(deltau)
    u = u+deltau;
    u1(kk) = u;
    y1(kk) = y;
    xm_old = xm;
    xm = Am*xm+Bm*u;
    y = Cm*xm;
    Xf = [xm-xm_old;y];  % Full state feedback is delta
end

k = 0:(N_sim-1);
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


% Plant Model
Amp = [1 1; 0.2 1];
Bmp = [0.8; 1];
Cmp = [1 0];
Dmp = 0;
