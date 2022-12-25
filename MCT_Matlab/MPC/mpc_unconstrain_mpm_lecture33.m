
Am = [1 1; 0 1];
Bm = [0.5; 1];
Cm = [1 0];
Dm = 0;

% Plant Model
Amp = [1 1; 0.2 1];
Bmp = [0.8; 1];
Cmp = [1 0];
Dmp = 0;

sys_p = ss(Am, Bm, Cm, Dm, 1);
%step(sys_p);        % Step response of the system

Np = 20;            % Prediction horizon
Nc = 4;             % Control horizon
R = .5;              % Input weight/Importance to control moves

% Calculate the matrices
[Phi_Phi, Phi_F, Phi_R,A_e,B_e,C_e] = mpcgain(Am, Bm, Cm, Nc, Np);

[n, n_in] = size(B_e);
xm = [0;0];
Xf = zeros(n,1);
N_sim = 100;
r = ones(N_sim, 1);
u = 0;
y = 0;
for kk=1:N_sim
    DeltaU = inv(Phi_Phi+R*eye(Nc,Nc))*(Phi_R*r(kk)-Phi_F*Xf);
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