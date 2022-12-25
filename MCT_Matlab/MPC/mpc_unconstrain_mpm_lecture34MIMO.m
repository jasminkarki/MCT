clear
Am1 = [-7.923 7.923; 9.781 -12.97];   
Bm1 = [5.093 0; 0 6.288];
Cm1 = [1 0; 0 1];
Dm1 = [0 0; 0 0];
Ts = 0.05;
sys_cp = ss(Am1, Bm1, Cm1, Dm1);   % Continuous System
step(sys_cp);

% step(sys_cp)    % step response of the system
sys_d = c2d(sys_cp, Ts);        % Discrete System
step(sys_d);            
Am = sys_d.A;
Bm = sys_d.B;
Cm = sys_d.C;
Dm = sys_d.D;
% Prediction horizon
Np = 10;            % Prediction horizon
% Control Horizon
Nc = 3;             % Control horizon

%R = .5;              % Input weight/Importance to control moves

% Calculate the matrices
[Phi_Phi, Phi_F, Phi_R, Phi, F, BarRs,A_e,B_e,C_e] = mpcgainMIMO(Am, Bm, Cm, Nc, Np);
q = size(Cm,1);
[n, m] = size(B_e);
nm = size(Am,1);
xm = zeros(nm,1);
Xf = zeros(n,1);
N_sim = 50;
r = ones(N_sim, q);
u = zeros(m,1);  %u(k-1)=0
y = zeros(q,1);
%% Input Weight and Setpoint
%uwt = [0.1, 0.1];       % Experiment with this.
uwt = [10, 5];         
ysp = [2 -2];
Rs = repmat(ysp', Np,1);
R = zeros(m*Nc, m*Nc);
R(1:m, 1:m) = diag(uwt);
for i = 2:Nc
    R((i-1)*m+1:i*m, (i-1)*m+1:i*m) = diag(uwt);
end
for kk=1:N_sim
    pred_err = (Rs-F*Xf);
    DeltaU = inv(Phi_Phi+R)*Phi'*(pred_err);
    deltau = DeltaU(1:m,1);
    u = u+deltau;
    u1(:,kk) = u;       % Store to plot later
    y1(:,kk) = y;       % Store to plot later
    xm_old = xm;
    xm = Am*xm+Bm*u;
    y = Cm*xm + 0.05*rand(q,1);   % White noise
    Xf = [xm-xm_old;y];  % Full state feedback is delta
end

k = 0:Ts:(N_sim-1)*Ts;
figure
subplot(211)
plot(k,y1')
xlabel('Sampling Instant')
legend('Output')
subplot(212)
plot(k,u1')
xlabel('Sampling Instant')
legend('Control')
% title(num2str(R))

%% Interpretation
%{ 
System is highly interactive
Plot 1 gives "How well the system is interactive".
    - Both input affects both outputs in a very similar way.
    - When we want to control them, it is not going to be straightforward.
    - But, with MPC, it has prediction horizon and control horizon. 
    - [2 and -2] as added noise. How nicely the output settles down at 5.
Plot 2 gives Controllability. 
    Each input affects the other output. System is not that easy to control.

Step Response
%}