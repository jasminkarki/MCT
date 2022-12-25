clear
% Include constraints
% Low Limit: U1, U2, delu1, delu2
% High Limit: U1, U2, delu1, delu2
Am1 = [-7.923 7.923; 9.781 -12.97];   
Bm1 = [5.093 0; 0 6.288];
Cm1 = [1 0; 0 1];
Dm1 = [0 0; 0 0];
Ts = 0.05;
sys_cp = ss(Am1, Bm1, Cm1, Dm1);   % Continuous System
%step(sys_cp);

% step(sys_cp)    % step response of the system
sys_d = c2d(sys_cp, Ts);        % Discrete System
%step(sys_d);            
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
N_sim = 100;
r = ones(N_sim, q);
u = zeros(m,1);  %u(k-1)=0
y = zeros(q,1);
%% Input Weight and Setpoint
%uwt = [0.1, 0.1];       % Experiment with this.
uwt = [1 1];             % R matrix
ysp = [2 1];             % Setpoint of 2 and 1
Rs = repmat(ysp', Np,1);
R = zeros(m*Nc, m*Nc);
R(1:m, 1:m) = diag(uwt);
for i = 2:Nc
    R((i-1)*m+1:i*m, (i-1)*m+1:i*m) = diag(uwt);
end

% Hessian Matrix
H = Phi_Phi+R;

%% Explanations
% Input Constraints  DUmin -0.1 to DUmax 0.2 and Umin -2 to Umax 2
% Output % U First move @ each % DU First incremental control move @ each step 
% Output 2 and 1 ma janu parne gaeracha  OK
% DU hits constraint initially 0.2 and -0.1 and tending to zero. 

% EFFECT OF TUNING : Relaxing limit to not hit the constraint i.e DU 0.5 0.5 
% Allowing system to move faster. In system performance
% DU thorai huda constaint hit garirako thyo, OUTPUT integrated +1 +1 in each step and it was taking long
%--Change-- DU lai dherai change huna allow garda, OUTPUT moved faster to the set point. 
% FREEDOM in input huda, MPC is using it to move output faster.
%%
DUmin = 0.5*[-1;-1];
DUmax = 0.5*[1;1];
Umin = 2*[-1;-1];
Umax = 2*[1;1];

ukprev = u;

for kk=1:N_sim
    %pred_err = (Rs-F*Xf);
    %DeltaU = inv(Phi_Phi+R)*Phi'*(pred_err);
    f = -(Phi')*(Rs-F*Xf);          % f function   F is the C CA matrix
    % Finding Linear Constraint    % M coming from input and output constr
    % M3 not included here in order to let students explore it.
    [M, gamma] = mpc_constraint_MIMO(Umin, Umax, DUmin, DUmax, ukprev, Nc);
    DeltaU = quadprog(H,f,M,gamma);
    deltau = DeltaU(1:m,1);
    u = u+deltau;
    deltau1(:,kk) = deltau;
    ukprev = u;
    u1(:,kk) = u;
    y1(:,kk) = y;
    xm_old = xm;
    xm = Am*xm+Bm*u;
    y = Cm*xm + 0.05*rand(q,1);   % White noise
    Xf = [xm-xm_old;y];  % Full state feedback is delta
end

k = 0:Ts:(N_sim-1)*Ts;
figure
subplot(311)
plot(k,y1')
xlabel('Sampling Instant')
legend('Output')
subplot(312)
plot(k,u1')
xlabel('Sampling Instant')
legend('Control')
subplot(313)
plot(k,deltau1')
xlabel('Sampling Instant')
legend('First incremental control')

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