%% Constants
A1 = 28;    A2 = 32;    A3 = 28;    A4 = 32;    
a1 = 0.071; a2 = 0.057; a3 = 0.071; a4 = 0.057;
kc = 0.50;  g=981;                                          

%% Minimum Phase Characteristics
h1o = 12.4;    h2o = 12.7;  h3o = 1.8;    h4o = 1.4;            
v1o = 3;    v2o = 3;   
Xo = [h1o; h2o; h3o; h4o];
Uo = [v1o; v2o];
T1 = 62;    T2 = 90;    T3 = 23;    T4 = 30;  
k1 = 3.33;  k2 = 3.35;  gamma1 = 0.7;  gamma2 = 0.6;    
Ac = [-1/T1     0           A3/(A1*T3)      0
       0        -1/T2       0               A4/(A2*T4) 
       0        0           -1/T3           0
       0        0           0               -1/T4];
Bc = [(gamma1*k1)/A1        0 
      0                     (gamma2*k2)/A2
      0                     (1-gamma2)*k2/A3
      (1-gamma1)*k1/A4      0];
Cc = [kc    0   0   0;
      0     kc  0   0];
Dc = 0;
msrmnts = xlsread('Measurements'); % Original Measurements
X_post = [1;1;1;1];      
P_post = 10^5 * eye(4);
Q = 20*eye(4);
R = 2*eye(2);

Ts=0.1;                         
continuousSys = ss(Ac,Bc,Cc,Dc);        % Ct time SS model          
discreteSys =   c2d(continuousSys, Ts); % Ct to Dt time
[A, B, H, D] = ssdata(discreteSys);     % SS Data          

v1_store = [];              % flow pump 1
v2_store = [];              % flow pump 2
X_prior_store = [];         % store prior X
P_prior_store = [];         % store prior covariance
X_post_store = [];          % store post X
P_post_store = [];          % store post covariance
K_store = [];               % store Kalman Gain
Innov_store = [];           % Msmnts error X prior
Residual_store = [];        % Msmnts error X post
Trace_P_prior = [];         % Store trace of P prior 
Trace_P_post = [];          % Store trace of P post
flag = false;

for i=1:1000  
% Compute u for each step using forward difference eq
dh1= (msrmnts(i+1,1)-msrmnts(i,1))/ Ts;         % h1 dot
dh2= (msrmnts(i+1,2)-msrmnts(i,2))/ Ts;         % h2 dot
h1 = msrmnts(i,1);    h2 = msrmnts(i,2);      
h3 = msrmnts(i,3);    h4 = msrmnts(i,4);      

v1 = (dh1 + (a1*sqrt(2*g*h1))/A1 - (a3*sqrt(2*g*h3))/A1)*(A1/(gamma1*k1));
v2 = (dh2 + (a2*sqrt(2*g*h2))/A2 - (a4*sqrt(2*g*h4))/A4)*(A2/(gamma2*k2));
v1_store = [v1_store v1];                      % Input U(1)
v2_store = [v2_store v2];                      % Input U(2)
U = [v1 v2]';                                  % Input U

% Prediction equations
X_prior = A*(X_post - Xo) + B*(U-Uo) + Xo;
P_prior = A*P_post*transpose(A) + Q;

% Kalman Gain calculation
K=P_prior*transpose(H)*(inv(H*P_prior*transpose(H) + R));

% True Z value
Z_true = H*transpose(msrmnts(i,1:4));

% Error between estimated measurement and true measurement
Innovation = Z_true - H*X_prior;

% Update equations
X_post = X_prior + K*Innovation;
P_post = P_prior-K*H*P_prior;

Residual = Z_true - H*X_post;

X_prior_store = [X_prior_store X_prior];
P_prior_store = [P_prior_store P_prior];
X_post_store = [X_post_store X_post];
P_post_store = [P_post_store P_post];
K_store = [K_store  K];
Innov_store = [Innov_store Innovation];
Residual_store = [Residual_store Residual];
Trace_P_post = [Trace_P_post trace(P_post)];
Trace_P_prior = [Trace_P_prior trace(P_prior)];
norm_X = norm(X_post-X_prior);
norm_P_prior = abs(norm(trace(P_prior)));
norm_P_post = abs(norm(trace(P_post)));
if flag == false & norm_X < 5e-3
    convergence = i;
    flag = true;
    sprintf("Converges %d", i)
end     
end