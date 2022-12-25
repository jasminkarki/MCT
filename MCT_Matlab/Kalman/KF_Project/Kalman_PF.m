%% Particle Filter %%
%% 4 Tank Problem

format shortG
step = 100;
A1 = 28;  A2 = 32;  A3 = 28;  A4 = 32;    % cm2
kc = 0.50;      % V/cm
T1 = 62;  T2 = 90;  T3 = 23;  T4 = 30;        % s
a1 = 0.071;     % cm2
a2 = 0.057;     % cm2
a3 = 0.071;     % cm2
a4 = 0.057;     % cm2
g = 981;        % cm/s2
gamma1 = 0.7;   % dimensionless
gamma2 = 0.6 ;  % dimensionless
k1 = 3.33;      % cm3/Vs
k2 = 3.35;      % cm3/Vs
kc = 0.5;       % V/cm
A = [-1/T1 0 A3/(A1*T3) 0; 0 -1/T2 0 A4/(A2*T4); 0 0 -1/T3 0; 0 0 0 -1/T4];
B = [gamma1*k1/A1 0 ; 0 gamma2*k2/A2; 0 (1-gamma2)*k2/A3; (1-gamma1)*k1/A1 0];
z10K = readmatrix('Measurements.xlsx');
C = [1/kc 0 0 0; 0 1/kc 0 0];
v1 = 3; v2 = 3; uk = [v1; v2];
P = 10^5*eye(4);

N = 500;            % particle size
n = 4;              % states
L = chol(P);        % cholesky factor of P
% Adding covariance to each particle
x = (z10K(1,:)'*ones(1,N))'+ randn(N,n)*L; 
x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);
x4 = x(:,4);
%% Process Noise
Q = 100*eye(4);             % process noise
w = chol(Q)*randn(n,N);     % roughening of the prior
w1 = w(1,:);
w2 = w(2,:);
w3 = w(3,:);
w4 = w(4,:);
%% Roughening the Prior
x1 = x1 + w1';
x2 = x2 + w2';
x3 = x3 + w3';
x4 = x4 + w4';

% Initialize Prediction values
xpred = zeros(N,n);
x1pred = xpred(:,1);
x2pred = xpred(:,2);
x3pred = xpred(:,3);
x4pred = xpred(:,4);

%% Prediction Step
for i = 1:N
     x1pred(i) = -a1/A1*sqrt(2*g*x1(i)) + a3/A1*sqrt(2*g*x3(i)) + (gamma1*k1*v1)/A1 + w1(i);
     x2pred(i) = -a2/A2*sqrt(2*g*x2(i)) + a4/A2*sqrt(2*g*x4(i)) + (gamma2*k1*v2)/A2 + w2(i);
     x3pred(i) = -a3/A3*sqrt(2*g*x3(i)) + (1 - gamma2)*k2*v2/A3  +  w3(i);
     x4pred(i) = -a4/A4*sqrt(2*g*x4(i)) + (1 - gamma1)*k1*v1/A4  +  w4(i);
end
xpred  = abs(xpred);
x1pred = abs(x1pred);
x2pred = abs(x2pred);
x3pred = abs(x3pred);
x4pred = abs(x4pred);

% Importance Weights (Likelihood Function)
z1 = z10K(1,1);
z2 = z10K(1,2);
z = [z1; z2];
z_true = z * ones(1,N);
R = 100 * eye(2);
z_est  =  C*xpred';
v = z_true - z_est;
for i = 1:N
    q(i) = exp(-0.5 * (v(:,i)' * inv(R) * v(:,i)));
end

%% Normalizing the weights
 for i = 1:N
    wt(i) = q(i)/sum(q);
end

%% Resampling
M = length(wt);
Q = cumsum(wt);
indx = zeros(1, N);
T = linspace(0,1-1/N,N) + rand/N;
i = 1; j = 1;
while(i<=N && j<=M)
    while Q(j) < T(i)
        j = j + 1;
    end
    indx(i) = j;
     x1(i) = x1pred(j);
     x2(i) = x2pred(j);
     x3(i) = x3pred(j);
     x4(i) = x4pred(j);
    i = i + 1;
end

x1_pri_cum = cumsum(x1pred);
x2_pri_cum = cumsum(x2pred);
x3_pri_cum = cumsum(x3pred);
x4_pri_cum = cumsum(x4pred);
  
x1_cum = cumsum(x1);
x2_cum = cumsum(x2);
x3_cum = cumsum(x3);
x4_cum = cumsum(x4);
 
x_post = [x1 x2 x3 x4];

x1_est = mean(x1);
x2_est = mean(x2);
x3_est = mean(x3);
x4_est = mean(x4);

k= 1:N;
plot(k, x1_est(k,:))
title("predictedX")

