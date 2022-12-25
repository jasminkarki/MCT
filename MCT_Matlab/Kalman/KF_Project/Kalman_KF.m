%format longg
format short
%% 4 Tank Problem
step = 100;
a1 = 0.071;
a2 = 0.057;
a3 = 0.071;
a4 = 0.057;
A1 = 28;        % cm2
A2 = 32;        % cm2
A3 = 28;        % cm2
A4 = 32;        % cm2
kc = 0.50;      % V/cm
T1 = 62;        % s
T2 = 90;        % s
T3 = 23;        % s
T4 = 30;        % s
g = 981;
gamma1 = 0.7;   % dimensionless
gamma2 = 0.6 ;  % dimensionless
k1 = 3.33;      % cm3/Vs
k2 = 3.35;      % cm3/Vs
kc = 0.5;       % V/cm
A = [-1/T1 0 A3/(A1*T3) 0; 0 -1/T2 0 A4/(A2*T4); 0 0 -1/T3 0; 0 0 0 -1/T4];
B = [gamma1*k1/A1 0 ; 0 gamma2*k2/A2; 0 (1-gamma2)*k2/A3; (1-gamma1)*k1/A4 0];
z10K = readmatrix('Measurements.csv');
C = [kc 0 0 0; 0 kc 0 0];
P = 10^5*eye(4);
delt = 0.1;
%Q = 0.005 * eye(4);    % Process Noise Cov Matrix  
%R = 14* eye(2);        % Measurement Noise Cov Matrix
Q = 100* eye(4);
R = 2* eye(2);

X = [1 1 1 1]';
Uks = [3; 3];
K_gains = zeros(4,2);
predictedX = zeros(4,1);
predictedP = zeros(4,4);
trace_predictedP = 0;
updatedX = zeros(4,1);
updatedP = zeros(4,4);
trace_updatedP = 0;
innov_X = [0 0]';
residual_X = [0 0]';
for i=1:step
    %% Predictor Loop
    % Make predictions
    X = A*X + B*Uks(:,i);
    P = A*P*transpose(A) + Q;
    trace_predictedP = [trace_predictedP trace(P)];
    predictedX = [predictedX X];
    predictedP = [predictedP P];
    norm1 = norm(X);

    %% Estimator Loop
    % Calculate Kalman Gain
    K_gain = P*transpose(C)*(inv(C*P*transpose(C)+R));
    K_gains = [K_gains K_gain];

    % Update estimates
%    X = X + K_gain*(kc*z10K(i,1:2)'-C*X);
    inn_X = kc*z10K(i,1:2)'-C*X;
    innov_X = [innov_X inn_X];
    X = X + K_gain*(inn_X);
    resid_X = kc*z10K(i,1:2)'-C*X;
    residual_X = [residual_X resid_X];
    updatedX = [updatedX X];
    norm2 = norm(X);
    if abs(norm2-norm1) <= 5e-3
        sprintf('Converging %d  %f', i, abs(norm2-norm1))    
    %else
    %    sprintf('Not converging %d  %f', i, abs(norm2-norm1))    
    end
    P=(eye(size(P,2))-K_gain*C)*P;
    updatedP=[updatedP P];
    trace_updatedP = [trace_updatedP trace(P)];
    dh1bydt = (z10K(i+1,1)-z10K(i,1))'/delt;
    dh2bydt = (z10K(i+1,2)-z10K(i,2))'/delt;
    inp1 = A1/(k1*gamma1)*(dh1bydt +a1/A1*sqrt(2*g*z10K(i,1))-a3/A1*sqrt(2*g*z10K(i,3)));
    inp2 = A2/(k2*gamma2)*(dh2bydt +a2/A2*sqrt(2*g*z10K(i,2))-a4/A2*sqrt(2*g*z10K(i,4)));
    uk = [inp1; inp2];
    Uks = [Uks uk];
end

k= 1:step;
plot(1:10, predictedX(1,2:11))
legend('predictedX')
%hold on;
%plot(1:10, updatedX(1,2:11))
%title("predictedX vs updatedX")
%legend('predictedX','updatedX')

%{

plot(k, predictedX(:,2:step+1))
title("predictedX X_post")
figure;
plot(k, innov_X(:,2:step+1))
title("innov_X")
figure;
plot(k, residual_X(:,2:step+1))
title("resid_X")
figure;
%}
%plot(k, updatedX(:,2:step+1))
%title("updatedX")
%legend()
%figure;
%plot(k, trace_predictedP(k+1))
%title("predictedP")
%{
plot(k,trace_updatedP(k+1))
hold on
plot(k,trace_predictedP(k+1))
hold off
title("Predicted and Updated Covariance Matrix Plot");
figure;
%}
%{
plot(k,predictedX(:,2:step+1))
title("predictedX")pr
plot(k,trace_updatedP(k+1))
title("updatedX")
%}





