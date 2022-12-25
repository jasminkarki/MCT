% Data from the paper:

steps = 10;
A1 = 28;        % cm2
A2 = 32;        % cm2
A3 = 28;        % cm2
A4 = 32;        % cm2
kc = 0.50;      % V/cm
g=981;

a1 = 0.071;
a2 = 0.057;
a3 = 0.071;
a4 = 0.057;
dt = 0.1; % Sampling time (s)
t = dt*(1:steps); % time vector (s)

%Minimum Phase Characteristics Values:

gamma1 = 0.7;   % dimensionless
gamma2 = 0.6 ;  % dimensionless
k1 = 3.33;      % cm3/Vs
k2 = 3.35;      % cm3/Vs

T1 = 62;        % s
T2 = 90;        % s
T3 = 23;        % s
T4 = 30;

A = [-1/T1 0 A3/(A1*T3) 0; 0 -1/T2 0 A4/(A2*T4); 0 0 -1/T3 0; 0 0 0 -1/T4];
B = [(gamma1*k1)/A1 0 ; 0 (gamma2*k2)/A2; 0 ((1-gamma2)*k2)/A3; ((1-gamma1)*k1)/A4 0];
H = [kc 0 0 0; 0 kc 0 0];

measurementData = xlsread('Measurements');

predictedX = [];
pPredict = [];
updatedX = [];
updatedCovariance = [];
kalmanGain = [];

X_post = [24;24;2;2];
%X_post=[1;1;1;1];


P_post = 100 * eye(4);
Q = 20 * eye(4);
R = 0.001 * eye(2);

Ts = 0.1;      %Sampling Time

v1_store = [];
v2_store = [];
Residual = [];
Innovation = [];


for i=1:steps
    
% Computation of input u using forward difference equation

h1=measurementData(i,1);
h2=measurementData(i,2);
h3=measurementData(i,3);
h4=measurementData(i,4);

devH1= ( measurementData(i+1,1)-measurementData(i,1) )/ Ts;
devH2= ( measurementData(i+1,2)-measurementData(i,2) )/ Ts;


v1 = (devH1 + (a1*sqrt(2*g*h1))/A1 - (a3*sqrt(2*g*h3))/A1)*(A1/(gamma1*k1));
v2 = (devH2 + (a2*sqrt(2*g*h2))/A2 - (a4*sqrt(2*g*h4))/A4)*(A2/(gamma2*k2));

v1_store = [v1_store v1];
v2_store = [v2_store v2];

U=[v1;v2];

%Prediction
    X_prior = A*X_post + B*U;
    P_prior = A*P_post*transpose(A) + Q;

    
%Calculate Kalman Gain
    KG = P_prior*transpose(H)*(inv(H*P_prior*transpose(H) + R));
    
%Measurement Data (True data)

    Z_true = H*transpose(measurementData(i,1:4));

%Estimated Measurement

    Z_est = H*X_prior;   % innovation

%Error between estimated measurement and true measurement

    Z_error = Z_true - Z_est;

    Innovation = [Innovation Z_est];

%Update with kalman gain
    
    X_post = X_prior + KG*Z_error;

    Z_res = H*X_post;   % innovation

    Residual = [Residual Z_res];

    P_post = P_prior-KG*H*P_prior;
    
    predictedX = [predictedX X_prior];
    pPredict = [pPredict P_prior];
    updatedX = [updatedX X_post];
    updatedCovariance = [updatedCovariance P_post];
    kalmanGain = [kalmanGain KG];
    
    norm1 = norm(X_prior);
    norm2 = norm(X_post);
    if abs(norm2-norm1) <= 5e-3
        sprintf('Converging %d  %f', i, abs(norm2-norm1))    
    %else
    %    sprintf('Not converging %d  %f', i, abs(norm2-norm1))    
    end
    
end

%k = 1:step;
%plot(k, predictedX)
%hold on
%plot(k,updatedX)


%% Plot the results
% Plot the states
figure(1);
subplot(211);
%plot(t, Innovation(1,1:10), 'g-', t, predictedX(1,:), 'b--', 'LineWidth', 2);
plot(t, predictedX(1,:), 'b--', 'LineWidth', 2);
hold on; plot(t, updatedX(1,:), 'r:', 'LineWidth', 1.5)
xlabel('t (s)'); ylabel('x_1'); grid on;
legend('Predicted','Updated');
%legend('Innovation','Predicted','Updated');
subplot(212);
plot(t, Innovation(2,:), 'b--', 'LineWidth', 2);
hold on; plot(t, Residual(2,:), 'r:', 'LineWidth', 1.5)
xlabel('t (s)'); ylabel('x_2'); grid on;
legend('Estimated','Residual');
