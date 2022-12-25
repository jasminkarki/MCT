% Constants
P = [10 0; 0 1];
Q = [1 0; 0 1];
R = 1;
x = [95;1];
A = [1 1; 0 1];
B = [-0.5 1];
H = [1 0];
U = 1;
z = [0 100 97.9 94.4];
step = 3;
KalmanGain=[0;0];
xPredicted=[0;0];
pPredicted=[0 0; 0 0];
updatedStates=[0;0];
updatedCovariance=[0 0; 0 0];

for k=1:step
    xest = A*x + B*U;               % Project state ahead
    xPredicted=[xPredicted xest];   % Estimated/Predicted states
    Pest = A*P*transpose(A)+ Q;     % Project error covariance ahead
    pPredicted=[pPredicted Pest];   % Estimated/Predicted covariance
    K = Pest*transpose(H)*inv(H*Pest*transpose(H)+R); % Compute Kalman Gain
    KalmanGain = [KalmanGain K];    % Kalman
    x = xest + K*(z(k+1)-H*xest);   % Update the estimate zk
    updatedStates=[updatedStates x];% Estimated/Predicted covariance
    P = (eye(size(P,2))-K*H)*Pest;  % Update error covariance
    updatedCovariance=[updatedCovariance P];  % Estimated/Predicted covariance
end
disp("Kalman Gain:")
KalmanGain(:,2:4)




%function [K,x,P] = Kalman_nsteps() %x,B,A,Q,P,H,R,K,z

