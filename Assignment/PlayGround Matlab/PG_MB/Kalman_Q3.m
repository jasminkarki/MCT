% Constants
% Suppose Q value as 1
Q = [1 0; 0 1];
R = 1;
U = 1;
A = [1 1; 0 1];
B = [-0.5 1];
H = [1 0];
KalmanGain=[0;0];

xPredicted=[0;0];
pPredicted=[0 0; 0 0];

updatedStates=[0;0];
updatedCovariance=[0 0; 0 0];

z = [100 97.9 94.4];
step = 3;

%function [K,x,P] = Kalman_nsteps() %x,B,A,Q,P,H,R,K,z
for k=1:step
    % Project state ahead
    xest = A*x + B*U;
    % Estimated/Predicted states
    xPredicted=[xPredict xest];

    % Project error covariance ahead
    Pest = A*P*transpose(A) + Q;
    % Estimated/Predicted covariance
    pPredict=[pPredict Pest];
    
    % Compute Kalman Gain
    K = Pest*transpose(H)*inv(H*Pest*transpose(H)+R);
    % Kalman
    KalmanGain = [KalmanGain K];

    % Update the estimate zk
    x = xest + K*(z(k+1)-H*xest);
    % Estimated/Predicted covariance
    updatedStates=[updatedStates x];

    % Update error covariance
    P = (eye(size(P,2))-K*H)*Pest;
    % Estimated/Predicted covariance
    updatedCovariance=[updatedCovariance P];
end
