% Constants
% Suppose Q value as 1
Q = 0.1;
R =0.1;
U = 0;
A = 1;
B = 1;
H = 1;
step = 10;
x=zeros(1,step+1);
P=zeros(1,step+1);
P(1) = 1;
K=zeros(1,step+1);
z = [0, 0.5, 0.52, 0.49, 0.53, 0.45, 0.54, 0.48, 0.47, 0.46, 0.53];

%function [K,x,P] = Kalman_nsteps() %x,B,A,Q,P,H,R,K,z
for k=1:step
    % Project state ahead
    xest = x(k) + B*U;
    
    % Project error covariance ahead
    Pest = A*P(k)*transpose(A) + Q;
    
    % Compute Kalman Gain
    K(k+1) = Pest*transpose(H)*inv(H*Pest*transpose(H)+R);
    
    % Update the estimate zk
    x(k+1) = xest + K(k+1)*(z(k+1)-H*xest);
    
    % Update error covariance
    P(k+1) = (eye(size(K(k+1)*H,1))-K(k+1)*H)*Pest;
end
