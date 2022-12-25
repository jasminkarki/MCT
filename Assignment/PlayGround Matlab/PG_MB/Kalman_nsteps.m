% Constants
Q = 0.001;
R =0.1;
U = 0;
A = 1;
B = 1;
H = 1;
x=zeros(1,5);
P=zeros(1,5);
disp(P)
P(1) = 1000;
disp(P)
K=zeros(1,5);
z = [0, 0.9, 0.8, 1.1, 1];

%function [x,P,K] = Kalman_nsteps() %x,B,A,Q,P,H,R,K,z
for k=1:4
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

