% Constants
%step = 4; % Qno 8
step = 3;   % Qno 9
P = 1000;
Q = 0.001;
R = 0.01;
U = 0;
A = 1;
B = 1;
H = 1;
x=zeros(1,step+1);
K=zeros(1,step+1);
%z = [0, 0.5, 0.5085, 0.5123, 0.5234];  % Qno 8
z = [0, 0.99, 0.94, 0.96];   % Qno 9
for k=1:step
    xest = x(k) + B*U;                            % Project state ahead
    Pest = A*P(k)*transpose(A) + Q;               % Project error covariance ahead
    K(k+1) = Pest*transpose(H)*inv(H*Pest*transpose(H)+R); % Compute Kalman Gain
    x(k+1) = xest + K(k+1)*(z(k+1)-H*xest); % Update the estimate zk
    P(k+1) = (eye(size(K(k+1)*H,1))-K(k+1)*H)*Pest; % Update error covariance
end
%x(2:5)  % Qno 8
x(2:4)  % Qno 9

%function [K,x,P] = Kalman_nsteps() %x,B,A,Q,P,H,R,K,z
