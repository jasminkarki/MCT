%Input Data
A = [-0.5 1; .5 0];
B = [0;1];
C = [2 1];
U = [0];
P = [1000 0; 0 1000];  % Assumed
Q = 5;
R = 1;
%Y = [.50 .52 .49 .53 .45 .54 .48 .47 .46 .53];  % Assumed
X = [0; 0];
KG = [0 0; 0 0];

xPredict=[0 ;0];
pPredict=[0 0; 0 0];
updatedStates=[0 ; 0];
updatedCovariance=[ 0 0; 0 0];

for i=1:10
    % Make predictions
    X=A*X + B*U;
    P=A*P*transpose(A) + Q;
    xPredict=[xPredict X];
    pPredict=[pPredict P];

    % Calculate Kalman Gain
    K=P*transpose(C)*(inv(C*P*transpose(C)+R));
    KG=[KG K];

    % Update estimates
    %X= X + K*(Y(i)-C*X);
    %updatedStates=[updatedStates X];
    P=(eye(size(P,2))-K*C)*P;
    updatedCovariance=[updatedCovariance P];
end
%disp('State Predictions')
%5xPredict(:,2:11)
%disp('Error Covariance Predictions')
%pPredict(:,3:22)
disp('Kalman Gain')
KG(:,3:11)
%disp('State Estimates')
%updatedStates(:,2:11)
%disp('Error Covariance Estimates')
%updatedCovariance(:,3:22)