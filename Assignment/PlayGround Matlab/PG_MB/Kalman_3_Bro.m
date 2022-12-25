%Input Data
A=[1 1; 0 1];
B=[-0.5;1];
C=[1 0];
U=[1];

P=[10 0; 0 1];
Q=[1 0; 0 1];
R=[1];

Y=[100 97.9 94.4];
X=[0];

KG=[0;0];

xPredict=[0;0];
pPredict=[0 0; 0 0];

updatedStates=[0;0];
updatedCovariance=[0 0; 0 0];


for i=1:3
%Prediction
    X=A*X + B*U;
    P=A*P*transpose(A) + Q;
    xPredict=[xPredict X];
    pPredict=[pPredict P];

%Calculate Kalman Gain

    K=P*transpose(C)*(inv(C*P*transpose(C)+R));

    KG=[KG K];

%Update with kalman gain

    X= X + K*(Y(i)-C*X);

    updatedStates=[updatedStates X];

    P=(eye(size(P,2))-K*C)*P;

    updatedCovariance=[updatedCovariance P];
end
xPredict
pPredict
updatedStates
updatedCovariance
KG