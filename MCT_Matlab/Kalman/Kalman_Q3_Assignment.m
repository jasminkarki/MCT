A=[1 1; 0 1];
B=[-0.5;1];
C =[1 0];
U=1;
XkHat = [95;1];
PkHat = [10 0; 0 1];
Q = [1 0; 0 1];
R =1;
Yk = 100;
Ykplus1 = 97.9;
Ykplus2 = 94.4;

%function result = kalmanq3(XkHat,PkHat,) 
% Kalman Gain for step 1 i.e. k+1   
% Step 1
XkHat1mid = A*XkHat + B*U

% Step 2
PkHat1mid = A* PkHat*transpose(A) + Q

% Step 3 
Kk = PkHat1mid*transpose(C)*inv(C*PkHat1mid*transpose(C)+R)
disp("Kalman Gain k")
disp(Kk)

% Step 4
XkHat1 = XkHat1mid + Kk*(Yk-C*XkHat1mid)

% Step 5
PkHat1 = (eye(2)-Kk*C)*PkHat1mid


% For Kalman Gain k+2
% Step 1
XkHat2mid = A*XkHat1 + B*U

% Step 2
PkHat2mid = A* PkHat1*transpose(A) + Q

% Step 3 
Kkplus1 = PkHat2mid*transpose(C)*inv(C*PkHat2mid*transpose(C)+R)
disp("Kalman Gain k+1")
disp(Kkplus1)

% Step 4
XkHat2 = XkHat2mid + Kk*(Ykplus1-C*XkHat2mid)

% Step 5
PkHat2 = (eye(2)-Kkplus1*C)*PkHat2mid


% For Kalman Gain K+3
% Step 1
XkHat3mid = A*XkHat2 + B*U;

% Step 2
PkHat3mid = A* PkHat2*transpose(A) + Q;

% Step 3 
Kkplus2 = PkHat3mid*transpose(C)*inv(C*PkHat3mid*transpose(C)+R);
disp("Kalman Gain k+2")
disp(Kkplus2)

% Step 4
XkHat3 = XkHat3mid + Kk*(Ykplus2-C*XkHat3mid);
disp([XkHat1, XkHat2, XkHat3])

% Step 5
PkHat3 = (eye(2)-Kkplus2*C)*PkHat3mid;
%disp([PkHat1, PkHat2, PkHat3])


