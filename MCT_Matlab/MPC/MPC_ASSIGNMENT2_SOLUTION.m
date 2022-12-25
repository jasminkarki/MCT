% MPC Assignment 1
clc; clear all; close all;

% SOLUTION: CHANGE NP,NC, WEIGHTS IN ANY COMBINATION TO GET THE DESIRED o/p
% ONE SUCH COMBINATION IS USED HERE

% STEP 1 : FORMULATING DISCRITIZED STATE SPACE MODEL
% 
Am=[2 0; 6 1]; 
Bm=[-1 -3]';
Cm=[0 -1];
Dm=0;
% Am =[0 1 0 ; 0.25 0 1; 0 0 0 ]
% Bm=[2 0 1]';
% Cm=[1 0 0; 1 -1 1];
% Dm=0;
sys_p = ss(Am,Bm,Cm,Dm,1);
step(sys_p); % step response of the system
% STEP 2 : COMPUTE TOPELITZ & OTHER GAIN MATRICES
Np = 7; 								% Prediction horizon	
Nc = 3; 						        % Control horizon 

% Compute matrices [Phi_Phi,Phi_F,Phi_R,A_e,B_e,C_e] 
[q,n1] = size(Cm);  					% No. of outputs = q = rows of C matrix 
[n1,m] = size(Bm);  					% No. of inputs  = m = columns of B matrix 

% Compute Augumented matrices A_e,B_e,C_e
A_e = eye(n1+q,n1+q);
A_e(1:n1,1:n1) = Am;
A_e(n1+1:n1+q,1:n1) = Cm*Am;

B_e = zeros(n1+q,m);
B_e(1:n1,:) = Bm;
B_e(n1+1:n1+q,:) = Cm*Bm;

C_e = zeros(q,n1+q);
C_e(:,n1+1:n1+q) = eye(q,q);

% Compute Topelitz matrix 
n = n1+q;
h(1:q,:) = C_e;
F = zeros(q*Np,n1+q);
F(1:q,:) = C_e*A_e;

for kk = 2:Np 
    h((kk-1)*q+1:kk*q,:) = h((kk-2)*q+1:(kk-1)*q,:)*A_e;
    F((kk-1)*q+1:kk*q,:) = F((kk-2)*q+1:(kk-1)*q,:)*A_e; 
end

v = h*B_e;
Phi = zeros(q*Np,m*Nc); 				% declare dimensions of Phi matrix 
Phi(:,1:m) = v; 						% first m columns of Phi matrix 

for i = 2:Nc
    Phi(:,(i-1)*m+1:i*m) = [zeros((i-1)*(q),m);v(1:q*(Np-i+1),1:m)]; %Toeplitz matrix 
end

Phi_Phi = Phi'*Phi;

% STEP 3 : INITIALIZE CONDITIONS FOR MPC 

q = size(Cm,1); 					% no. of outputs
[n,m] = size(B_e);	
nm = size(Am,1); 

xm = [0 0]'; 			% initial tank heights or inital state vector 

Xf = zeros(n,1); 					% Declarae Xf matrix where Xf = [x(k)-x(k-1) ; y(k)]   

N_sim = 50;						% No. of simulations 
 
u = zeros(m,1);    					% u(k-1) = 0 
y = Cm*xm;    					% inital output measurement = 0 

uwt = [0.2];						% penalty on control moves 

ysp = [1.5]; 					% desired set point for tracking output 
Rs = repmat(ysp',Np,1); 		    % Declare set point matrix 

R = zeros(m*Nc,m*Nc); 
R(1:m,1:m) = diag(uwt);				% Declare Penalty matrix 
for i = 2:Nc
    R((i-1)*m+1:i*m,(i-1)*m+1:i*m) = diag(uwt);
end

H = Phi_Phi+R;						% Compute Hessian matrix
Phi_R = Phi'*Rs;
Phi_F = Phi'*F;

uk_prev = u;						% previous control input   
r=ones(N_sim,1);

% STEP 4: KALMAN FILTER (Observer & Estimator) Initialisation

for kk=1:N_sim;
    % Free response
    yo = F*Xf;
    DeltaU=inv(Phi_Phi+R*eye(Nc,Nc))*(Phi_R*r(kk)-Phi_F*Xf);
    KMpc = (inv(Phi_Phi+R*eye(Nc,Nc))*Phi_F);
    KMPc = KMpc(1,:);
   
    % Closed loop Eigen values
    MPC_CL = eig(A_e - B_e*KMPc)
    deltau=DeltaU(1,:);
   
    % Forced response
    yf = Phi*DeltaU;
    %plot([yo+yf])
    u=u+deltau;
    u1(:,kk)=u;
    y1(:,kk)=y;

    % plant simulation
    xm_old=xm;
    xm=Am*xm+Bm*u;
    y=Cm*xm;
    Xf=[xm-xm_old;y]; % Full state feedback in delta
end
  
% Energy calculation
sum(u1.^2)
k=0:(N_sim-1);
figure
subplot(211)
plot(k,y1)
xlabel('Sampling Instant')
legend('Output')
subplot(212)
plot(k,u1)
xlabel('Sampling Instant')
legend('Control')
title(num2str(R))
