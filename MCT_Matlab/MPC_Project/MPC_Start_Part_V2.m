close all; clear all; clc

%% Parameter Initialization
%Parameter values of the laboratory process
A1 = 28;    A2 = 32;    A3 = 28;    A4 = 32;    % cm^2
a1=0.071;   a2=0.057;   a3=0.071;   a4=0.057;   % cm^2

kc = 1;      % V/cm
g=981;          %cm/s^2

gamma1 = 0.7;   gamma2 = 0.6 ;
k1 = 3.33;      k2 = 3.35;      % cm3/Vs

T1 = 62;    T2 = 90;    T3 = 23;    T4 = 30;    %sec

%Operating points of Linearized System
h10 = 12.4;     h20 = 12.7;     h30 = 1.8;     h40 = 1.4;
v10=3;  v20=3;

%Operating point vectors
Xo=[h10;h20;h30;h40];
Uo=[v10;v20];

%State Matrix of Linearized System in Continuous Time
Ac = [-1/T1 0 A3/(A1*T3) 0; 0 -1/T2 0 A4/(A2*T4); 0 0 -1/T3 0; 0 0 0 -1/T4];

%Input Matrix of Linearized System in Continuous Time
Bc = [(gamma1*k1)/A1 0 ; 0 (gamma2*k2)/A2; 0 ((1-gamma2)*k2)/A3; ((1-gamma1)*k1)/A4 0];


%Direct Transmission Matrix of Linearized System in Continuous Time
Dc=[0];

Ts=0.1;      %Sampling Time in second

% Switch variable for cases 
sysType=1        % =1 Measured: 3,4 and Controlled: 1,2
                % =2 Measured: 1,2 and Controlled: 3,4
                % =3 Measured: 1,4 and Controlled: 2,3
                % =4 Measured: 2,4 and Controlled: 1,3

switch sysType
    case 1
        Cc=[0 0 1 0; 0 0 0 1];
        Cm_ctrl=[1 0 0 0; 0 1 0 0];
        rki=[13.4;13.7];
        Xf = [0 0 0 0 12.4 12.7]';
        % Np=12;
        % Nc=4;
        Np = 12;
        Nc = 4;
        N_sim=1000;
        uwt = 0.5*ones(1,2);  % Input weight
%         ysp=rki-Xo(1:2,:);
        
    case 2
        Cc=[1 0 0 0; 0 1 0 0];
        Cm_ctrl=[0 0 1 0; 0 0 0 1];
        rki=[2.8;2.4];
        Xf = [0 0 0 0 1.8 1.4]';
        Np=20;
        Nc=8;
        N_sim=1000;
        uwt = 0.5*ones(1,2);  % Input weight
%         ysp=rki-Xo(3:4,:);
        
    case 3
        Cc=[1 0 0 0; 0 0 0 1];
        Cm_ctrl=[0 1 0 0; 0 0 1 0]; 
        rki=[13.7;2.8];
        Xf = [0 0 0 0 12.7 1.8]';
        Np=200;
        Nc=40;
        N_sim=2000;
        uwt = 0.5*ones(1,2);  % Input weight
%         ysp=rki-Xo(2:3,:);
        
    case 4
        Cc=[0 1 0 0; 0 0 0 1];
        Cm_ctrl=[1 0 0 0; 0 0 1 0];
        rki=[13.7;2.4];
        Xf = [0 0 0 0 12.4 1.8]';
        Np=12;
        Nc=4;
        N_sim=400;
        uwt = 0.5*ones(1,2);  % Input weight
%         ysp=rki-Xo([1,3],:);        
end

ysp=rki;

%Forming Continuous Time State Space Model in Matlab From Ac, Bc, Cc, Dc
%Matrix
cSys=ss(Ac,Bc,Cc,Dc);

%Discretizing the continuous time state space model
dSys=c2d(cSys,Ts);

%State Matrix, Input Matrix, Output Matrix and Direct Transmission Matrix
%of discrete time system
[Am Bm Cm Dm]=ssdata(dSys);

n=size(Am,1);   %No of States
q=size(Bm,2);   %No of Inputs
m=size(Cm,1);   %No of Outputs

%Define the initial condition of system:
ymki = [12.4 12.7 1.8 1.4]';
%xmki=                  %Initial condition for the original model


%Define the Simulation Parameters:

u=zeros(q,1); % u(k-1) =0
%y=zeros(m,1);

%Define Reference Parameter(Set Point)
Rs = repmat(ysp,Np,1);

%Define Design Parameters:(Input Penalty)
R = zeros(q*Nc, q*Nc);
R(1:q, 1:q) = diag(uwt);

for i = 2:Nc
    R((i-1)*q+1:i*q, (i-1)*q+1:i*q) = diag(uwt);
end

%Noise for implementing Kalman Filter
P_post_0 =  10*eye(4);
Q = 0.02*eye(n);  %Variance of process Noice
R_m = 0.001*eye(m);   %Variance of Measurement Noice

% Initial values for Kalman Filter
x_prev = [12.4 12.7 1.8 1.4]';      %zeros(n,1);  % X_post_0 = [0 0 0 0 ]' 
p_prev = P_post_0;
xk_post = x_prev;
pk_post = P_post_0;

%Define the System Constraints:
DUmin = 5*[-1;-1];
DUmax = 5*[1;1];
Umin = 0*[-1;-1];
Umax = 20*[1;1];

ukprev = u;

%Get the F and phi matrix:
[Phi,F]= mpcgainMIMO(Am,Bm,Cm_ctrl,Nc,Np);

% Hessian Matrix 
H = Phi'*Phi+R;
f = -(Phi')*(Rs-F*Xf);          % f function   F is the C CA matrix

%% Simulation Step
for kk=1:N_sim       
    f = -(Phi')*(Rs-F*Xf);          % f function   F is the C CA matrix
    % Finding Linear Constraint    % M coming from input and output constr
    % M3 not included here in order to let students explore it.
    [M, gamma] = mpc_constraint_MIMO(Umin, Umax, DUmin, DUmax, ukprev, Nc);
    DeltaU = quadprog(H,f,M,gamma);
    deltau = DeltaU(1:m,1);
    u = u+deltau;
    deltau1(:,kk) = deltau;
    ukprev = u;
    u1(:,kk) = u;
    xm= Am*xk_post + Bm*u + ones(n,q)*0.0001*rand(q,1);
    yk1=Cm*xm + 0.001*randn(m,1);
    [xk1_post, pk1_post] = kalmanFilter(Am,Bm,Cm,u,yk1,xk_post,pk_post,Q,R_m,Xo,Uo);
    
    yk1_ctrl=Cm_ctrl*xk1_post + 0.0001*rand(m,1);
    y1(:,kk) = yk1_ctrl;

    Xf=[xk1_post-xk_post;yk1_ctrl];
    
    xk_post=xk1_post;
    pk_post=pk1_post;
end


%% Visualizations
linewidth=1.2;
tickFontSize=11;
legendFontSize=11;
labelFontSize=11;
k = 0:Ts:(N_sim-1)*Ts;
figure
subplot(211)
plot(k,y1','LineWidth',linewidth)
xlabel('Sampling Instant')
ylabel('Height (cm)','FontSize',labelFontSize);
legend({'h_1','h_3'},'Location','best','FontSize',legendFontSize)
subplot(212)
plot(k,u1','LineWidth',linewidth)
xlabel('Sampling Instant')
ylabel('Input Voltage (Volt)','FontSize',labelFontSize)
legend({'u_1','u_2'},'Location','best','FontSize',legendFontSize)
titl = sprintf('Case %d Np=%d Nc=%d', sysType, Np, Nc);
sgtitle(titl) 



%% 3 subplot Visualization
linewidth=1.2;
tickFontSize=10.5;
legendFontSize=10;
labelFontSize=10;
k = 0:Ts:(N_sim-1)*Ts;
figure
subplot(311)
plot(k,y1','LineWidth',linewidth)
xlabel('Sampling Instant')
ylabel('Height (cm)','FontSize',labelFontSize);
legend({'h_1','h_2'},'Location','best','FontSize',legendFontSize)
subplot(312)
plot(k,u1','LineWidth',linewidth)
xlabel('Sampling Instant')
ylabel('Input Volt (V)','FontSize',labelFontSize)
legend({'u_1','u_2'},'Location','best','FontSize',legendFontSize)
subplot(313)
plot(k,deltau1','LineWidth',linewidth)
xlabel('Sampling Instant')
ylabel('Input Volt Change (dV)','FontSize',labelFontSize)
legend({'deltau_1','deltau_2'},'Location','best','FontSize',legendFontSize)
titl = sprintf('Case %d Np=%d Nc=%d', sysType, Np, Nc);
sgtitle(titl) 

% ax=gca;
% ax.FontSize=tickFontSize;
% set(gcf, 'position', [540, 1000, 750, 550]); %Size of the figure 


% %% Visualizations
% linewidth=1.2;
% lineCol={[0.3010 0.7450 0.9330]; [0 0.4470 0.7410] };
% tickFontSize=11;
% legendFontSize=12;
% labelFontSize=13;
% 
% 
% k = 0:Ts:(N_sim-1)*Ts;
% figure
% subplot(211)
% h=plot(k,y1','LineWidth',linewidth)
% xlabel('Time(sec)','FontSize',labelFontSize)
% ylabel('Height in cm','FontSize',labelFontSize);
% legend({'h_2','h_3'},'Location','best','FontSize',legendFontSize)
% 
% set(h, {'color'}, lineCol);
% ax=gca;
% ax.FontSize=tickFontSize;
% 
% lineCol={[1 0 0]; [0.6350 0.0780 0.1840]};
% 
% subplot(212)
% h=plot(k,u1','LineWidth',linewidth)
% xlabel('Time(sec)','FontSize',labelFontSize)
% ylabel('Input Voltage(Volt)','FontSize',labelFontSize)
% legend('u_1','u_2','Location','best','FontSize',legendFontSize)

% set(h, {'color'}, lineCol);
% ax=gca;
% ax.FontSize=tickFontSize;

% set(gcf, 'position', [540, 1000, 750, 550]); %Size of the figure 

% subplot(313)
% plot(k,deltau1')
% xlabel('Sampling Instant')
% legend('First incremental control')