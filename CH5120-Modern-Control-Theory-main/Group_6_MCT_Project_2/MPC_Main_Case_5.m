clear

%% Given Cases
Measured = [];
Controlled = [];
Y_set_point = [];
%Case_1
Measured{1} = [3, 4]; Controlled{1} = [1, 2]; Y_set_point{1} = [13.4, 13.7];
%Case_2
Measured{2} = [1, 2]; Controlled{2} = [3, 4]; Y_set_point{2} = [2.8, 2.4];
%Case_3
Measured{3} = [1, 2]; Controlled{3} = [1, 2]; Y_set_point{3} = [13.4, 13.7];
%Case_4
Measured{4} = [3, 4]; Controlled{4} = [3, 4]; Y_set_point{4} = [2.8, 2.4];
%Case_5
Measured{5} = [1, 4]; Controlled{5} = [2, 3]; Y_set_point{5} = [13.7, 2.8];
%Case_6
Measured{6} = [2, 4]; Controlled{6} = [1, 3]; Y_set_point{6} = [13.7, 2.4];


%% Constraints
Constraints = struct();
Constraints.DUmin = 5*[-1 -1]' ;
Constraints.DUmax = 5*[1 1]' ;
Constraints.Umin = 0*[-1 -1]' - [3 3]' ;  % u_0 = [3 3]', liniarized State Space Model
Constraints.Umax = 20*[1 1]' - [3 3]';    % u_0 = [3 3]', liniarized State Space Model

%%

% Initialization of Kalman_Filter
A = []; B = []; H_ctrl = []; H_measure = []; D = []; Ts = 0.1;

%Required outputs
Heights = []; Inputs = []; Sensor_Measurement =[];  Controlled_heights = [];

% MPC initialization
Np = 8;
Nc = 3  ;
N_sim=1000;



%%
for i = 5:5
    
    fprintf("Case : %d => ",i);
    % Finding State Space Model
    [A{i},B{i},H_measure{i},H_ctrl{i},D{i},h_0] = Kalman_init_(Measured{i}, Controlled{i});
    
    no = size(H_measure{i},1); % n0 - Number of Outputs/Measured States
    [ns,ni] = size(B{i}); % ns - Number of States; ni = Number of Inputs
   
    %Weights and Noises
    uwt = 0.5*ones(1,ni);  % Input weight
    Q = 0.02*eye(ns);  %Variance of process Noice
    R_m = 0.001*eye(no);   %Variance of Measurement Noice
    
    
    % Finding required Matrices for Cost Function
    [Phi{i},F{i}]= Mpcgain_MIMO(A{i},B{i},H_ctrl{i},Nc,Np);
    
   
    P_post_0 =  10*eye(4);
    X_prev = zeros(ns,1);  % X_post_0 = [0 0 0 0 ]' 
    Xf{i}=zeros(ns+no,1); %initial Xf = [0 0 0 0 0 0 ]' ns = 4, no = 2
    P_prev = P_post_0;
    u=zeros(ni,1); % u(k-1) =0
    y=zeros(no,1);
    
    ysp{i} = Y_set_point{i}-h_0(Controlled{i});
    Rs{i} = repmat(ysp{i}',Np,1);
    
    
    R = []; %penality matrix for deltaU
    R{i} = zeros(ni*Nc, ni*Nc);

    for j = 1:Nc
        R{i}((j-1)*ni+1:j*ni, (j-1)*ni+1:j*ni) = diag(uwt);
    end
    
    for tc=1:N_sim %tc = current time
        
        H = 0.5*(Phi{i}'*Phi{i} + R{i});
        f = (-Phi{i}')*(Rs{i} - F{i}*Xf{i});
        
        % Finding Constraint Matrices
        [M,N] = Constraint_Matrices(Constraints,u,X_prev,Phi{i},F{i}, ni,Nc);
        
        %Finding change in Input by optimizing the Cost Function using quadprog function
        options =  optimset('Display', 'off');

        DeltaU = quadprog(H,f,M,N,[],[],[],[],[], options);
        deltau=DeltaU(1:ni,1);
        u=u+deltau;
        X_m=A{i}*X_prev+B{i}*u + ones(ns,ni)*0.0002*rand(ni,1) ;
        y = H_measure{i} * X_m + 0.001*randn(no,1);  %Measured tank heights
        
        [X_post,P_post] = Kalman_Filter(A{i}, B{i}, H_measure{i},u,y, X_prev , P_prev,Q,R_m);

        y_c = H_ctrl{i}*X_post +0.001*randn(no,1); 
        
        y1{i}(:,tc)=y_c;
        Xf{i}=[X_post-X_prev;y_c]; % Full state feedback in delta
        
        X_prev = X_post;
        P_prev = P_post;
        
        Heights{i}(:,tc) = X_prev + h_0';
        Inputs{i}(:,tc) = u +[3,3]' ;
        Sensor_Measurement{i}(:,tc) = y + h_0(Measured{i})';
        Controlled_heights{i}(:,tc) = y_c + h_0(Controlled{i})';
        
        
    end
    k= Ts:Ts:(N_sim)*Ts;
    figure()
    
    subplot(2,2,1)
    plot(k, Sensor_Measurement{i})
    xlabel('Time')
    ylabel('Height')
    txt = ['h_',num2str(Measured{i}(1));'h_',num2str(Measured{i}(2)) ];
    legend(txt)
    title('Measured Tank Heights')
    grid()
 
    subplot(2,2,2)
    plot(k,Controlled_heights{i})
    xlabel('Time')
    ylabel('Height')
    txt = ['h_',num2str(Controlled{i}(1));'h_',num2str(Controlled{i}(2)) ];
    legend(txt)
    title('Controlled Tank Heights')
    grid()
    
    subplot(2,2,[3,4])
    plot(k,Inputs{i})
    xlabel('Time')
    ylabel('Voltage')
    legend('u_1','u_2')
    title('Input Voltage')
    grid()
    
    txt = ['Case : ',num2str(i)];
    sgtitle(txt, Fontweight = 'bold')
    %title(num2str(R))
end
    

