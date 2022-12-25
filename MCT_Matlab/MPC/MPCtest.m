u=[0]; X_prev=[0 0 0 0]'; ni=2; Nc=10;
DUmin = 5*[-1 -1]' ;
DUmax = 5*[1 1]' ;
Umin = 0*[-1 -1]' - [3 3]' ;  % u_0 = [3 3]', liniarized State Space Model
Umax = 20*[1 1]' - [3 3]';    % u_0 = [3 3]', liniarized State Space Model

    
for i=1:6
    Nsim = 20;  
    Ncoeffs= 20;
    J = 5;
    for k = 1:Nsim
        val=(k+J-Ncoeffs);
    end
    
    %[Phi{i},F{i}]= Mpcgain_MIMO(A,B,H_ctrl{i},Nc,Np);
    

%     DUmin = 5*[-1 -1]' ;
%     DUmax = 5*[1 1]' ;
%     Umin = 0*[-1 -1]' - [3 3]' ;  % u_0 = [3 3]', liniarized State Space Model
%     Umax = 20*[1 1]' - [3 3]';    % u_0 = [3 3]', liniarized State Space Model
%     

    % Constraints, u=[0,0], X_prev=[0,0,0,0], ni=2, Nc=global variable
    % s 
    [M,N] = Constraint_Matrices(Constraints,u,X_prev,ni,Nc);  %Phi{i},F{i}, ni,Nc);
        disp("Print")
    M = [];
    N = [];
    C1 = zeros(ni*Nc,ni);
    C2 = zeros(ni*Nc,ni*Nc);
    
    for i=1:Nc
        C1(2*i-1:2*i,:) = eye(ni);
        for j = 1:i
            C2(2*i-1:2*i,2*j-1:2*j) = eye(ni);
        end
    end

    if isfield(Constraints,"Umin")
         Umin = repmat(Constraints.Umin,Nc,1);
         M = [M;-C2];
         N = [N; -Umin + C1*Uprev];
    end

    if isfield(Constraints,"Umax")
         Umax = repmat(Constraints.Umax,Nc,1);
         M = [M;C2];
         N = [N; Umax - C1*Uprev];
    end

    if isfield(Constraints,"DUmin")
         DUmin = repmat(Constraints.DUmin,Nc,1);
         M = [M; -eye(ni*Nc)];
         N = [N; -DUmin];
    end
    disp(i)
    disp("M")
    disp(M)
    disp("N")
    disp(N)
end
