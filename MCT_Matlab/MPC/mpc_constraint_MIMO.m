% Computation of constraint matrix
% From slide of Lecture 
function [M,gamma] = mpc_constraint_MIMO(Umin, Umax, DUmin, DUmax, Ukprev, Nc)
m = size(Ukprev, 1);        % inputs
% q = size(Ymin, 1);          % outputs

% Max u (Constraint)
%u(ki)   = u(ki-1) + deltau(ki)               = u(ki-1) + [1 0 0 0]deltau
%u(ki+1) = u(ki-1) + deltau(ki)+ deltau(ki+1) = u(ki-1) + [1 1 0 0]deltau

C1 = repmat(Ukprev, Nc, 1);        % Corresponds to u(ki-1)    
C2i = repmat(eye(m,m), Nc, 1);     % Corresponds to [1 0...0; 1 1 ...0;...]
C2 = zeros(m*Nc, m*Nc);
C2(:,1:m) = C2i;
UMIN = repmat(Umin, Nc, 1);
UMAX = repmat(Umax, Nc, 1);
for i=2:Nc                  % Toeplitz Matrix, also called C*B.
    C2((i-1)*m+1:Nc*m,(i-1)*m+1:i*m) = C2i((i-1)*m+1:Nc*m,1:m);
end

M1 = [-C2;  C2];               % LHS of  % -C2delU<= -Umin + C1u(ki-1)
N1 = [-UMIN+C1; UMAX-C1];      % RHS of  % +C2delU<= +Umax - C1u(ki-1)

% Max rate of change deltaU constraint
deltaMin = repmat(DUmin,Nc,1);
deltaMax = repmat(DUmax,Nc,1);

% delUmin <= delU <= delUmax        This term is rearranged to get M2 N2
M2 = [-eye(m*Nc,m*Nc);eye(m*Nc,m*Nc)];      % LHS of -I delU <= -delUmin
N2 = [-deltaMin; deltaMax];                 % RHS of +I delU <= -delUmin
disp(size(M1))
disp(size(M2))
M =[M1;M2];
gamma = [N1;N2];
% BarRs = ones(q*Np,1);
% Phi_Phi = Phi'*Phi;
% Phi_F = Phi'*F;
% Phi_R = Phi'*BarRs;
% disp("Phi_F")
% disp(Phi_F)
% disp("Phi_R")
% disp(Phi_R)
%disp("F")
%disp(F)


