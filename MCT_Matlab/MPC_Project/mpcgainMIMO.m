% Computation of MPC Gains
%{ Tutorial 1.2. The objective of this tutorial is to produce a MATLAB function for calculating ΦTΦ, ΦTF, ΦTR¯s. 
% The key here is to create F and Φ matrices. Φ matrix is a Toeplitz matrix, which is created by defining its first column, 
% and the next column is obtained through shifting the previous column.
%}
function [Phi, F] = mpcgainMIMO(Am, Bm, Cm, Nc, Np);
[q, n1] = size(Cm);         % Size of C matrix is n1; q is no of outputs coming from C matrix
[n1, m] = size(Bm);         % m is  the no of inputs; n1 is no of states
A_e = eye(n1+q, n1+q);      % Appending integrators in deltaU form: AUGMENTED MATRIX
A_e(1:n1, 1:n1) = Am;       % Add original Am in Augmented Matrix
A_e(n1+1:n1+q, 1:n1) = Cm*Am;   % Add Cm*Am in Augmented Matrix
B_e= zeros(n1+q, m);        % Appending B_e to Augmented Matrix
B_e(1:n1, :) = Bm;  
B_e(n1+1:n1+q, :) = Cm*Bm;
C_e = zeros(q, n1+q);
C_e(:, n1+1:n1+q) = eye(q,q);
n = n1+q;
h(1:q,:) = C_e;
F = zeros(q*Np, n1+q);      % [CA, CA^2, CA^3,... CA^Np]' matrix.
F(1:q,:) = C_e*A_e;

for kk=2:Np
    h((kk-1)*q+1:kk*q,:) = h((kk-2)*q+1:(kk-1)*q,:)*A_e;
    % Go through the loop upto Np times, Take previous and multiply with A
    F((kk-1)*q+1:kk*q,:) = F((kk-2)*q+1:(kk-1)*q,:)*A_e;
end
v = h*B_e;
Phi = zeros(q*Np,m*Nc);     % declare the dimension of Phi
Phi(:,1:m) = v;             % first column of Phi
for i=2:Nc                  % Toeplitz Matrix, also called C*B.
    % Compute first column by shifting the first column one at at time
    Phi(:,(i-1)*m+1:i*m) = [zeros((i-1)*(q), m); v(1:q*(Np-i+1),1:m)];  % Toeplitz matrix
end
BarRs = ones(q*Np,1);
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;
Phi_R = Phi'*BarRs;
% disp("Phi_F")
% disp(Phi_F)
% disp("Phi_R")
% disp(Phi_R)
%disp("F")
%disp(F)


