% Computation of MPC Gains
%{ Tutorial 1.2. The objective of this tutorial is to produce a MATLAB function for calculating ΦTΦ, ΦTF, ΦTR¯s. 
% The key here is to create F and Φ matrices. Φ matrix is a Toeplitz matrix, which is created by defining its first column, 
% and the next column is obtained through shifting the previous column.
%}
function [Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e] = mpcgain(Am, Bm, Cm, Nc, Np);
[m1, n1] = size(Cm);
[n1, n_in] = size(Bm);
A_e = eye(n1+m1, n1+m1);
A_e(1:n1, 1:n1) = Am;
A_e(n1+1:n1+m1, 1:n1) = Cm*Am;
B_e= zeros(n1+m1, n_in);
B_e(1:n1, :) = Bm;
B_e(n1+1:n1+m1, :) = Cm*Bm;
C_e = zeros(m1, n1+m1);
C_e(:, n1+1:n1+m1) = eye(m1,m1);
n = n1+m1;
h(1,:) = C_e;
F(1,:) = C_e*A_e;

for kk=2:Np
    h(kk,:) = h(kk-1,:)*A_e;
    F(kk,:) = F(kk-1,:)*A_e;
end
v = h*B_e;
Phi = zeros(Np,Nc);     % declare the dimension of Phi
Phi(:,1) = v;           % first column of Phi
for i=2:Nc
    Phi(:,i) = [zeros(i-1, 1); v(1:Np-i+1,1)];  % Toeplitz matrix
end
BarRs = ones(Np,1);
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;
Phi_R = Phi'*BarRs;
disp("Phi_F")
disp(Phi_F)
disp("Phi_R")
disp(Phi_R)
%disp("F")
%disp(F)


