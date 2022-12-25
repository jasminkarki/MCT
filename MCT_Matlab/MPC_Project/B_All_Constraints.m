% It is applied when all three kinds of contraints are given. Those
% constraints are input, amplitude and output contraints.

% compute the constraint matrix
function [M, gamma] = B_All_Constraints(Phi, Umin, Umax, DUmin, DUmax, Ymin, Ymax, F, ukprev, xprev, Nc, Np)
m= size(ukprev, 1);  %input size

%% For Input Contraint
C1= repmat(ukprev, Nc, 1);
C2i= repmat(eye(m,m), Nc,1);
C2= zeros(m*Nc, m*Nc);
C2(:, 1:m)= C2i;

UMIN= repmat(Umin, Nc, 1);
UMAX= repmat(Umax, Nc, 1);

for i = 2: Nc
    C2((i-1)*m+1: Nc*m, (i-1)*m+1: i*m)= C2i((i-1)*m+1: Nc*m, 1:m);
end

M1= [-C2; C2];
N1= [-UMIN + C1; UMAX- C1];

%% For Amplitude Constraint

deltaMin= repmat(DUmin, Nc, 1);
deltaMax= repmat(DUmax, Nc, 1);

M2= [-eye(m*Nc, m*Nc); eye(m*Nc, m*Nc)];
N2= [-deltaMin; deltaMax];

%% For Output Constraint

yMin= repmat(Ymin, Np, 1);
yMax= repmat(Ymax, Np, 1);
M3= [-Phi; Phi];
N3= [-yMin + F*xprev; yMax- F*xprev];


%% Final return

M= [M1; M2; M3];
gamma= [N1; N2; N3];


