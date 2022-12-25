%Enter Matrix
M=[5 1 1; -2 8 1; 2 -2 5]
%Get eigen matrix and diagonal matrix
[E,D] =eig(M)

syms t;
%create diagonal matrix with diagonal elements being syms variable
diag_t=[t 0 0; 0 t 0; 0 0 t]
 
%modify diagnoal_t matrix and formulate formulate exp(Dt) matrix
for i=1:3
diag_t(i,i)=exp(t*D(i,i));
end

%calculate e^At matrix as:
eAt_value=E*diag_t*inv(E)
 