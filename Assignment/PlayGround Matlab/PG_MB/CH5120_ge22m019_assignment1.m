A =[-3 -4 7 -4 9;-8 -9 3 8 -2;-4 -6 7 -6 -5];

%% Rank of A
r = rank(A);
sprintf("The rank of the matrix is %d",r)

%% Column Space of A
col = colspace(sym(A));
disp("Column space of the matrix is:")
disp(col)

%% Null Space of A
special_sols = null(A);
disp("Null space of the matrix is:")
disp(special_sols())

%% Qno 4: Sum of eigen values of 3d matrix.
%%% The sum of eigenvalues is the trace of matrix i.e 0+1+7=8
M_4 = [0 3 2.5; 3 1 0.5; 2.5 0.5 7];
[P,D] = eig(M_4);   % P eigen vector and D diagonal eigen values
disp(sum(D))
disp(sum(eig(M_4)))

%% Qno 5: Product of eigen values of 3d matrix.
%%% The product of eigenvalues is the determinant of matrix
M_5 = [-2 -4 -6.5; -4 -4 -4.5; -6.5 -4.5 -2];
disp(det(sym(M_5)))

%% Qno 6: Left eigen vector of 3d matrix
%%% [V,D,W] = eig(A,B) returns full matrix W.
%%% Columns of W are the corresponding left eigenvectors
%%% so that W'*A = D*W'*B.
M_6 = [0 -3 3; 3 5 5; -6 8 2];
[V,D,W]=eig(M_6);
disp(W)

%% Qno 7. Singular values of matrix
%%% Singular values of A are the square roots of the nonzero eigenvalues of ATA or AAT.
M_7 = [-8 -5 4 2; -5 -5 6 8; -4 3 4 -8];
disp(sqrt(eig(M_7*transpose(M_7))))

%% Qno 8. Clockwise rotation of matrix
M_8 = [13 -7; 9 1]
rotated_M_8 = transpose(rot90(transpose(M_8)));
disp(rotated_M_8)

%% Qno 10. Eigen value decomposition of M and 
%% calculate the trace of the inverse of the eigen vectors of matrix M.
M_10 = [128 32 120; 32 187 47; 120 47 129];
[V,D,W] = eig(M_10);
disp(V)
disp(trace(inv(V)))
disp(W)
trace(inv(W))

