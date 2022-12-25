A = [1 2; 1 1];
B = [1;0];
C = [1 1];
D = [0];
poles =[-5,-4];
k =place(A,B,poles)

% Q2
P = [1 2 3; 1 0 0 ; 0 1 0];
Q = [1;0;0];
pols = [0.1, 0.4+j*0.4, 0.4-j*0.4];
k = place(P,Q,pols)

% Q3
A = [0 1 0; 0 0 1; -4 -2 -1];
B = [0;0;1];
poles = [0, -0.5-j*0.5, -0.5+j*0.5];
k = place(A,B,poles)

%Q4 Chegg
A = [0 1 0; 0 0 1; -6 -11 -6];
B = [0;0;1];
C = [1 0 0];
poles = [-1+2*j, -1-2*j, -5];
k = place(A,B,poles)


disp('Mid Semester Exam')
A = [-6 0; -1 -7];
B = [1;0];
poles = [-3,-5];
k = place(A,B,poles)