% Transfer Function to State Space 

%F(s) = (s+2)/s^3+3s^2+4s+ 6
num = [2 7 1];
den = [24 26 9 1];

[A,B,C,D] = tf2ss(num,den);


% State Space to Transfer Function
A = [0 1; -6 -5]
B = [0; 1];
C = [8 1];
D = 0;

num, den = ss2tf(A,B,C,D);

A = [0 1 0; 0 0 1; -5 -6 -4];
B = [0 ; 0; 1];
C = [1 3 2];
D = [0]

num, den = ss2tf(A,B,C,D);

A = [ -1 -1 -1; 0 1 -1; 1 -1 1];
B = [0;1;0];
C = [ 0 0 1];
D = 0;

% Wrap numerator and denominator in [ ] brackets for proper output
[num, den] = ss2tf(A,B,C,D,1);
% tf representation of the numerator and denominator
tf(num,den)

disp('Mid Sem Exam')
% State Space to Transfer Function
A = [0 1 0; 0 0 1;-3 -1 2]
B = [0; 0; 4];
C = [1 0 0];
D = [1];

[num, den] = ss2tf(A,B,C,D);


disp('Mid Sem Exam Q2')
num = [0 0 1];
den = [11.2 22 5];

[A,B,C,D] = tf2ss(num,den);
