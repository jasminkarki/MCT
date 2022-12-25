% s = tf('s')
syms s
A = [-3 -6; 2 5];
B = [5;8];
C = [-4 5];
M = [(s -5 )/(s^2 - 2*s - 3) 2/(s^2 - 2*s - 3) ; -6/(s^2 - 2*s - 3)  (s+3)/(s^2 - 2*s - 3)];

for i =1:2
    for j =1:2
        discreteM(i,j) = ilaplace(M(i,j));
    end
end

