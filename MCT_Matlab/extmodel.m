Ac = [0 1 0; 3 0 1; 0 1 0 ];
Bc= [1; 1; 3];
Cc=[0 1 0];
Dc=zeros(1,1);
Delta_t=1;
[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);

[m1,n1]=size(Cd);
[n1,n_in]=size(Bd);

A_e=eye(n1+m1,n1+m1);
disp(A_e)
A_e(n1+1:n1+m1,1:n1)=Cd*Ad;
disp(A_e)
B_e=zeros(n1+m1,n_in);
disp(B_e)
B_e(1:n1,:)=Bd;
disp(B_e)
B_e(n1+1:n1+m1,:)=Cd*Bd;
disp(B_e)
C_e=zeros(m1,n1+m1);
disp(C_e)
C_e(:,n1+1:n1+m1)=eye(m1,m1);
disp(C_e)