function [Phi,F]= Mpcgain_MIMO(Am,Bm,Cm,Nc,Np)
[q,n1]=size(Cm);
[n1,m]=size(Bm);
A_e=eye(n1+q,n1+q);
A_e(1:n1,1:n1)=Am;
A_e(n1+1:n1+q,1:n1)=Cm*Am;
B_e=zeros(n1+q,m);
B_e(1:n1,:)=Bm;
B_e(n1+1:n1+q,:)=Cm*Bm;
C_e=zeros(q,n1+q);
C_e(:,n1+1:n1+q)=eye(q,q);
n=n1+q;
h(1:q,:)=C_e;
F = zeros(q*Np,n1+q);
F(1:q,:)=C_e*A_e;
for kk=2:Np
h((kk-1)*q+1:kk*q,:)=h((kk-2)*q+1:(kk-1)*q,:)*A_e;
F((kk-1)*q+1:kk*q,:)= F((kk-2)*q+1:(kk-1)*q,:)*A_e;
end
v=h*B_e;
Phi=zeros(q*Np,m*Nc); %declare the dimension of Phi
Phi(:,1:m)=v; % first column of Phi
for i=2:Nc
Phi(:,(i-1)*m+1:i*m)=[zeros((i-1)*q,m);v(1:q*(Np-i+1),1:m)]; %Toeplitz matrix
end
