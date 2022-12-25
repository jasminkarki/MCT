function [xk1_post, pk1_post] = kalmanFilter(A, B , H, U, zk1, xk_post, pk_post, Q, R,Xo,Uo)


%Prediction
    X_prior=A*(xk_post - Xo) + B*(U-Uo) + Xo;
    P_prior=A*pk_post*transpose(A) + Q;

%Calculate Kalman Gain
    KG=P_prior*transpose(H)*(inv(H*P_prior*transpose(H) + R));


 %Update with kalman gain
    
    xk1_post= X_prior + KG*(zk1 - H*X_prior);
    pk1_post=P_prior-KG*H*P_prior;

    end