function [X_k_post,P_k_post] = Kalman_Filter(A, B, H,u_k,z_k, X_k1_post, P_k1_post,Q,R_m)

    X_k_pri = A* X_k1_post + B*u_k;
    P_k_pri = A *P_k1_post * A.' + Q;


    K_k = (P_k_pri*H.')*(H*P_k_pri*H.' + R_m)^-1;
    X_k_post = X_k_pri + K_k*(z_k - H*X_k_pri);
    P_k_post = (eye(size(K_k*H)) - K_k*H)*P_k_pri;
end

