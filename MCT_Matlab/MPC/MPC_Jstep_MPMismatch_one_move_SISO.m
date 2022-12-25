% to develop the 2 step prediction with one move controller
clear
format shortG
% create transfer function 1/(5s+5)^5
den_p1 = conv([5,1],[5,1]);
den_p1 = conv(den_p1, den_p1);
den_p = conv([5,1],den_p1);
% 0.8/(5s+5)^5   MODEL  esko lagi 0.8 is 100%
sys_ex_m = tf(0.8, den_p);

% Controller goes to 1.25. becoz 1.25*0.8 = 1 in controller's mind 
% But plant takes 1.25 input & has gain of 1 & goes to 1.25 in steady state
% 

% 1/(5s+5)^4(2s+1)    PLANT  % 0.8 is 1 then is 1 is 1.25 
den_p_1 = conv([2,1], den_p1);   
sys_ex_p = tf(1, den_p_1);
Ts = 5;
sys_d_m = c2d(sys_ex_m, Ts);
sys_d_p = c2d(sys_ex_p, Ts);

% No of step response coefficients to represent the process
Ncoeffs= 20;
t = 0:5:Ncoeffs*Ts;

[ym, tim] = step(sys_d_m, t);
[yp, tim] = step(sys_d_p, t);

%plot(tim,y,'*')
%title("Step response coefficients: s(1)...s(20)");

% Step response coefficients array
S = [ym(2:20);ym(20)];
Sp = [yp(2:20);yp(20)];

J = 5; % Prediction Horizon 

M = 1; % control horizon

ysp = 1; % set point at the Jth step

% Control move calculation
% du = (ysp - yfree(J))/Sj

Nsim = 20;

t_sim = 20*Ts;

t_act = 0;
du_steps = [];
time_sim = t_act;
err_k = 0;

for k = 1:Nsim  % Simulation loop - in real - life its forever
    
    % 1. Compute free response at the Jth step
    % 2. Compute control move
    % 3. Implement first move (here there is only one move)
        % Whatever we compute, we are implementing. for two moves we need
        % to do optimization
    
    % 1 . Get free response
    
    free_resp = 0;
    for i = J+1:Ncoeffs-1   % du from J+1 to N-1
        ind = k+J - i;  % below k-N doesn't influence
        if (ind) > 0   
            free_resp = free_resp + S(i)*du_steps(ind);
        else % all prev del u are zero
            break;
        end
    end
    if (k+J-Ncoeffs) > 0 % Effect of move Ncoeff before.
        free_resp = free_resp + S(Ncoeffs)*u_steps(k+J-Ncoeffs);
    end

    %% Adjust for model plant mismatch
    free_resp = free_resp + err_k(k);
    %%
    
    % Control move calculation du = (ysp - yfree(J))/ Sj
    du = (ysp - free_resp)/ S(J);
    du_steps = [du_steps du];
    du_steps(k) = du;
    if k-1 > 0
        u_k = u_steps(k-1) + du;
    else
        u_k = du;
    end
    u_steps(k) = u_k;%[u_k;u_steps(1:19)];

    % Plant measurement here is one step
    J1 = 1;
    ymeas = 0;
    for i = J1+1:Ncoeffs-1
        ind = k+J1 - i;
        if (ind) > 0
            ymeas = ymeas + Sp(i)*du_steps(ind);
        else % all prev delta u are zero
            break;
        end
    end
    if (k+J1-Ncoeffs) > 0
        ymeas = ymeas + Sp(Ncoeffs)*u_steps(k+J1-Ncoeffs);
    end
    yplant(k) = ymeas;
    
    J1 = 1;
    k1 = k + 1;
    ypred_1 = 0;
    for i = J1+1:Ncoeffs-1
        if (k1+J1 - i)>0
            ypred_1 = ypred_1 + S(i)*du_steps(k1+J1-i);
        else
            break;
        end
    end
    if (k1+J1-Ncoeffs) > 0
        ypred_1 = ypred_1 + S(Ncoeffs)*u_steps(k+1+J1-Ncoeffs);
    end

    err_k(k + 1) = ymeas - ypred_1;
    t_act = t_act + Ts;
    time_sim = [time_sim t_act];
end
figure
subplot(311)
stairs(time_sim,[0 u_steps])
title('u - control moves')
subplot(312)
plot(time_sim,[0 yplant]);  
xlabel('Time in minutes')
title('y - process output')
subplot(313)
plot(time_sim, err_k);
xlabel('Time in minutes')
title('error')