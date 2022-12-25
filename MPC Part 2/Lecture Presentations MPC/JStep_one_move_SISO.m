% to develop the J step prediction with one move controller
clear
den_p1 = conv([5 1],[5 1]);
den_p1 = conv(den_p1, den_p1);
den_p = conv([5 1], den_p1);
sys_ex = tf(2,den_p);
Ts = 5;
sys_d = c2d(sys_ex,Ts);

% No of step response coefficients to represent the process
Ncoeffs = 20;
t = 0:5:Ncoeffs*Ts;

[y,tim] = step(sys_d,t);

plot(tim,y,'*')

% Step response coefficients array
S = [y(2:20);y(20)];

J = 4; % Prediction Horizon 

M = 1; % control horizon

ysp = 1; % set point at the Jth step

% Control move calculation
% du = (ysp - yfree(J))/ Sj

Nsim = 20;

t_sim = 20*Ts; % no of time steps this simulation to run

t_act = 0;

%u_steps = zeros(Nsim,1);
%du_steps = zeros(Nsim-1,1);
time_sim = t_act;

for k = 1:Nsim  % Simulation loop - in real - life its forever
    
    % 1. Compute free response at the Jth step
    % 2. Compute control move
    % 3. Implement first move (here there is only one move)
    
    % 1 . Get free response
    
    free_resp = 0;
    for i = J+1:Ncoeffs-1
        ind = k+J - i;
        if (ind) > 0
            free_resp = free_resp + S(i)*du_steps(ind);
        else % all prev del u are zero
            break;
        end
    end
    if (k+J-Ncoeffs) > 0 % Effect of move Ncoeff before.
        free_resp = free_resp + S(Ncoeffs)*u_steps(k+J-Ncoeffs);
    end
    % Control move calculation
    % du = (ysp - yfree(J))/ Sj
    du = (ysp - free_resp)/ S(J);
%     dumax = 0.1;
%     if abs(du) > dumax
%         du = sign(du)*dumax;
%     end
    %du_steps = [du;du_steps(1:18)];
    du_steps(k) = du;
    if k-1 > 0
        u_k = u_steps(k-1) + du;
    else
        u_k = du;
    end
    u_steps(k) = u_k;%[u_k;u_steps(1:19)];
    
    
    % plant measurement here is one step
%     J1 = 1;
%     ymeas = 0;
%     for i = J1+1:Ncoeffs-1
%         ind = k+1+J1 - i;
%         if (ind) > 0
%             ymeas = ymeas + S(i)*du_steps(ind);
%         else % all prev delta u are zero
%             break;
%         end
%     end
%     if (k+1+J1-Ncoeffs) > 0
%         ymeas = ymeas + S(Ncoeffs)*u_steps(k+1+J1-Ncoeffs);
%     end
%   yplant(k) = ymeas;
    u = [0 u_steps u_k];
    t = [0:5:(k+1)*Ts];
    [y] = lsim(sys_d,u,t,0,'foh');
    yplant(k) = y(end);
    t_act = t_act + Ts;
    time_sim = [time_sim t_act];
end
figure
subplot(211)
stairs(time_sim,[0 u_steps])
title('u - control moves')
subplot(212)
plot(time_sim,[0 yplant]);  
xlabel('Time in minutes')
title('y - process output')