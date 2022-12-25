% to develop the 2 step prediction with one move controller
% One control horizon. ysp setpoint to go to one. 1/(5s+1)^5 so Gain=1
% In OPEN LOOP RESPONSE, 50-60 minutes to settle down. 
% Sk: unit step response coefficients. Multiply with deltaU 
% J = 5, started with higher du(inputs) to reach setpoint faster    
% J = 3, started with huge du(inputs) to reach setpoint faster
%% Can we design a controller such that y(k+J) reaches setpoint?
clear
format shortG
% create transfer function 1/(4s+1)^4
den_p1 = conv([4,1],[4,1]);
den_p = conv(den_p1, den_p1);
sys_ex = tf(1, den_p);
Ts = 5;
sys_d = c2d(sys_ex, Ts);

% No of step response coefficients to represent the process
Ncoeffs= 50;
t = 0:Ts:Ncoeffs*Ts;

[y, tim] = step(sys_d, t);

plot(tim,y,'*')  % tells OPEN LOOP SETTLING TIME
title("Step response coefficients: s(1)...s(20)");

% Step response coefficients array
S = [y(2:Ncoeffs);y(Ncoeffs)];

J = 1; % Prediction Horizon 

M = 1; % control horizon

%% set point at the Jth step, no of steps
ysp = 1;
Nsim = 15;   % 15 instances
%%

% Control move calculation
% du = (ysp - yfree(J))/ Sj

t_sim = Nsim*Ts; % no of time steps this simulation to run

t_act = 0;

%u_steps = zeros(Nsim,1);
%du_steps = zeros(Nsim-1,1);
du_steps = [];
time_sim = t_act;

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
        disp(free_resp)
        if free_resp == 0
            disp(k)
        end
    end
    % Control move calculation du = (ysp - yfree(J))/ Sj
    du = (ysp - free_resp)/ S(J);   % Control Law
    du_steps = [du_steps du];
    du_steps(k) = du;
    if k-1 > 0
        u_k = u_steps(k-1) + du;
    else
        u_k = du;           % First, Zero
    end
    u_steps(k) = u_k;       %[u_k;u_steps(1:19)];
    
    
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
    t = 0:5:(k+1)*Ts;
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