function [A,B,H_measure,H_ctrl,D,h_0] = Kalman_init_(Measurement_index, Control_index)
    Area = [28 32 28 32]; %cm^2; A_1 = Area(1), A_2 = Area(2),...
    area = [0.071 0.057 0.071 0.057]; %cm^2; a_1 = area(1), a_2 = area(2),...
    h_0 = [12.4 12.7 1.8 1.4]; %cm^2; h1_0 = h_0(1), h2_0= h_0(2),...
    k_c = 1; %V/cm
    g = 981; %cm/s^2
    T = zeros(1,4);
    for i=1:4
        T(i) = (Area(i)/area(i))*sqrt((2*h_0(i))/g);
    end

    gamma_1 = 0.7; gamma_2 = 0.6;
    k_1 = 3.33; k_2 = 3.35;

    A_c = [-1/T(1) 0 Area(3)/(Area(1)*T(3)) 0; 0 -1/T(2) 0 Area(4)/(Area(2)*T(4)); 0 0 -1/T(3) 0; 0 0 0 -1/T(4)];
    B_c = [(gamma_1*k_1)/Area(1) 0; 0 (gamma_2*k_2)/Area(2); 0 ((1-gamma_2)*k_2)/Area(3); ((1-gamma_1)*k_1)/Area(4) 0];
    H_m = zeros(length(Measurement_index),4);
    H_ctrl = zeros(length(Control_index),4);
    for i = 1:4
        if i<= length(Measurement_index)
            H_m(i,Measurement_index(i)) = k_c;
        end
        if i<= length(Control_index)
            H_ctrl(i,Control_index(i)) = k_c;
        end
        
    end
        
    %H_m = [k_c 0 0 0;0 k_c 0 0] or ..
    %H_c = [0 0 k_c 0;0 0 0 k_c] or ..
    
    D_c = 0;
    ssc = ss(A_c,B_c,H_m,D_c); %Statespace Model in Continuous domine
    sys = c2d(ssc,0.1);

    [A,B,H_measure,D] = ssdata(sys); %Statespace Model in Discrete domine
    
    Co_M = ctrb(A,B); %Controlability Matrix
    Ob_M = obsv(A,H_measure); %Observability Matrix
    unco = length(A) - rank(Co_M);
    unob = length(A) - rank(Ob_M);
    
    if unco == 0
        fprintf("The System is Controllable and ");
    else
        fprintf("The System is Uncontrollable and ");
    end


    if unob == 0
        fprintf("Observable \n");
    else
        fprintf("Unobservable \n");
    end
end