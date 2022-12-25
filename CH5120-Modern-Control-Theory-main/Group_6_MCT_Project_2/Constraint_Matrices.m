function [M,N] = Constraint_Matrices(Constraints,Uprev,xk,Phi,F, ni,Nc)

    M = [];
    N = [];
    C1 = zeros(ni*Nc,ni);
    C2 = zeros(ni*Nc,ni*Nc);
    
    for i=1:Nc
        C1(2*i-1:2*i,:) = eye(ni);
        for j = 1:i
            C2(2*i-1:2*i,2*j-1:2*j) = eye(ni);
        end
    end
    if isfield(Constraints,"Umin")
         Umin = repmat(Constraints.Umin,Nc,1);
         M = [M;-C2];
         N = [N; -Umin + C1*Uprev];
    end
    
    if isfield(Constraints,"Umax")
         Umax = repmat(Constraints.Umax,Nc,1);
         M = [M;C2];
         N = [N; Umax - C1*Uprev];
    end
    
    if isfield(Constraints,"DUmin")
         DUmin = repmat(Constraints.DUmin,Nc,1);
         M = [M; -eye(ni*Nc)];
         N = [N; -DUmin];
    end
    
    if isfield(Constraints,"DUmax")
         DUmax = repmat(Constraints.DUmax,Nc,1);
         M = [M; eye(ni*Nc)];
         N = [N; DUmax];
    end
    
    if isfield(Constraints,"Ymin")
        Ymin = repmat(Constraints.Ymin,Nc,1);
        M = [M; -C1*Phi];
        N = [N; -C1*Ymin + F * xk];
    end
    
    if isfield(Constraints,"Ymax")
        Ymax = repmat(Constraints.Ymax,Nc,1);
        M = [M; C1*Phi];
        N = [N; C1*Ymiax - F * xk];
    end

end