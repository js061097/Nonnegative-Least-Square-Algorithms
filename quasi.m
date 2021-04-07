function [x_quasi, tend] = quasi(A,b)
    clearvars -except A b xo;
    %A = [1 2 -1; -1 1 0;1 2 3; 1 0 1; 1 2 4]; b=[1 2 3 4 5]';
    %Initial Conditions
    
    
    tstart = cputime;
    n = size(A,2);
    x(:,1) = zeros(n,1);    %change it to zeros(n,1) later on
 
    x_ind = 1:n;
    S(:,:,1) = eye(n);
    k = 1;
    
    %Fixed and Free variable sets
    x_zero = find(x(:,1)==0);  %what if this is empty?
    F = A'*(A*x(:,1)-b);  %what if this is empty?
    z_ind = intersect(x_zero,find(F>=0));  %what about equal to condition?
    z = zeros(size(z_ind));
    y_ind = setdiff(x_ind, z_ind);      %Is there a better option to find remaining indices, cos setdiff(A,B) when A<B returns empty
    y = x(y_ind,1);
    %grad_y = F(y_ind); 
    

    %refer gam.m, alp.m, BFGS.m files for the functions
    
    An = A(:,y_ind);
    Sn = S(y_ind,y_ind,1);
    grad_y = An'*(An*y-b);
    beta = 1e-6;
    %beta = 0.000001/10;
    %beta = beta*2;
    %beta = armijo(An,b,Sn,y);
    gamma = gam(y,beta,Sn,grad_y);
    
    alpha = alp(y,gamma,An,b);
    
    if (alpha>1)
        alpha=1;
    elseif (alpha<0)
        alpha=0;
    end
    
    y = y+ alpha*(gamma - y);
    x(z_ind,2) = z;
    x(y_ind,2) = y;

    k = k+1;
    
    while (~(norm(x(:,k)-x(:,k-1))<1e-10)) %is there a better convergence condition? than ~(all(x(:,k))>0 && norm(x(:,k)-x(:,k-1))<0.01)
        u = x(:,k) - x(:,k-1);    %(all(x(:,k))>=0 && norm(x(:,k)-x(:,k-1))<0.000001)
        if(all(u==0))
            S(:,:,k) = S(:,:,k-1);
        else
            S(:,:,k) = BFGS(A,u,S(:,:,k-1));
        end
        

        x_zero = find(x(:,k)==0);  %what if this is empty?
        F = A'*(A*x(:,k)-b);  %what if this is empty?
        

        z_ind = intersect(x_zero,find(F>=0));  %what about equal to condition (F>=0)?
        z = zeros(size(z_ind));
        y_ind = setdiff(x_ind, z_ind);      %Is there a better option to find remaining indices, cos setdiff(A,B) when A<B returns empty
        y = x(y_ind,k);
        
        % grad_y = F(y_ind); MOVED TO NEW LINE BELOW An
        %{
        if(all(grad_y<0.001))
            break;
        end
        %}

        Sn = S(y_ind,y_ind,k);
        An = A(:,y_ind);
        grad_y = An'*(An*y-b);
        %beta = 0.01;
        %beta = armijo(An,b,Sn,y,grad_y);
        gamma = gam(y,beta,Sn,grad_y);
        
        alpha = alp(y,gamma,An,b);
        if (alpha>1)
            alpha=1;
        elseif (alpha<0)
            alpha=0;
        end
        
        k=k+1;

        y = y+ alpha*(gamma - y);
        x(z_ind,k) = z;
        x(y_ind,k) = y;
        disp(norm(A*x(:,k)-b));
 
    end
    disp(k);
    x_quasi = x(:,k);
    tend = cputime - tstart;
end
