tl = cputime;
%Initialization
n = size(A,2);

P = [];
N = 1:n;

x = zeros(n,1);
w = A'*(b - A*x);   %Negative Gradient
wN = w(N);
iters = 0;

%Main Loop
while(~(isempty(N) || all(wN<=1e-12)))
    %Moving active and passive set indices
    wmax = max(wN);
    t = find(w==wmax);
    P = [P t];
    N(N==t) = [];
    
    %Solving the LS problem
    z = zeros(n,1);
    clear("Ap");
    Ap = A(:,P);
    tls = cputime;
    zp = Ap\b;
    tls = cputime - tls;
    disp(tls);
    z(N) = zeros(length(N),1);
    z(P) = zp;
    
    %Inner loop handling the non-negativity constraint
    while(min(zp)<=0)
        xp = x(P);
        ratio = zeros(length(xp),1);
        for i=1:length(xp)
            ratio(i) = xp(i)/(xp(i)-zp(i));
        end
        qind = zp<=0;
        [alpha,pq] = min(ratio(qind));
        q = P(pq);
        x = x + alpha*(z-x);
        temp = P(x(P)==0);
        N = [N temp];
        P(P==temp) = [];
        
        clear("Ap");
        Ap = A(:,P);
        zp = Ap\b;
        z(N) = zeros(length(N),1);
        z(P) = zp;      
    end
    x=z;
    w = A'*(b - A*x);
    wN = w(N);
    iters = iters + 1;
end
tlhnnls = cputime - tl;


%{
while(1)
    w = A'*(b - A*x);   
    if(isempty(N) || all(wN<=0))
        break;
    end
    [wmax,t] = max(w);
    P = [P t];
    N(N==t) = [];
    %REMEMBER TO CLEAR THE RIGHT VARIABLES
    z = zeros(n,1);
    Ap = A(:,P);
    zp = lsqr(Ap,b);
    z(N) = zeros(length(N),1);
    z(P) = zp;
    if(all(zp>0))
        x = z; 
        continue
    end
    xp = x(P);
    [alpha,pq] = min(xp/(xp-zp));
    q = P(pq);
    x = x + alpha*(z-x);
    temp = (xp==0);
    N = [N P(temp)];
    P(temp) = [];
end
%}
