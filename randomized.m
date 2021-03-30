function [x_rand, tprep, tsmall] = randomized(A,b)

    clearvars -except A b r;
    [m,n] = size(A);
    pad = ceil(log2(m));
    m2 = 2^pad;
    A(m+1:m2,:) = zeros(m2-m,n);
    b(m+1:m2) = zeros(m2-m,1);
    [n,d] = size(A);
    r = d+20;  %Setting r value (r = d+20 is best as per paper)
    
    tstart = cputime;
    %D matrix
    dd = 2*randi(2,1,n)-3;   %assigns +1 or -1 

    %S matrix
    s = zeros(1,n);
    s(randsample(n,r)) = sqrt(n/r);
    ind = find(s~=0); %used to make mat mult easier


    %Hadamard Submatrix
    H = (1/sqrt(n))*hadamard(n); %multiply with * to normalize it

    i = 1:n;
    i = i(ind);

    for i=i
        H(i,:) = H(i,:)*s(i);
    end
    H = H(ind,:);

    HD = zeros(size(H));
    for j=1:size(H,2)
        HD(:,j) = H(:,j)*dd(j);
    end


    An = HD*A;
    bn = HD*b;
    tprep = cputime-tstart;
    
    tnew = cputime;
    x_rand = quasi(An,bn);     %Using quasi newton fnnls for the subproblem
    tsmall = cputime - tnew;

end