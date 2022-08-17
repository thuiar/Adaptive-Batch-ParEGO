function ParetoSet = GetPS(PDec,d,lower,upper,problem)
 % This function is to obtain the rescaled Pareto Set which lies in the
 % range [-1,1]

    lower = lower(:,1:d);
    upper = upper(:,1:d);
    bounds = [lower' upper'];

    D = size(lower,2);
    PSet = zeros(D,1);
    if strcmp(problem, 'DTLZ1')||strcmp(problem,'DTLZ6')
        PSet = zeros(D,1);
    elseif strcmp(problem,'DTLZ2') ||strcmp(problem,'DTLZ3')||strcmp(problem,'DTLZ5')
        PSet = ones(D,1)*0.5;
    elseif strcmp(problem,'WFG1')||strcmp(problem,'WFG2')||strcmp(problem,'WFG3')||strcmp(problem,'WFG4')
        PSet = GetWFGPS(PDec);
        
    elseif strcmp(problem, 'UF1')
        PSet = GetUF1PS(PDec);
        
    elseif strcmp(problem, 'UF2')
        PSet = GetUF2PS(PDec);
        
    elseif strcmp(problem, 'UF3')
        PSet = GetUF3PS(PDec);
        
    elseif strcmp(problem, 'UF4')
        PSet = GetUF4PS(PDec);
        
    elseif strcmp(problem, 'UF5')
        PSet = GetUF5PS(PDec);
        
    elseif strcmp(problem, 'UF6')
        PSet = GetUF6PS(PDec);
        
    elseif strcmp(problem, 'UF7')
        PSet = GetUF7PS(PDec);
    end
    
    PSet = bsxfun(@minus, PSet, bounds(:, 1));
    PSet = bsxfun(@rdivide, PSet, bounds(:, 2) - ...
            bounds(:, 1))*2-1;
    ParetoSet = PSet;
end

function PS = GetWFGPS(PopDec)
% For WFG problems, we use tests with 3 objectives and 10 variables
    [N,D] = size(PopDec);
    M = 3;
    K = M-1;
    L = D - K;
    D = 1;
    S = 2 : 2 : 2*M;
    A = ones(1,M-1);

    z01 = PopDec./repmat(2:2:size(PopDec,2)*2,N,1);
    z01(:,K+1:D) = 0.35;

    t1 = zeros(N,K+L);
    t1(:,1:K)     = z01(:,1:K);
    t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

    t2 = zeros(N,K+L);
    t2(:,1:K)     = t1(:,1:K);
    t2(:,K+1:end) = b_flat(t1(:,K+1:end),0.8,0.75,0.85);

    t3 = zeros(N,K+L);
    t3 = b_poly(t2,0.02);

    t4 = zeros(N,M);
    for i = 1 : M-1
        t4(:,i) = r_sum(t3(:,(i-1)*K/(M-1)+1:i*K/(M-1)),2*((i-1)*K/(M-1)+1):2:2*i*K/(M-1));
    end
    t4(:,M) = r_sum(t3(:,K+1:K+L),2*(K+1):2:2*(K+L));

    x = zeros(N,M);
    for i = 1 : M-1
        x(:,i) = max(t4(:,M),A(i)).*(t4(:,i)-0.5)+0.5;
    end
    x(:,M) = t4(:,M);
    PS = x(round(rand(1,1)*N),:);
end

function PS = GetUF1PS(PopDec)
    [N,D] = size(PopDec);
    X = sin(6*pi*repmat(PopDec(:,1),1,D)+repmat(1:D,size(PopDec,1),1)*pi/D);
    PS = X(round(rand(1,1)*N),:);
end

function PS = GetUF2PS(PopDec)
    [N,D] = size(PopDec);
    J1 = 3 : 2 : D;
    J2 = 2 : 2 : D;
    X = PopDec;
    X1  = repmat(X(:,1),1,length(J1));
    X(:,J1) = (0.3*X1.^2.*cos(24*pi*X1+4*repmat(J1,size(X,1),1)*pi/D)+0.6*X1).*cos(6*pi*X1+repmat(J1,size(X,1),1)*pi/D);
    X1 = repmat(X(:,1),1,length(J2));
    X(:,J2)=(0.3*X1.^2.*cos(24*pi*X1+4*repmat(J2,size(X,1),1)*pi/D)+0.6*X1).*sin(6*pi*X1+repmat(J2,size(X,1),1)*pi/D);
    
    PS = X(round(rand(1,1)*N),:);
end

function PS = GetUF3PS(PopDec)
    [N,D] = size(PopDec);
    X = repmat(PopDec(:,1),1,D).^(0.5*(1+3*(repmat(1:D,size(PopDec,1),1)-2)/(D-2)));        
    PS = X(round(rand(1,1)*N),:);
end

function PS = GetUF4PS(PopDec)
    [N,D] = size(PopDec);
    X = sin(6*pi*repmat(PopDec(:,1),1,D)+repmat(1:D,size(PopDec,1),1)*pi/D);    
    PS = X(round(rand(1,1)*N),:);
end

function PS = GetUF5PS(PopDec)
    [N,D] = size(PopDec);
    X = sin(6*pi*repmat(PopDec(:,1),1,D)+repmat(1:D,size(PopDec,1),1)*pi/D);   
    PS = X(round(rand(1,1)*N),:);
end

function PS = GetUF6PS(PopDec)
    [N,D] = size(PopDec);
     X = sin(6*pi*repmat(PopDec(:,1),1,D)+repmat(1:D,size(PopDec,1),1)*pi/D)
    
    PS = X(round(rand(1,1)*N),:);
end

function PS = GetUF7PS(PopDec)
    [N,D] = size(PopDec);
     X = sin(6*pi*repmat(PopDec(:,1),1,D)+repmat(1:D,size(PopDec,1),1)*pi/D);  
     PS = X(round(rand(1,1)*N),:);
end

function Output = s_linear(y,A)
    Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = b_flat(y,A,B,C)
    Output = A+min(0,floor(y-B))*A.*(B-y)/B-min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
    Output = roundn(Output,-6);
end

function Output = b_poly(y,a)
    Output = y.^a;
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end
