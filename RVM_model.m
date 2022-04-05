global N M D
% N x M is size of design matrix and D are no.of nonzero entries in w.
N=20;
M=40;
D=7;

%2.
%Dict is design/dictionary matrix.
Dict = random('Normal',0,1,[N M]);

%w is sparse weight vector w with D randomly selected nonzero entries.
w = zeros(M,1);
w(randperm(M*1,D)) = 1;
w = w.*random('Normal',0,1,[M 1]);

x= [-20 -15 -10 -5 0]; %noise variances in dB.
var2 = 10.^-(x/20); %noise variances.
nmse = zeros(1,5); % stores final nmse corresponding to each noise variance.
p=1; % for iterating through each variance.
for(var = var2)   
    NMSE = 0;
    %doing for 500 iterations to get average value.
    for(count = 1:500)
        e = random('Normal',0,sqrt(1/var),[N 1]); %noice entries  
        
        %3.
        % vector t for each noice variance.
        t=Dict*w + e; %
        A=100*ones(M,1);
        A=diag(A); % A matrix with diagonal entries a_i(alpha_i)
        cov = (var*(Dict.')*Dict + A)^-1; %posterior covariance
        mean = var*cov*(Dict.')*t; %posterior mean
        Y=zeros(M,1); %matrix to store Y_i(Gamma_i)
        for(i=1:M)
            Y(i,1) = 1- A(i,i)*cov(i,i);
        end
        
        %4. SBL.
        % loop for getting Wmp.
        while(1)
            for(i=1:M)  %upadating a_i and Y_i.
                Y(i,1) = 1- A(i,i)*cov(i,i);
                A(i,i) = Y(i,1)/(mean(i,1)^2);
            end
            cov = (var*(Dict.')*Dict + A)^-1;
            mean_n = var*cov*(Dict.')*t;
            stop = norm(mean_n - mean)^2 / norm(mean)^2;
            mean = mean_n;
            if(stop<=0.001) %stopping condition
                break;
            end
        end
        % mean at the end of loop is Wmp.
        NMSE = NMSE + (norm(mean - w)^2 / norm(w)^2);
    end
    NMSE = NMSE/500; %average nmse
    nmse(1,p) = NMSE;
    p=p+1;
end

%5. Plot
semilogy(x,nmse);
title('NMSE vs noise variance');
xlabel('noise variance (in dB)');
ylabel('NMSE (log scale)');
