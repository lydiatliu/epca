function [X_c,ind,num_sel] = corr_prune(X,cor,normalization, remove_const)

%normalize and prune X for correlations

if ~exist('cor','var')
    cor = 1;
end

if ~exist('normalization','var')
    normalization = 'SD';
end

if ~exist('remove_const','var')
    remove_const = 0;
end

[n,p] = size(X);
num_miss = 0;
for i=1:p
    missing = find(X(:,i)==9);
    X(missing,i)=nanmean(X(:,i));
    num_miss = num_miss+length(missing);
end

%center data
mx = mean(X);
X = X - ones(n,1)*mx;
%standardize data
switch normalization
    case 'SD'
        D_n = var(X);
    case 'HW'
        D_n = mx.*(1-mx./2); %2*p*(1-p), where p = my/2;  
end
num_const = 0;
for i=1:p
    if D_n(i)>0
        X(:,i) = (n-1)^(-1/2)*X(:,i)*D_n(i)^(-1/2);
    else
        switch remove_const
            case 0
                X(:,i)   =  0;
                num_const = num_const+1;
            case 1 %NaNs lead to removing constants in the next few lines   
                X(:,i)   =  NaN;
        end 
    end
end
ind = zeros(p,1);
ind(1) = 1;

i = 1;
j = 2;
while (j<=p)
    if abs(X(:,i)'*X(:,j))<cor
        ind(j)=1;
        i = j;
    end
    j = j+1;
end

X_c  = X(:,ind==1);
num_sel = sum(ind);
