%%Low rank clean images simulation
%%p = dimension of images
%%r = rank of clean image subspace
%%n = number of noisy image samples to generate
function[X, R] = simulate_poisson(n,p,r)

%%% generate random low rank subspace
%% first generate v_1, v_2, ..., v_r whose coordinates are i.i.d uniformly 
%% distributed in [0,1]
%% Then, we find an orthonormal basis for this subspace

V = rand(p,r)+0.01*ones(p,r); %perturbation to avoid singularities

%V = orth(raw_V);
%%Normalize each v by L1 norm = sum of entries

L1norm = sum(V,1);
V = V./repmat(L1norm,p,1);

%%generate coefficients
A = rand(n, r);

%normalize coefficients
A = A./repmat(sum(A,2),1,r);
gamma = p/n;
%A= 200*((1+sqrt(gamma))^(1/2))*A; %scale up to push spike outside of bulk
A= 25*(1+sqrt(gamma))^2*A;
%% get clean images in original p-dim space: pxn matrix
X = V*A';
% clean covariance matrix;
%true_cov = cov(X'); %pxp matrix
%[U,D,V] = svd(sigma,'econ');
%true_mean = sum(X,2)/n;

%tic
%%poisson noise
R = poissrnd(X,p,n);
%R = random('Poisson',X,p,n);
%noisy_mean = sum(R,2)/n;
%noisy_cov = cov(R');
%toc
end


