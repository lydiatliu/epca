%process HGDP dataset: transform from PLINK to Matlab format
%dataset provided by Nick Patterson
cd('C:\Git\poissonpca\experiments\snp_data\process_HGDP')
beep off
%% load data
fid=fopen('yfhx.geno','r');
r = textscan(fid,'%s');
r = r{:};
fclose(fid);
p  = length(r);
n  = 60;
snp = zeros(n,p);
t = tic;
for i=1:p
    x = r{i};
    for j=1:n
        snp(j,i) = str2double(x(j));
    end
end
toc(t);
%% 
%original
ind1 = [(1:10)';(21:30)'; (41:50)'];
X = snp(ind1,:);
%new sample
ind2 = [(11:20)';(31:40)'; (51:60)'];
Y = snp(ind2,:);
%%
pops = [ones(10,1); 2*ones(10,1); 3*ones(10,1)];
pop_code = cell(3,1);
pop_code{1} = 'French'; 
pop_code{2} = 'Han'; 
pop_code{3} = 'Yoruba';
%%
save('HGDP_small.mat','X','Y','pops','pop_code');

