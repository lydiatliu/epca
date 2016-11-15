function shr_eig = op_norm_shrink2( eig, gamma )
%%operator norm shrinkage
%from page 13 (4.4) of Donoho et al 2013 "Optimal Shrinkage..."
%available at https://arxiv.org/abs/1311.0851

l = length(eig);
shr_eig = zeros(l,1);
for i=1:l
    if eig(i)>(1+sqrt(gamma))^2
        shr_eig(i)=(eig(i)+1-gamma+sqrt((eig(i)+1-gamma).^2-4.*eig(i)))/2;
    else
        shr_eig(i)=1+sqrt(gamma);
    end
end
end


