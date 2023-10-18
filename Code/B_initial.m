function [B_ini]=B_initial(a,b,m,n,Xt,Yt,Zt,j,k)
cvx_begin quiet
    variable B_ini;
    variable Beta0;
    neglog = [];
    for i = 1:n
        neglog = [neglog; log(1+exp(Beta0 + B_ini*Xt{i,m}(j,k)))- Yt{i,m}*(Beta0 + B_ini*Xt{i,m}(j,k))];
    end
    minimize( sum(neglog) );
    cvx_end
end