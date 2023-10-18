function [U1,Beta0]=OptU(V1,a,r,m,n,Xt,Yt,Zt,lamda)

cvx_begin quiet
    variable Uf(a,r)
    expression wsr(a)
    for k = 1:a
        wsr(k) = norm(Uf(k,1:r),2);
    end
    variable Beta0
    wsc = zeros(b);
    for l = 1:b
        wsc(l) = norm(V1(1:r,l),2);
    end
    neglog_U = [];
    for i = 1:n
        neglog_U = [neglog_U; log(1+exp(Beta0+sum(dot(Uf*V1, Xt{i,m}))))- Yt{i,m}*(Beta0+sum(dot(Uf*V1, Xt{i,m})))];
    end
    minimize( sum(neglog_U) + lamda*(sqrt(r)' * (sum(wsr())+sum(wsc()))))
    cvx_end
U1=Uf;

end