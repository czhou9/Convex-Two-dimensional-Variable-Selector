function [V1,cvx_optval,Beta0]=OptV(U1,b,r,m,n,Xt,Yt,Zt,lamda)
U1 = rand(10,5)
lamda = 0.2
cvx_begin quiet
    variable Vf(r,b)
    expression wsc(b)
    for k = 1:b
        wsc(k) = norm(Vf(1:r,k),2);
    end
    variable Beta0
    wsr = zeros(a);
    for l = 1:a
        wsr(l) = norm(U1(l,1:r),2);
    end
    neglog_V = [];
    for i = 1:n
        neglog_V = [neglog_V; log(1+exp(Beta0 + sum(dot(U1*Vf, Xt{i,m}))))- Yt{i,m}*(Beta0 + sum(dot(U1*Vf, Xt{i,m})))];
    end
    minimize( sum(neglog_V) + lamda*(sqrt(r)' * (sum(wsc())+sum(wsr()))))
    cvx_end
V1=Vf; 
end