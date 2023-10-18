function [B9,cvx_optval]=newOpt(a,b,m,n,Xt,Yt,lamda)
% m = 3
% lamda = 0.1;
cvx_begin quiet
    variable B9(a,b);
    expression wsr(a) 
    for j = 1:a
        wsr(j) = norm(B9(j,1:b),2);
    end
    expression wsc(b)
    for k = 1:b
        wsc(k) = norm(B9(1:a,k),2);
    end
    variable Beta0
    neglog = [];
    for i = 1:n
        neglog = [neglog; log(1+exp(Beta0 + sum(dot(B9, Xt{i,m}))))- Yt{i,m}*(Beta0 + sum(dot(B9, Xt{i,m})))];
    end
    minimize( sum(neglog) + wt * (lamda*(sqrt(b)' * sum(wsr())+sqrt(a)' * sum(wsc()))...
        + (1 - lamda)*sum(sum(abs(B9)))))
    %minimize( sum(neglog))
    
    cvx_end
end