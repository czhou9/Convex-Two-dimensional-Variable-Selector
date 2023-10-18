function [Beta0,B9cv]=cv_new(a,b,m,n,Xt,Yt,Zt,lamda,wt,fd,n_fold)

cvx_begin quiet
variable B9cv(a,b);
expression wsr(a)
for j = 1:a
    wsr(j) = norm(B9cv(j,1:b),2);
end
expression wsc(b)
for k = 1:b
    wsc(k) = norm(B9cv(1:a,k),2);
end
variable Beta0
if (fd == 1)
    neglog = [];
    for i =  (fd *(n/n_fold)+1) : (n_fold *(n/n_fold))
        neglog = [neglog; log(1+exp(Beta0 + sum(dot(B9cv, Xt{i,m}))))- Yt{i,m}*(Beta0 + sum(dot(B9cv, Xt{i,m})))];
    end
else if (fd > 1) && (fd < n_fold)
        neglog = [];
        for i = [1 : ((fd-1) *(n/n_fold))    ((fd)*(n/n_fold)+1) : (n_fold *(n/n_fold))]
            neglog = [neglog; log(1+exp(Beta0 + sum(dot(B9cv, Xt{i,m}))))- Yt{i,m}*(Beta0 + sum(dot(B9cv, Xt{i,m})))];
        end
    else
        neglog = [];
        for i = 1 : ((n_fold-1) *(n/n_fold))
            neglog = [neglog; log(1+exp(Beta0 + sum(dot(B9cv, Xt{i,m}))))- Yt{i,m}*(Beta0 + sum(dot(B9cv, Xt{i,m})))];
        end
    end
end

minimize(  sum(neglog) +(wt/100)*( (lamda/100)*(sqrt(b)' * sum(wsr())+sqrt(a)' * sum(wsc()))...
    + (1 - lamda/100)*sum(sum(abs(B9cv)))))
cvx_end
%% 根据B9CV的值选择Xt{i,m},即Xt{i,m}中凡是对应B9CV里0或者接近0的元素都被去掉。然后重新拟合一个常规的logistic regression（没有被惩罚的），获得的B0和B9CV作为输出，参与外面的CV的预测和Error计算。
end