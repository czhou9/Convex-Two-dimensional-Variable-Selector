function [Beta0_2,B9cv2]=cv_new2(a,b,m,n,Xt,Yt,Zt,lamda,wt,fd,n_fold,B9cv)

cvx_begin quiet

variable B9cv2(a,b);
variable Beta0_2;


if (fd == 1)
    neglog = [];
    for i =  (fd *(n/n_fold)+1) : (n_fold *(n/n_fold))
        neglog = [neglog; log(1+exp(Beta0_2 + sum(dot(B9cv2, logical(round(B9cv,4)).*Xt{i,m}))))- Yt{i,m}*(Beta0_2 + sum(dot(B9cv2, logical(round(B9cv,4)).*Xt{i,m})))];
    end
else if (fd > 1) && (fd < n_fold)
        neglog = [];
        for i = [1 : ((fd-1) *(n/n_fold))    ((fd)*(n/n_fold)+1) : (n_fold *(n/n_fold))]
            neglog = [neglog; log(1+exp(Beta0_2 + sum(dot(B9cv2, logical(round(B9cv,4)).*Xt{i,m}))))- Yt{i,m}*(Beta0_2 + sum(dot(B9cv2, logical(round(B9cv,4)).*Xt{i,m})))];
        end
    else
        neglog = [];
        for i = 1 : ((n_fold-1) *(n/n_fold))
            neglog = [neglog; log(1+exp(Beta0_2 + sum(dot(B9cv2, logical(round(B9cv,4)).*Xt{i,m}))))- Yt{i,m}*(Beta0_2 + sum(dot(B9cv2, logical(round(B9cv,4)).*Xt{i,m})))];
        end
    end
end

minimize(  sum(neglog) )
end