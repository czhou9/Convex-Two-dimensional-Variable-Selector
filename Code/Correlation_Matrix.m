%% new correct method
% generate correlation matrix and cholesky decomposition
cor = zeros(10,10);
for i = 1:10
    for j = i:10
        cor(i,j) = 1/(2^(j-i));
        cor(j,i) = cor(i,j);
    end
end
L = chol(cor);
% Generae row-wise correlation Xt matrix 
n = 200;
rep = 100;
Xt_row = cell(n,rep);
for j = 1:rep
    for i = 1:n
        flag = 1;
        while flag > 0
              X = randn(10, 10);
              X = bsxfun(@minus,X, mean(X));
              [R,flag] = chol(cov(X));       
        end
X = X * inv(chol(cov(X)));
X = X * L;
%X = [X(:,1:9) randn(10,1)];
Xt_row{i,j} = X';
    end
end
% check if Xt satisfy the correlation matrix
corrcoef(Xt_row{5,1}')