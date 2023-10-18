%% input data
a = 10;
b = 10;
rep = 100;
sigma = 0.05;
NSR = 0;
n = 300;
k = 10;
%m_initial = 20;

%r = 5;

%% Generate B matrix 
B = 2 * rand(a,b) - 1;
B(1:2:end,:) = 0;
B(:,2:2:end) = 0;
row = [2 2 4 4 6 6 6 8 8 10 10 10];
col = [1 9 3 7 1 5 9 5 7 3  5  9];
indices = sub2ind(size(B), row, col);
B(indices) = 0;
N0 = sum(B(:)==0);
N1 = numel(B) - N0;
N = numel(B);

%% Generate Xt, Yt, Zt
for j = 1:rep
    %sigma = 0;
    for i = 1:n
    Xt{i,j} = normrnd(0,1,a,b);
    %sigma = sigma+ NSR*norm(Xt{i,j},2)/n;
    end
end

for j = 1:rep
    for i = 1:n
    Yt{i,j} = 1/(1+exp(-(sum(dot(Xt{i,j}+normrnd(0,sigma,a,b), B)))));
    if Yt{i,j} >= 0.5
        Zt{i,j} = 1;
    else
        Zt{i,j} = 0;
    end  
    end
end

% transfer Yt, Zt from cell to matrix
Yt_m = cell2mat(Yt);
Zt_m = cell2mat(Zt);

% generate balanced/unbalanced data (5:5/9:1)
index_1_90 = zeros(90,rep);
index_1_50 = zeros(50,rep);
index_0_10 = zeros(10,rep);
index_0_50 = zeros(50,rep);

for i = 1:rep
    index_Zt_1 = find(Zt_m(:,i)); 
    index_1_90(:,i) = index_Zt_1(1:90);
    index_1_50(:,i) = index_Zt_1(1:50);
    index_Zt_0 = find(Zt_m(:,i) == 0);
    index_0_10(:,i) = index_Zt_0(1:10);
    index_0_50(:,i) = index_Zt_0(1:50);
end

index_unbalanced = [index_1_90; index_0_10];
index_balanced = [index_1_50; index_0_50];

Xt_balanced = cell(100,100);
Xt_unbalanced = cell(100,100);
for i = 1:rep
    for j = 1:100
        Xt_balanced{j,i} = Xt{index_balanced(j,i),i};
        Xt_unbalanced{j,i} = Xt{index_unbalanced(j,i),i};
    end
end

Yt_m_balanced = zeros(100,100);
Yt_m_unbalanced = zeros(100,100);
Yt_m_balanced = Yt_m(index_balanced);
Yt_m_unbalanced = Yt_m(index_unbalanced);

Yt_balanced = cell(100,100);
Yt_unbalanced = cell(100,100);
for i = 1:100
    for j = 1:100
        Yt_balanced{i,j} = Yt_m_balanced(i,j);
        Yt_unbalanced{i,j} = Yt_m_unbalanced(i,j);
    end
end





% % Check the number of 0
% ZNt = zeros(1,100);
% for m = 1:100
%    ZNt(1,m) = sum(Zt_m(:,m)==1);
% end



% add noise directly to Y
% for j = 1:rep
%     for i = 1:n
%     Yt{i,j} = 1/(1+exp(-(sum(dot(Xt{i,j}, B))))) + +normrnd(0,sigma,1,1);
%     if Yt{i,j} >= 0.5
%         Zt{i,j} = 1;
%     else
%         Zt{i,j} = 0;
%     end  
%     end
% end
