%% Data generation
a = 10;
b = 10;
lamda = 1;
m = 20;
%B = 2* rand(a,b) -1
%B = [0.9  0   -0.6   0    -0.5; 
     %0    0    0     0     0  ;
    %-1    0    0.9   0     0.5;
     %0    0   -0.8   0    -0.4;
     %0.5  0    0.9   0     0  ;
     %0    0    0     0     0  ;
     %0.3  0    0     0     0.9;
     %0    0    0     0     0  ;
    %-0.6  0   -1     0     0.6;
     %1    0    0.4   0    -1  ]
B = 2 * rand(a,b) - 1;
B(1:2:end,:) = 0;
B(:,2:2:end) = 0;
row = [2 2 4 4 6 6 6 8 8 10 10 10];
col = [1 9 3 7 1 5 9 5 7 3  5  9];
indices = sub2ind(size(B), row, col);
B(indices) = 0;

for j = 1:m
    for i = 1:100
    Xt{i,m} = 2* rand(a,b) -1;
    Yt{i,m} = sum(dot(Xt{i,m}, B)) + randn;
    end  
end
N0 = sum(B(:)==0);
N1 = numel(B) - N0;
N = numel(B);
m = 1
%% cvx with no penalty
cvx_setup
cvx_begin
    variable B1(a,b);
     neglog = [];
    for i = 1:n
        neglog = [neglog; log(1+exp(sum(dot(B1, Xt{i,1}))))- Yt{i,1}*( sum(dot(B1, Xt{i,1})))];
    end
    minimize( sum(neglog) )
cvx_end
TP1 = sum(sum(round(B1,3)&B))/N1
TN1 = (N - sum(sum(round(B1,3)|B)))/N0
Ac1 = (sum(sum(round(B1,3)&B)) + (N - sum(sum(round(B1,3)|B))))/N
table1 = [TP1 TN1 Ac1]
B1
%% cvx grouping column elements (10 groups)
cvx_begin
    variable B2(a,b);
    expression wsc(b)
    for j = 1:b
        wsc(j) = norm(B2(1:a,j),2);
    end
    neglog = [];
    for i = 1:n
        neglog = [neglog; log(1+exp(sum(dot(B1, Xt{i,1}))))- Yt{i,1}*( sum(dot(B1, Xt{i,1})))];
    end
    minimize( sum(neglog) + sqrt(10)' * sum(wsc) )
cvx_end
TP2 = sum(sum(round(B2,3)&B))/N1
TN2 = (N - sum(sum(round(B2,3)|B)))/N0
Ac2 = (sum(sum(round(B2,3)&B)) + (N - sum(sum(round(B2,3)|B))))/N
table2 = [TP2 TN2 Ac2]
B2
%% cvx grouping row elements (10 groups)
cvx_begin
    variable B3(a,b);
    expression wsr(a)
    for k = 1:a
        wsr(k) = norm(B3(a,1:b),2);
    end
    Err3 = []
    for i = 1:100
        Err3 = [Err3; Y{i} - sum(dot(B3, X{i}))]
    end
    minimize( norm(Err3,2) + lamda*(sqrt(10)' * sum(wsr)))
cvx_end
TP3 = sum(sum(round(B3,3)&B))/N1
TN3 = (N - sum(sum(round(B3,3)|B)))/N0
Ac3 = (sum(sum(round(B3,3)&B)) + (N - sum(sum(round(B3,3)|B))))/N
table3 = [TP3 TN3 Ac3]

%% cvx grouping row and column elements together
cvx_begin
    variable B4(a,b);
    expression wsr(a) 
    for j = 1:a
        wsr(j) = norm(B4(j,1:b),2);
    end
    expression wsc(b)
    for k = 1:b
        wsc(k) = norm(B4(1:a,k),2);
    end
    Err4 = []
    for i = 1:100
        Err4 = [Err4; Y{i} - sum(dot(B4, X{i}))]
    end
    minimize( norm(Err4,2) + 0.4*(sqrt(b)' * sum(wsr())+sqrt(a)' * sum(wsc())))
cvx_end
TP4 = sum(sum(round(B4,3)&B))/N1
TN4 = (N - sum(sum(round(B4,3)|B)))/N0
Ac4 = (sum(sum(round(B4,3)&B)) + (N - sum(sum(round(B4,3)|B))))/N
table4 = [TP4 TN4 Ac4]

%% cvx within group
cvx_begin
    variable B5(a,b);
    Err5 = []
    for i = 1:100
        Err5 = [Err5; Y{i} - sum(dot(B5, X{i}))]
    end
    minimize( norm(Err5,2) + lamda * sum(sum(abs(B5))))
cvx_end
TP5 = sum(sum(round(B5,3)&B))/N1
TN5 = (N - sum(sum(round(B5,3)|B)))/N0
Ac5 = (sum(sum(round(B5,3)&B)) + (N - sum(sum(round(B5,3)|B))))/N
table5 = [TP5 TN5 Ac5]

%% cvx between group and within group
cvx_begin
    variable B6(a,b);
    expression wsr(a) 
    for j = 1:a
        wsr(j) = norm(B6(a,1:b),2);
    end
    expression wsc(b)
    for k = 1:b
        wsc(k) = norm(B6(1:a,b),2);
    end
    Err6 = []
    for i = 1:100
        Err6 = [Err6; Y{i} - sum(dot(B6, X{i}))]
    end
    minimize( norm(Err6,2) + 0.1*(sqrt(b)' * sum(wsr())+sqrt(a)' * sum(wsc()))...
        + 0.9*sum(sum(abs(B6))))
cvx_end
TP6 = sum(sum(round(B6,3)&B))/N1
TN6 = (N - sum(sum(round(B6,3)|B)))/N0
Ac6 = (sum(sum(round(B6,3)&B)) + (N - sum(sum(round(B6,3)|B))))/N
table6 = [TP6 TN6 Ac6]

%% lamda selection -- from 0.1 to 1, repeat 20 times
Xt = cell(100,20);
Yt = cell(100,20);
for m = 1:20;
    for i = 1:100;
    Xt{i,m} = 2* rand(a,b) -1;
    Yt{i,m} = sum(dot(Xt{i,m}, B)) + randn;
    end  ;
end;

    TP_mean = []
    TN_mean = []
    Ac_mean = []
for n = 0.1:0.1:1
    TP = []
    TN = []
    Ac = []
    for m = 1:20  
    cvx_begin
    variable B7(10,5);
    expression wsr(10) 
    for j = 1:10
        wsr(j) = norm(B7(j,1:5),2);
    end
    expression wsc(5)
    for k = 1:5
        wsc(k) = norm(B7(1:10,k),2);
    end
    Err7 = []
    for i = 1:100
        Err7 = [Err7; Yt{i,m} - sum(dot(B7, Xt{i,m}))]
    end
    minimize( norm(Err7,2) + n*(sqrt(5)' * sum(wsr())+sqrt(10)' * sum(wsc()))...
        + (1 - n)*sum(sum(abs(B7))))
    cvx_end
        TP = [TP; sum(sum(round(B7,4)&B))/18]
        TN = [TN; (50 - sum(sum(round(B7,4)|B)))/32]
        Ac = [Ac; (sum(sum(round(B7,4)&B)) + (50 - sum(sum(round(B7,4)|B))))/50 ]
    end
    TP_mean = [TP_mean; mean(TP)]
    TN_mean = [TN_mean; mean(TN)]
    Ac_mean = [Ac_mean; mean(Ac)]
end  
lamda1 = 0.1:0.1:1
title1 = {'lamda', 'TP', 'TN', 'Accuracy'}
table_lamda1 = table(lamda1.', TP_mean, TN_mean, Ac_mean,'VariableNames',title1)

%% lamda selection -- from 0.01 to 0.2, repeat 20 times
    TP_means = []
    TN_means = []
    Ac_means = []
for ns = 0.01:0.01:0.2
    TPs = []
    TNs = []
    Acs = []
    for m = 1:20  
    cvx_begin quiet
    variable B7(10,5);
    expression wsr(10) 
    for j = 1:10
        wsr(j) = norm(B7(j,1:5),2);
    end
    expression wsc(5)
    for k = 1:5
        wsc(k) = norm(B7(1:10,k),2);
    end
    Err7 = []
    for i = 1:100
        Err7 = [Err7; Yt{i,m} - sum(dot(B7, Xt{i,m}))]
    end
    minimize( norm(Err7,2) + ns*(sqrt(5)' * sum(wsr())+sqrt(10)' * sum(wsc()))...
        + (1 - ns)*sum(sum(abs(B7))))
    cvx_end
        TPs = [TPs; sum(sum(round(B7,4)&B))/N1]
        TNs = [TNs; (N - sum(sum(round(B7,4)|B)))/N0]
        Acs = [Acs; (sum(sum(round(B7,4)&B)) + (N - sum(sum(round(B7,4)|B))))/N ]
    end
    TP_means = [TP_means; mean(TPs)]
    TN_means = [TN_means; mean(TNs)]
    Ac_means = [Ac_means; mean(Acs)]
end  
lamda2 = 0.01:0.01:0.2
title2 = {'lamda', 'TP', 'TN', 'Accuracy'}
table_lamda2 = table(lamda2.', TP_means, TN_means, Ac_means,'VariableNames',title2)
n = 0.13
%% group (row/column) selection (lamda from 0.05 to 0.15)

    TP_group = [];
    TN_group = [];
    Ac_group = [];
for n = 0.05:0.01:0.15
    TP = [];
    TN = [];
    Ac = [];
    for m = 1:20  
    cvx_begin
    variable B9(a,b);
    expression wsr(a) 
    for j = 1:a
        wsr(j) = norm(B9(j,1:b),2);
    end
    expression wsc(b)
    for k = 1:b
        wsc(k) = norm(B9(1:a,k),2);
    end
    Err9 = [];
    for i = 1:100
        Err9 = [Err9; Yt{i,m} - sum(dot(B9, Xt{i,m}))];
    end
    minimize( norm(Err9,2) + n*(sqrt(b)' * sum(wsr())+sqrt(a)' * sum(wsc()))...
        + (1 - n)*sum(sum(abs(B9))))
    cvx_end
        TP = [TP; (nnz(sum(round(B9,4)&B) ~= 0)+nnz(sum(round(B9,4)&B,2) ~= 0))/((a+b)/2)];
        TN = [TN; (nnz(sum(round(B9,4)|B) == 0)+nnz(sum(round(B9,4)|B,2) == 0))/((a+b)/2)];
        Ac = [Ac; (nnz(sum(round(B9,4)&B) ~= 0)+nnz(sum(round(B9,4)&B,2) ~= 0)...
                 + nnz(sum(round(B9,4)|B) == 0)+nnz(sum(round(B9,4)|B,2) == 0))/(a+b)];
    end
    TP_group = [TP_group; mean(TP)];
    TN_group = [TN_group; mean(TN)];
    Ac_group = [Ac_group; mean(Ac)];
end  
lamda_group = 0.05:0.01:0.15;
title_group = {'lamda', 'TP', 'TN', 'Accuracy'};
table_lamda_group = table(lamda_group.', TP_group, TN_group, Ac_group,'VariableNames',title_group)
%% within group selection (lamda from 0.01 to 0.1)
    TP_in = []
    TN_in = []
    Ac_in = []
    B_in = B([2 4 6 8 10],[1 3 5 7 9])
    N0_in = sum(B_in(:)==0)
    N1_in = numel(B_in) - N0_in
    N_in = numel(B_in)
for n = 0.01:0.01:0.1
    TP = []
    TN = []
    Ac = []
    for m = 1:20  
    cvx_begin
    variable B10(a,b);
    expression wsr(a) 
    for j = 1:a
        wsr(j) = norm(B10(j,1:b),2);
    end
    expression wsc(b)
    for k = 1:b
        wsc(k) = norm(B10(1:a,k),2);
    end
    Err10 = []
    for i = 1:100
        Err10 = [Err10; Yt{i,m} - sum(dot(B10, Xt{i,m}))]
    end
    minimize( norm(Err10,2) + n*(sqrt(b)' * sum(wsr())+sqrt(a)' * sum(wsc()))...
        + (1 - n)*sum(sum(abs(B10))))
    cvx_end
        B10_in = B10([2 4 6 8 10],[1 3 5 7 9])
        TP = [TP; sum(sum(round(B10_in,4)&B_in))/N1_in]
        TN = [TN; (N_in - sum(sum(round(B10_in,4)|B_in)))/N0_in]
        Ac = [Ac; (sum(sum(round(B10_in,4)&B_in)) + (N_in - sum(sum(round(B10_in,4)|B_in))))/N_in]
    end
    TP_in = [TP_in; mean(TP)]
    TN_in = [TN_in; mean(TN)]
    Ac_in = [Ac_in; mean(Ac)]
end  
lamda_in = 0.01:0.01:0.1
title_in = {'lamda', 'TP', 'TN', 'Accuracy'}
table_lamda_in = table(lamda_in.', TP_in, TN_in, Ac_in,'VariableNames',title_in)
