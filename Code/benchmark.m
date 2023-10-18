%% Genrate initialization of U and V
% transfer cell to matrix
Xtt = rand(a,b,n,rep);
Ytt = rand(a,b,n,rep);
Ztt = rand(a,b,n,rep);
for m = 1:rep
    for i = 1:n
        for j = 1:a
            for k = 1:b
                
                Xtt(j,k,i)= Xt{i,m}(j,k);
                Ytt(j,k,i)= Yt{i,m}(j,k);
                Ztt(j,k,i)= Zt{i,m}(j,k);
            end
        end
    end
end
                
 
B_initial_rep = cell(rep);
for m = 1:rep
    B_initial_run = zeros(a,b);
    parfor j = 1:a
        for k = 1:b
    [B_ini]=B_initial(a,b,m,n,Xt,Yt,Zt,j,k)
     B_initial_run(j,k) = B_ini;;
        end
    end
     B_initial_rep{m} = B_initial_run;
end
%B_initial_avg = mean();

    
    

%% Select lamda and rank based on AIC
AIC_UV = [];
for lamda = 0.1:0.1:1
    AIC_initial_best = [];
for m = 1:rep
    sc = 10;
    opt_0 = 100;
    iterations = 0;
    Beta0 = 0;
    AIC_initial = [];
    
    for i = 1: m_initial
        V1 = V{i};
while sc > 10^(-5)
    [U1,Beta0]=OptU(V1,a,r,m,n,Xt,Yt,Zt,lamda);
    
    [V1,cvx_optval]=OptV(U1,b,r,m,n,Xt,Yt,Zt,lamda);
    
    sc = opt_0 - cvx_optval;
    opt_0 = cvx_optval;
    iterations = iterations + 1;
end
AIC_initial = [AIC_initial; 2*k - 2*opt_0];
    end
    AIC_initial_best= [AIC_initial_best,min(AIC_initial)];
end
    AIC_UV = [AIC_UV; mean(AIC_initial_best)];
end
lamda0_AIC = 0.1:0.1:1;
title0_AIC = {'lamda', 'AIC_UV'};
table0_AIC = table(lamda0_AIC.', AIC_UV,'VariableNames',title0_AIC)

%% group selection and within selection of benchmark

lamda = 0.4;
TP_UV_group = [];
TN_UV_group = [];
Ac_UV_group = [];
TP_UV_in = [];
TN_UV_in = [];
Ac_UV_in = [];
B_in = B([2 4 6 8 10],[1 3 5 7 9]);
N0_in = sum(B_in(:)==0);
N1_in = numel(B_in) - N0_in;
N_in = numel(B_in);
parfor m = 1:rep
    sc = 10;
    opt_0 = 100;
    iterations = 0;
    TP_group_initial = [];
    TN_group_initial = [];
    Ac_group_initial = [];
    TP_in_initial = [];
    TN_in_initial = [];
    Ac_in_initial = [];
    
    for i = 1: m_initial
        V1 = V{i};
while sc > 10^(-5)
    [U1,Beta0] = OptU(V1,a,r,m,n,Xt,Yt,Zt,lamda);
    
    [V1,cvx_optval]=OptV(U1,b,r,m,n,Xt,Yt,Zt,lamda);
    
    sc = opt_0 - cvx_optval;
    opt_0 = cvx_optval;
    iterations = iterations + 1;
end
B0 = U1 * V1;
TP_group_initial = [TP_group_initial; (nnz(sum(round(B0,4)&B) ~= 0)+nnz(sum(round(B0,4)&B,2) ~= 0))/((a+b)/2)];
TN_group_initial = [TN_group_initial; (nnz(sum(round(B0,4)|B) == 0)+nnz(sum(round(B0,4)|B,2) == 0))/((a+b)/2)];
Ac_group_initial = [Ac_group_initial; (nnz(sum(round(B0,4)&B) ~= 0)+nnz(sum(round(B0,4)&B,2) ~= 0)...
                 + nnz(sum(round(B0,4)|B) == 0)+nnz(sum(round(B0,4)|B,2) == 0))/(a+b)];
B0_in = B0([2 4 6 8 10],[1 3 5 7 9]);
TP_in_initial = [TP_in_initial; sum(round(B0_in,4)&B_in,'all')/N1_in];
TN_in_initial = [TN_in_initial; (N_in - sum(round(B0_in,4)|B_in,'all'))/N0_in];
Ac_in_initial = [Ac_in_initial; (sum(round(B0_in,4)&B_in,'all') + (N_in - sum(round(B0_in,4)|B_in,'all')))/N_in];
    end
    [val,idx]= max(Ac_group_initial);
    TP_UV_group = [TP_UV_group; TP_group_initial(idx)];
    TN_UV_group = [TN_UV_group; TN_group_initial(idx)];
    Ac_UV_group = [Ac_UV_group; Ac_group_initial(idx)];
     TP_UV_in = [TP_UV_in; TP_in_initial(idx)];
     TN_UV_in = [TN_UV_in; TN_in_initial(idx)];
     Ac_UV_in = [Ac_UV_in; Ac_in_initial(idx)];
end
    TP_UV_group_avg =  mean(TP_UV_group)
    TN_UV_group_avg =  mean(TN_UV_group)
    Ac_UV_group_avg =  mean(Ac_UV_group)
    TP_UV_in_avg =  mean(TP_UV_in)
    TP_UV_in_avg =  mean(TP_UV_in)
    TP_UV_in_avg =  mean(TP_UV_in)



    
    
