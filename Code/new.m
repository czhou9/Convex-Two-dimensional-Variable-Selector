%% Select lamda based on AIC
% AIC_new = [];
% parfor m = 1:rep
%     AIC_new_lamda = [];
%     for lamda = 0:0.01:0.1
%     [B9,cvx_optval]=newOpt(a,b,m,n,Xt,Yt,Zt,lamda);
%     AIC_new_lamda =  [AIC_new_lamda; 2*(nnz(sum(round(B9,4)) ~= 0)+nnz(sum(round(B9,4),2) ~= 0)) + 2 * cvx_optval];
%     end
%     AIC_new = [AIC_new AIC_new_lamda];
% end  
% lamda_AIC = 0:0.01:0.1;
% title_AIC = {'lamda', 'AIC_new'};
% table_AIC = table(lamda_AIC.', AIC_new,'VariableNames',title_AIC)

%% group selection and within selection
Xt = Xt_balanced;
Yt = Yt_balanced;
n = 100;
wt = 10;
lamda = 0.9;
TP_group = [];
TN_group = [];
Ac_group = [];
TP_in = [];
TN_in = [];
Ac_in = [];
B_in = B([2 4 6 8 10],[1 3 5 7 9]);
N0_in = sum(B_in(:)==0);
N1_in = numel(B_in) - N0_in;
N_in = numel(B_in);
    
for m = 1:rep  
    [B9,cvx_optval]=newOpt(a,b,m,n,Xt,Yt,lamda,wt);
        % calculate group TP, TN and Accuracy for each run
        TP_group = [TP_group; (nnz(sum(round(B9,4)&B) ~= 0)+nnz(sum(round(B9,4)&B,2) ~= 0))/((a+b)/2)];
        TN_group = [TN_group; (nnz(sum(round(B9,4)|B) == 0)+nnz(sum(round(B9,4)|B,2) == 0))/((a+b)/2)];
        Ac_group = [Ac_group; (nnz(sum(round(B9,4)&B) ~= 0)+nnz(sum(round(B9,4)&B,2) ~= 0)...
                 + nnz(sum(round(B9,4)|B) == 0)+nnz(sum(round(B9,4)|B,2) == 0))/(a+b)];
        % calculate within TP, TN and Accuracy for each run
        B9_in = B9([2 4 6 8 10],[1 3 5 7 9]);
        TP_in = [TP_in; sum(round(B9_in,4)&B_in,'all')/N1_in];
        TN_in = [TN_in; (N_in - sum(round(B9_in,4)|B_in,'all'))/N0_in];
        Ac_in = [Ac_in; (sum(round(B9_in,4)&B_in,'all') + (N_in - sum(round(B9_in,4)|B_in,'all')))/N_in];
end
TP_group_avg =  mean(TP_group);
TN_group_avg =  mean(TN_group);
Ac_group_avg =  mean(Ac_group);
TP_in_avg =  mean(TP_in);
TN_in_avg =  mean(TN_in);
Ac_in_avg =  mean(Ac_in);
TP_group_sd =  std(TP_group);
TN_group_sd =  std(TN_group);
Ac_group_sd =  std(Ac_group);
TP_in_sd =  std(TP_in);
TN_in_sd =  std(TN_in);
Ac_in_sd =  std(Ac_in);
