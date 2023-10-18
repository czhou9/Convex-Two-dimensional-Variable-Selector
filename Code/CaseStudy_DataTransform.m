%% Case study data transformation and generation
% Generate variable matrix 
n_bad = length(BadSample);
n_good = length(GoodSample);
bad_trans = cell(1,n_bad);
good_trans = cell(1,n_good);
X_bad = cell(1,n_bad);
X_good = cell(1,n_good);

for i = 1:n_bad
    bad_trans{1,i} = mean(BadSample{1,i});
    X_bad{1,i}(1,1:7) = bad_trans{1,i}(1,5:11);
    X_bad{1,i}(2,1:7) = bad_trans{1,i}(1,12:18);
    X_bad{1,i}(3,1:7) = [0 bad_trans{1,i}(1,19:24)];
    X_bad{1,i}(4,1:7) = bad_trans{1,i}(1,[25 28 31 34 37 40 43]);
    X_bad{1,i}(5,1:7) = bad_trans{1,i}(1,[26 29 32 35 38 41 44]);
    X_bad{1,i}(6,1:7) = bad_trans{1,i}(1,[27 30 33 36 39 42 45]);
    X_bad{1,i}(7,1:7) = bad_trans{1,i}(1,46:52);
    X_bad{1,i}(8,1:7) = [0 bad_trans{1,i}(1,59:64)];
    X_bad{1,i}(9,1:7) = linspace(bad_trans{1,i}(1,2), bad_trans{1,i}(1,3),7);
    X_bad{1,i}(10,1:7) = zeros(1,7);
    X_bad{1,i}(11,1:7) = nnz(BadSample{1,i}(:,1)>0) / 1500;
end
for i = 1:n_good
    good_trans{1,i} = mean(GoodSample{1,i});
    X_good{1,i}(1,1:7) = good_trans{1,i}(1,5:11);
    X_good{1,i}(2,1:7) = good_trans{1,i}(1,12:18);
    X_good{1,i}(3,1:7) = [0 good_trans{1,i}(1,19:24)];
    X_good{1,i}(4,1:7) = good_trans{1,i}(1,[25 28 31 34 37 40 43]);
    X_good{1,i}(5,1:7) = good_trans{1,i}(1,[26 29 32 35 38 41 44]);
    X_good{1,i}(6,1:7) = good_trans{1,i}(1,[27 30 33 36 39 42 45]);
    X_good{1,i}(7,1:7) = good_trans{1,i}(1,46:52);
    X_good{1,i}(8,1:7) = [0 good_trans{1,i}(1,59:64)];
    X_good{1,i}(9,1:7) = linspace(good_trans{1,i}(1,2), good_trans{1,i}(1,3),7);
    X_good{1,i}(10,1:7) = ones(1,7);
    X_good{1,i}(11,1:7) = nnz(GoodSample{1,i}(:,1)>0) / 1500;
end
% Generate the total 490 samples
sample_total = cat(3,reshape(cell2mat(X_bad),[11,7,n_bad]),reshape(cell2mat(X_good),[11,7,n_good]));
% normalize the samples
for i = 1:9
    for j = 1:7
        sample_total(i,j,:) = normalize(sample_total(i,j,:));
    end
end
sample_total(isnan(sample_total)) = 0;
            
% Select 400 samples from the total 490 samples
rndIDX = randperm(490); 
sample_sub = sample_total(:, :,rndIDX(1:400)); 

% initial value and generae Xt, Yt
a = 9;
b = 7;
rep = 100;
n = 400;
Xt = cell(n,rep);
Yt = cell(n,rep);
for j = 1:rep
         rndIDX = randperm(490); 
         sample_sub = sample_total(:, :,rndIDX(1:n)); 
     for i = 1:n
         Xt{i,j} = sample_sub(1:9,:,i);
         Yt{i,j} = sample_sub(10,1,i);
     end
end
Zt = Yt; 