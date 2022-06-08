addpath('../')
close all
clear all
clc

%% Prepare data
label_names = load('batches.meta.mat', 'label_names');
class_index = uint8([4, 5, 6, 7] - 1); %! cat,deer,dog,frog
class_index = uint8([2, 4, 9, 10] - 1); %! automobile,cat,ship,truck
classes = label_names.label_names(class_index+1);
N = length(class_index);
load('data_batch_1.mat', 'labels', 'data');
nsamples = length(labels);
max_per_class = sum(repmat(labels,[1 N])==repmat(class_index,[nsamples 1]),1);

n_pxl_x = 32;
n_pxl_y = 32;
n = n_pxl_x*n_pxl_y*3;

min_m_class = 600;
max_m_class = 900;
m = [888   890   647   892]; % randi([min_m_class,max_m_class],[1,N]);
mb = max_per_class - m;

A = zeros(sum(m), n);   % Training set
B = zeros(sum(mb), n);  % Test set
test_labels = zeros(sum(mb), 1);
for i = 1:N
    data_class = data(labels == class_index(i),:);
    A(1+sum(m(1:i-1)):sum(m(1:i-1))+m(i),:) = double(data_class(1:m(i),:))/255;
    B(1+sum(mb(1:i-1)):sum(mb(1:i-1))+mb(i),:) = double(data_class(m(i)+1:end,:))/255;
    test_labels(1+sum(mb(1:i-1)):sum(mb(1:i-1))+mb(i)) = i;
end

%% Call the HO-GSVD (which calls the HO-CSD after padding A)
[U, S, V, Q, R, Z, Tau, T, taumin, taumax, mpad] = hogsvd(A, N, m, n);

is_pad = length(mpad)>length(m);
if is_pad, classes{end+1}='padding'; end
if is_pad, N=N+1; end

%% Plot the generalized singular values and the eigenvalues of T
figure;
subplot(N+1,1,1);
plot(1:n, (diag(Tau)-taumin)./(taumax-taumin), 'k-');
title('Eigenvalues of T relative to \tau_{min} and \tau_{max}')
xlim([1 n]);
ylim([0 1]);
for i = 1 : N
    subplot(N+1, 1, i+1);
    Si = S(1+(i-1)*n : i*n, :);
    plot(1:n,diag(Si),'k-');
    xlim([1 n]);
    ylim([0 1]);
    title(sprintf('HO-GSVs of A_%d (%s)',i,classes{i}));
end

if true
    addpath('/home/idris/Documents/EngSci/Matlab/my_functions/')
    cols = {};
    data = zeros(n,N+2);
    data(:,1) = 1:n; cols{end+1} = 'n';
    data(:,2) = (diag(Tau)-taumin)./(taumax-taumin); cols{end+1} = 'tau';
    for i=1:N
        data(:,2+i) = diag(S(1+(i-1)*n : i*n, :));  cols{end+1} = classes{i};
    end
    csvwrite_with_headers(['gsvs_',classes{:},'.csv'],data,cols);
end

%% Visualize Reconstruction
j=10;
figure;
if is_pad, i_plot=length(classes)-1; else, i_plot=length(classes); end
for i=1:i_plot
    inds_Si = 1+(i-1)*n:i*n;
    inds_Ui = 1+sum(m(1:i-1)):sum(m(1:i));
    Si = S(inds_Si,:);
    Ui = U(inds_Ui,:);
    Ai = A(inds_Ui,:);
    
    inds_iso = diag(Si) >= 1-1e-12;
    Ui_iso = Ui(:,inds_iso);
    Vi_iso = V(:,inds_iso);
    Ai_iso = Ui_iso*Vi_iso';
    n_iso = size(Vi_iso,2);
    
    [~,inds] = sort(vecnorm(Vi_iso,2,1), 'descend');
    inds = inds(1:300);
    Ai_iso_maxnorm = Ui_iso(:,inds)*Vi_iso(:,inds)';
    
    subplot(4,length(classes),i);
    imshow(imrotate(abs(reshape(Ai(j,:),[n_pxl_x n_pxl_y 3])),-90))
    title(sprintf('%s',classes{i}));
    subplot(4,length(classes),i+length(classes));
    imshow(imrotate(abs(reshape(Ai_iso(j,:),[n_pxl_x n_pxl_y 3])),-90))
    title(sprintf('Isolated Subspace'));
    subplot(4,length(classes),i+2*length(classes));
    imshow(imrotate(abs(reshape(Ai_iso_maxnorm(j,:),[n_pxl_x n_pxl_y 3])),-90))
    title(sprintf('Isolated Subspace'));
end

%% Classify: Select some vectors from each isolated subspace
guesses = zeros(sum(mb),1);
score = zeros(sum(mb),N);

for i = 1:N
    Si = S(1+(i-1)*n : i*n, :);
    iso = diag(Si) >= 1-1e-12;
    V_iso = V(:,iso);
    Ui_iso = U(1+sum(m(1:i-1)):sum(m(1:i)),iso);
    n_iso = size(Vi_iso,2);
    
    % select by V_iso
    [~,inds] = sort(vecnorm(V_iso,2,1), 'descend');
    inds = inds(1:end);
    
    V_iso_inds = V_iso(:,inds);
    oracle = pinv(V_iso_inds')*V_iso_inds';
    for j=1:size(B,1) % B(j,:) = x*V_iso_inds'
        score(j,i) = norm(B(j,:)-B(j,:)*oracle);
    end
end

for j=1:sum(mb)
    guesses(j) = find(score(j,:) == min(score(j,:)));
end

sum(test_labels == guesses)/sum(mb)*100
