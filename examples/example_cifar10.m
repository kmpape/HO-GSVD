addpath('../')
close all
clear all
clc

%% Prepare data
label_names = load('batches.meta.mat', 'label_names');
%class_index = uint8([1, 2, 3, 4, 9, 10] - 1); 
class_index = uint8([2, 4, 9, 10] - 1); %! automobile,cat,ship,truck
classes = label_names.label_names(class_index+1);
N = length(class_index);
load('data_batch_1.mat', 'labels', 'data');
nsamples = length(labels);
m = sum(repmat(labels,[1 N])==repmat(class_index,[nsamples 1]),1);

nx = 32;
ny = 32;
n = nx*ny*3;

A = zeros(sum(m), n);
for i = 1:N
    data_class = data(labels == class_index(i),:);
    Ai = double(data_class(1:m(i),:))/255;
    fprintf("rank(A_%d)=%d, m_%d=%d\n", i, sum(svd(Ai)>1e-14), i, m(i));
    A(1+sum(m(1:i-1)):sum(m(1:i-1))+m(i),:) = double(data_class(1:m(i),:))/255;
end

%% Call the HO-GSVD 
[U, S, V, Tau, taumin, taumax, iso_classes] = hogsvd(A, m);

%% Plot the generalized singular values and the eigenvalues of T
nplot = n;
figure;
subplot(N+1,1,1);
plot(1:nplot, (diag(Tau(1:nplot,1:nplot))-taumin)./(taumax-taumin), 'k-');
title('Eigenvalues of T relative to \tau_{min} and \tau_{max}')
xlim([1 nplot]);
ylim([0 1]);
for i = 1 : N
    subplot(N+1, 1, i+1);
    Si = S(1+(i-1)*n : i*n, :);
    plot(1:nplot,diag(Si(1:nplot,1:nplot)),'k-'); % .*(vecnorm(V,2,1)')
    xlim([1 nplot]);
    ylim([0 1]);
    title(sprintf('HO-GSVs of A_%d (%s)',i,classes{i}));
end

%% Visualise right basis vectors with largest weight for image j
ijk_triples = [20,82,203,278];

for i=1:length(ijk_triples)
    figure
    imshow(imrotate(abs(reshape(V(:,ijk_triples(i)),[nx ny 3])),-90), 'Border','tight');
end

%% Visualize Reconstruction
ij_pairs = [16,19,40,50];
figure;
for i=1:N
    j = ij_pairs(i);
    inds_Si = 1+(i-1)*n:i*n;
    inds_Ui = 1+sum(m(1:i-1)):sum(m(1:i));
    Si = S(inds_Si,:);
    Ui = U(inds_Ui,:);
    Ai = A(inds_Ui,:);

    iso_inds_i = find(iso_classes == i);
    Ui_iso = Ui(:,iso_inds_i);
    Vi_iso = V(:,iso_inds_i);
    Ai_iso = Ui_iso*Vi_iso';

    subplot(2,length(classes),i);
    imshow(imrotate(abs(reshape(Ai(j,:),[nx ny 3])),-90))
    title(sprintf('%s',classes{i}));
    subplot(2,length(classes),i+length(classes));
    imshow(imrotate(abs(reshape(Ai_iso(j,:),[nx ny 3])),-90))
    title(sprintf('Isolated Subspace'));
end

%% Create a new A that has a common subspace
Anew = zeros(size(A));
for i = 1:N
    Ai = A(1+sum(m(1:i-1)):sum(m(1:i-1))+m(i),:);    
    Ai(1,1) = 0.5;
    Ai(2:end,1) = 0.0;
    Ai(1,2:end) = 0.0;    
    Anew(1+sum(m(1:i-1)):sum(m(1:i-1))+m(i),:) = Ai;
end

%% Call the HO-GSVD on the new dataset
[U_new, S_new, V_new, Tau_new, taumin_new, taumax_new, iso_classes_new] = hogsvd(Anew, m);

%% Plot the generalized singular values and the eigenvalues of T
nplot = n;
figure;
subplot(N+1,1,1);
plot(1:nplot, (diag(Tau_new(1:nplot,1:nplot))-taumin_new)./(taumax_new-taumin_new), 'k-');
title('Eigenvalues of modified T relative to \tau_{min} and \tau_{max}')
xlim([1 nplot]);
ylim([0 1]);
for i = 1 : N
    subplot(N+1, 1, i+1);
    Si = S_new(1+(i-1)*n : i*n, :);
    plot(1:nplot,diag(Si(1:nplot,1:nplot)),'k-'); % .*(vecnorm(V,2,1)')
    xlim([1 nplot]);
    ylim([0 1]);
    title(sprintf('HO-GSVs of A_%d (%s)',i,classes{i}));
end




