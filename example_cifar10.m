clear all
close all
clc

%% Prepare data
label_names = load('batches.meta.mat', 'label_names');
label_names = label_names.label_names;
classes = {'cat','deer','frog','dog'};
class_index = [4, 5, 6, 7] - 1; %! hard-coded
N = length(class_index);
load('data_batch_1.mat', 'labels', 'data');

n_pxl_x = 32;
n_pxl_y = 32;
n = n_pxl_x*n_pxl_y*3;

min_m_class = 100;
max_m_class = 200;
m = randi([min_m_class,max_m_class],[1,N]);

A = zeros(sum(m), n);
for i = 1:N
    data_class = data(labels == class_index(i),:);
    row_offset = sum(m(1:i-1));
    A(1+row_offset:row_offset+m(i),:) = data_class(1:m(i),:);
end

%% Call the HO-GSVD (which calls the HO-CSD after padding A)
ppi = 0.1;
[U, S, V, Tau, T, taumin, taumax] = hogsvd(A, N, m, n, ppi);

% Plot the generalized singular values and the eigenvalues of T
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
    title(sprintf('HO-GSVs of A_%d',i));
end