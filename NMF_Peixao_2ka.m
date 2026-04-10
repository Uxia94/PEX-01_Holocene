%% Description:
% Performs non-negative matrix factorization (NMF) on n-alkane
% distributions from the Peixão dataset (Peixao_2000_alkanes). Evaluates models with
% 1–4 endmembers, compares residual variance with PCA, and
% visualizes endmember distributions.

% Requirements: External NMF function download in: Polissar, P; Karp, A. Tyler; D'Andrea, William (2025), “Mixed messages: Unmixing sedimentary molecular distributions reveals source contributions and isotopic values”, Mendeley Data, V2, doi: 10.17632/3hymgv47jv.2

%% =========================
% 0. Initialize environment
clear
clc
close all

%% =========================
% 1. Load data

data_file = 'data/Peixao_2000_alkanes.xlsx';

if ~isfile(data_file)
    error('Input file not found. Make sure it is in the /data folder.')
end

P_2000 = readtable(data_file);
Pei_2000 = table2array(P_2000);

P_2000 = readtable(data_file);
Pei_2000 = table2array(P_2000);

%% =========================
% 2. Normalize data

Pei_norm_2000 = Pei_2000 ./ sum(Pei_2000, 2);

% Check normalization
sum_Pei_norm_2000_filas = sum(Pei_norm_2000,2); 
Matriz_Pei_2000 = Pei_norm_2000';
rng(1)
%% =========================
% 3. Run NMF models

[W1pP,H1pP,stats1pP] = NMF(Matriz_Pei_2000, 1, 1e6, 'random', 1e-9, true, 500, false);
[W2pP,H2pP,stats2pP] = NMF(Matriz_Pei_2000, 2, 1e6, 'random', 1e-9, true, 500, false);
[W3pP,H3pP,stats3pP] = NMF(Matriz_Pei_2000, 3, 1e6, 'random', 1e-9, true, 500, false);
[W4pP,H4pP,stats4pP] = NMF(Matriz_Pei_2000, 4, 1e6, 'random', 1e-9, true, 500, false);

%% =========================
% 4. Model constraints

sum_endmembers = round(sum(W3pP,1) - 1, 12);

figure;
plot(1:3, sum_endmembers, 'ko-', 'LineWidth',2, 'MarkerFaceColor','k')
yline(0,'k-')
ylabel('sum (3 endmembers or weights) minus 1')
xlabel('Endmember #')
set(gca,'FontSize',14)

sum_weights = sum(H3pP,1) - 1;

figure;
plot(1:length(sum_weights), sum_weights, 'k.-')
yline(0,'k--')
ylabel('sum (3 weights) minus 1')
xlabel('Sample #')

%% =========================
% 5. PCA comparison

[~, ~, latentPp] = pca(Pei_norm_2000);

resid_var_pca = zeros(length(latentPp),1);
for k = 1:length(latentPp)
    resid_var_pca(k) = sum(latentPp(k+1:end));
end

%% =========================
% 6. NMF residual variance

resid_var_nmf = zeros(4,1);

models = {W1pP,H1pP; W2pP,H2pP; W3pP,H3pP; W4pP,H4pP};

for i = 1:4
    W = models{i,1};
    H = models{i,2};
    recon = W * H;
    resid_var_nmf(i) = sum((Matriz_Pei_2000(:) - recon(:)).^2);
end

resid_var_nmf_scaled = resid_var_nmf / resid_var_nmf(1) * resid_var_pca(1);

%% =========================
% 7. Scree plot

figure; hold on
plot(resid_var_pca, 'ko-', 'LineWidth',1.5)
plot(1:4, resid_var_nmf_scaled, 'bs-', 'LineWidth',1.5)
xlabel('PC # or N factors')
ylabel('Residual variance')
legend('PCA','NMF')
hold off

%% =========================
% 8. Endmember distributions

alkane_names_Peixao = P_2000.Properties.VariableNames();

carbon_numbers = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), alkane_names_Peixao);
endmember_colors = {[1 0.65 0], [0 0.7 0], [0.4 0.7 1]}; % naranja, verde, azul, púrpura

nalk = length(carbon_numbers);
ylimit = max([W1pP(:); W2pP(:); W3pP(:)])*1.1;

figure('Color','w','Position',[100 100 1100 850]);
sgtitle('Peixao 2000 yrs', 'FontSize',14,'FontWeight','bold')

ncols = 3; 
nrows = 3; 

% ---- LOOP  ----
for col = 1:ncols
    switch col
        case 1, n_endm = 1; W = W1pP;
        case 2, n_endm = 2; W = W2pP;
        case 3, n_endm = 3; W = W3pP;
       
    end
    for row = 1:n_endm
        subplot(nrows, ncols, (row-1)*ncols + col)
        bar(carbon_numbers, W(:,row), 0.8, 'FaceColor', endmember_colors{row}, 'EdgeColor','k');
        ylim([0 ylimit])
        set(gca, 'XTick', carbon_numbers, 'XTickLabel', carbon_numbers, 'FontSize',12)
        if col == 1
            ylabel(['endmember ' num2str(row)])
        end
        if row == 1
            switch col
                case 1, title('A. 1 endmember')
                case 2, title('B. 2 endmembers')
                case 3, title('C. 3 endmembers')
                
            end
        end
        if row == n_endm
            xlabel('\itn-alkane carbon #')
        end
        box on
    end
end

%% =========================
% 9. Save results

save('results_workspace.mat')


