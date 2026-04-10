%% =========================
% 0. Initialize environment
clear
clc
close all

%% =========================
% 1. Load and preprocess n-alkane data. This script requires the file: 'Peixao_11ka_alkanes.xlsx'

filename = fullfile(pwd, 'data', 'Peixao_11ka_alkanes.xlsx')

if ~isfile(filename)
    error('Input file not found. Make sure it is in the /data folder.')
end

Pei = table2array(filename); 

%% =========================
% 2. Normalize data
Pei_norm = Pei ./ sum(Pei, 2);

% Check normalization
sum_Pei_norm_filas = sum(Pei_norm,2); 
Matriz_Pei = Pei_norm';

%% =========================
% 3. Run NMF models

[W1p,H1p,stats1p] = NMF(Matriz_Pei, 1, 1e6, 'random', 1e-9, true, 500, false);
[W2p,H2p,stats2p] = NMF(Matriz_Pei, 2, 1e6, 'random', 1e-9, true, 500, false);
[W3p,H3p,stats3p] = NMF(Matriz_Pei, 3, 1e6, 'random', 1e-9, true, 500, false);
[W4p,H4p,stats4p] = NMF(Matriz_Pei, 4, 1e6, 'random', 1e-9, true, 500, false);

%% =========================
% 4. Model constraints

sum_endmembers = round(sum(W3p,1) - 1, 12);
figure;
plot(1:3, sum_endmembers, 'ko-', 'LineWidth',2, 'MarkerFaceColor','k')
yline(0,'k-')
ylabel('sum (3 endmembers or weights) minus 1')
xlabel('A. endmember #')
set(gca,'FontSize',14)

sum_weights = sum(H3p,1) - 1; 
figure;
plot(1:length(sum_weights), sum_weights, 'k.-')
yline(0,'k--')
ylabel('sum (3 weights) minus 1')
xlabel('sample #')
set(gca,'FontSize',14);
xlim([1 length(sum_weights)])  

%% =========================
% 5. PCA comparison and scree plot

[coeffP, scoreP, latentP, tsquaredP, explainedP, muP] = pca(Pei_norm); 


var_total = sum(latentP);
resid_var_pca = zeros(length(latentP),1);
for k = 1:length(latentP)
    resid_var_pca(k) = sum(latentP(k+1:end));
end

resid_var_nmf = zeros(4,1);
for nfacs = 1:4
    switch nfacs
        case 1
            Wcurr = W1p; Hcurr = H1p;
        case 2
            Wcurr = W2p; Hcurr = H2p;
        case 3
            Wcurr = W3p; Hcurr = H3p;
        case 4
            Wcurr = W4p; Hcurr = H4p;
    end
    Matriz_Pei_hat = Wcurr * Hcurr;    
    resid_var_nmf(nfacs) = sum((Matriz_Pei(:) - Matriz_Pei_hat(:)).^2); % suma de cuadrados del residuo
end

resid_var_nmf_scaled = resid_var_nmf / resid_var_nmf(1) * resid_var_pca(1);

% scree plot
figure; hold on
plot(1:length(resid_var_pca), resid_var_pca, 'ko-', 'MarkerFaceColor','w', 'LineWidth',1.5, 'MarkerSize',6)
plot(1:4, resid_var_nmf_scaled, 'bs-', 'LineWidth',1.5, 'MarkerFaceColor','w', 'MarkerSize',6)
xlabel('PC # or N factors');
ylabel('resid. variance');
legend('pca','nmf','Location','northeast');
box on
set(gca,'FontSize',13)
hold off

%% =========================
% 6. Endmember distributions

alkane_names_Peixao = P.Properties.VariableNames();

carbon_numbers = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), alkane_names_Peixao);
endmember_colors = {[1 0.65 0], [0 0.7 0], [0.4 0.7 1]}; % naranja, verde, azul, púrpura

nalk = length(carbon_numbers);
ylimit = max([W1p(:); W2p(:); W3p(:)])*1.1;

figure('Color','w','Position',[100 100 1100 850]);
sgtitle('Peixao', 'FontSize',14,'FontWeight','bold')

ncols = 3; 
nrows = 3; 

% ---- LOOP ----
for col = 1:ncols
    switch col
        case 1, n_endm = 1; W = W1p;
        case 2, n_endm = 2; W = W2p;
        case 3, n_endm = 3; W = W3p;
        
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
% 7. Temporal evolution of endmember contributions

endmember_colors = [1 0.65 0; 0 0.7 0; 0.4 0.7 1; 0.7 0 0.4]; 

Pedad = readtable('Peixao_11ka_alkanes.xlsx', 'Sheet', 'Age');
Pedad = table2array(Pedad);
edad_Peixao = Pedad;

figure('Color','w','Position',[100 100 900 800]);
sgtitle('Peixao', 'FontSize',14,'FontWeight','bold')

nrows = 1;
margen_x = (max(edad_Peixao) - min(edad_Peixao)) * 0.03; 
x_num = max(edad_Peixao) + margen_x;            

font_axis = 8; 
font_label = 10;

subplot(nrows,1,1)
weights = H3p';
weights_norm = weights ./ sum(weights,2);
h = area(edad_Peixao, weights_norm, 'LineStyle', 'none'); 
for i=1:3, h(i).FaceColor = endmember_colors(i,:); end
ylim([0 1])
ylabel('C. weights (0-1) 3 facs','FontSize',font_label)
set(gca,'FontSize',font_axis)
xlim([min(edad_Peixao) max(edad_Peixao)+margen_x])
set(gca, 'XDir','reverse')
text(max(edad_Peixao)*0.98, 0.05, '1', 'Color', endmember_colors(1,:), 'FontSize',12, 'FontWeight','bold')
text(max(edad_Peixao)*0.98, 0.32, '2', 'Color', endmember_colors(2,:), 'FontSize',12, 'FontWeight','bold')
text(max(edad_Peixao)*0.98, 0.65, '3', 'Color', endmember_colors(3,:), 'FontSize',12, 'FontWeight','bold')
xlabel('Age (modeled y BP)','FontSize',font_label)

%% =========================
% 7. Projetction of NMF results into PCA space

[coeffP, scoreP, latentP, tsquaredP, explainedP, muP] = pca(Pei_norm);  

muP = mean(Pei_norm,1); 
endmembers_NMF = W3p'; % 
endmembers_centered = endmembers_NMF - muP; 

endmembers_proj = endmembers_centered * coeffP(:,1:2); 

% NMF scores
nmf_reconstructed = (W3p * H3p)'; 
nmf_reconstructed_centered = nmf_reconstructed - muP;   
nmf_scores = nmf_reconstructed_centered * coeffP(:,1:2);

% biplot
figure('Color','w'); hold on;

% a) Scores PCA
plot(scoreP(:,1), scoreP(:,2), 'ro', 'MarkerFaceColor','none', 'LineWidth',1.5, 'DisplayName','PCA scores');

% b) Vectores
for i = 1:size(coeffP,1)
    quiver(0,0, coeffP(i,1)*0.15, coeffP(i,2)*0.15, 0, 'Color','b','LineWidth',1.2, 'MaxHeadSize',0.4)
    if exist('alkane_names_Doninhos','var')
        text(coeffP(i,1)*0.17, coeffP(i,2)*0.17, strrep(alkane_names_Peixao{i},'_','\_'), 'FontSize',10, 'Color','b');
    end
end

% c) NMF scores
plot(nmf_scores(:,1), nmf_scores(:,2), 'k+','LineWidth',1.2,'MarkerSize',8,'DisplayName','NMF scores')

% d) NMF endmembers
plot(endmembers_proj(:,1), endmembers_proj(:,2), 'bs', 'MarkerFaceColor','none', 'LineWidth',2, 'MarkerSize',10, 'DisplayName','NMF endmembers')
for k=1:3
    text(endmembers_proj(k,1), endmembers_proj(k,2), ['  ' num2str(k)],'Color','b','FontWeight','bold','FontSize',13);
end

% e)
plot([endmembers_proj(:,1); endmembers_proj(1,1)], [endmembers_proj(:,2); endmembers_proj(1,2)], 'b-', 'LineWidth',1.5);

% f)
K = convhull(nmf_scores(:,1), nmf_scores(:,2));
plot(nmf_scores(K,1), nmf_scores(K,2), 'r--','LineWidth',1.2);

% 
xlabel(['PC 1 (' num2str(explainedP(1),'%.1f') '%)']);
ylabel(['PC 2 (' num2str(explainedP(2),'%.1f') '%)']);
axis tight
box on
set(gca,'FontSize',13)
xlim([-0.40 0.40]);
ylim([-0.2 0.4]);
title('A. PCA Peixao')
grid on

% Environmental indices (ACL, CPI, Paq)

idx = readtable('Peixao_11ka_alkanes.xlsx', 'Sheet', 'Alkane_indices');
ACL_P = idx.ACL17_35;   % Ajusta el nombre de columna si es diferente
CPI_P = idx.CPI17_35; % Ajusta el nombre de columna si es diferente
Paq_P = idx.Paq;        % Ajusta el nombre de columna si es diferente

% Correlation between PC1 and PC2
ACL_P = ACL_P(:);
Paq_P = Paq_P(:);
CPI_P = CPI_P(:);
rACL_PC1 = corr(ACL_P, scoreP(:,1));
rACL_PC2 = corr(ACL_P, scoreP(:,2));
rPaq_PC1 = corr(Paq_P, scoreP(:,1));
rPaq_PC2 = corr(Paq_P, scoreP(:,2));
rCPI_PC1 = corr(CPI_P, scoreP(:,1));
rCPI_PC2 = corr(CPI_P, scoreP(:,2));

% arrows
scale = 0.09;

quiver(0,0, rACL_PC1*scale, rACL_PC2*scale, 0, 'Color', [0 .7 0], 'LineWidth',1.5)
text(rACL_PC1*scale*1.1, rACL_PC2*scale*1.1, 'ACL', 'Color',[0 .7 0], 'FontWeight','bold','FontSize',10)

quiver(0,0, rPaq_PC1*scale, rPaq_PC2*scale, 0, 'Color', [0 .7 0], 'LineWidth',1.5)
text(rPaq_PC1*scale*1.1, rPaq_PC2*scale*1.1, 'Paq', 'Color',[0 .7 0], 'FontWeight','bold','FontSize',10)

quiver(0,0, rCPI_PC1*scale, rCPI_PC2*scale, 0, 'Color', [0 .7 0], 'LineWidth',1.5)
text(rCPI_PC1*scale*1.1, rCPI_PC2*scale*1.1, 'CPI', 'Color',[0 .7 0], 'FontWeight','bold','FontSize',10)

%% =========================
% 8. Hydrogen isotope analysis (δ2H)

% 1. endmember 3 vs δ2H
filename = 'Peixao_11ka_alkanes.xlsx';
opts = detectImportOptions(filename, 'Sheet', 'alkane_Hydrogen_isotopes');
P_H = readtable(filename, opts);

disp(P_H.Properties.VariableNames)

frac_acuatico3 = H3p(3,:) ./ sum(H3p, 1); 

C25 = P_H.D_C25;
C27 = P_H.D_C27;
C29 = P_H.D_C29;
C31 = P_H.D_C31;
C33 = P_H.D_C33;

homologos = {'D_C25', 'D_C27', 'D_C29', 'D_C31', 'D_C33'};
colores = {'bd','g^','mo', 'rv', 'cp'};

figure; hold on;
for i = 1:length(homologos)
    y = P_H.(homologos{i});
    plot(frac_acuatico3, y, colores{i}, 'MarkerFaceColor', colores{i}(1));
end

xlabel('endmember 3 / sum(all)');
ylabel('\delta^2H (‰ VSMOW)');
legend(homologos, 'Location', 'best');
grid on;
hold off;

% 2. Paq vs δ2H values

figure; hold on;
for i = 1:length(homologos)
    y = P_H.(homologos{i});
    plot(Paq_P , y, colores{i}, 'MarkerFaceColor', colores{i}(1));
end
xlabel('P_{aq}');
ylabel('\delta^2H (‰ VSMOW)');
legend({'C25', 'C27', 'C29', 'C31', 'C33'}, 'Location', 'best');
grid on;
hold off;

% 3. NMF Endmember 3 vs Paq & Regresion

frac_acuatico3 = H3p(3,:) ./ sum(H3p, 1);

figure;
scatter(Paq_P, frac_acuatico3, 60, 'ko', 'filled');
xlabel('P_{aq}');
ylabel('endmember 3/(1+2+3)');

p = polyfit(Paq_P, frac_acuatico3, 1);
xfit = linspace(min(Paq_P), max(Paq_P), 100);
yfit = polyval(p, xfit);
hold on;
plot(xfit, yfit, 'k-','LineWidth',1.5);

R = corrcoef(Paq_P, frac_acuatico3);
R2 = R(1,2)^2;

eqn = sprintf('y=%.2fx%+.2f\nR^2=%.2f', p(1), p(2), R2);
text(mean(Paq_P), max(frac_acuatico3)*0.95, eqn, 'FontSize', 12);

grid on;
hold off;

%% =========================
% 9. Carbon isotope analysis (δ13C)

% 1. δ13C data
opts = detectImportOptions(filename, 'Sheet', 'alkane_Carbon_isotopes');
C_Peixao = readtable(filename, opts);

disp(C_Peixao.Properties.VariableNames)
nombres_homologos = {'x13C_C25', 'x13_C27', 'x13_C29', 'x13_C31', 'x13_C33'};
C_matrix_Peixao = table2array(C_Peixao(:, nombres_homologos)); % [n_muestras x 5]
   
% 2. δ13C vs ACL
% a) endmember 1 (grasses)

frac_EM1 = H3p(1,:) ./ sum(H3p,1);
frac_EM1 = frac_EM1'; 

ACL = ACL_P;

figure; hold on;
scatter(ACL, frac_EM1, 50, 'o', 'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor','k');

p = polyfit(ACL, frac_EM1, 1);
xfit = linspace(min(ACL), max(ACL), 100);
yfit = polyval(p, xfit);
plot(xfit, yfit, 'k-', 'LineWidth', 1.5);

R = corrcoef(ACL, frac_EM1);
R2 = R(1,2)^2;

text(mean(ACL), max(frac_EM1)*0.8, ...
    sprintf('y = %.2fx %+ .2f\nR² = %.2f', p(1), p(2), R2), ...
    'FontSize', 12, 'BackgroundColor','w');

xlabel('ACL_{17–35}');
ylabel('Fraction EM1 (grasses)');
set(gca, 'FontSize', 14);
grid on;

% b) endmember 2 (woody)

frac_EM2 = H3p(2,:) ./ sum(H3p,1);
frac_EM2 = frac_EM2';   % vector columna

figure; hold on;
scatter(ACL, frac_EM2, 50, '^', 'MarkerFaceColor', [0.2 0.6 0.2], 'MarkerEdgeColor','k');

p = polyfit(ACL, frac_EM2, 1);
xfit = linspace(min(ACL), max(ACL), 100);
yfit = polyval(p, xfit);
plot(xfit, yfit, 'k-', 'LineWidth', 1.5);

R = corrcoef(ACL, frac_EM2);
R2 = R(1,2)^2;

text(mean(ACL), max(frac_EM2)*0.8, ...
    sprintf('y = %.2fx %+ .2f\nR² = %.2f', p(1), p(2), R2), ...
    'FontSize', 12, 'BackgroundColor','w');

xlabel('ACL_{17–35}');
ylabel('Fraction EM2 (woody vegetation)');
set(gca, 'FontSize', 14);
grid on;

%% =========================
% 9. Results δ2Hterr from NMF EM1+EM2

N = size(H3p,2);        
Nh = size(W3p,1);        
Nfuentes = 3;           

delta2H_modelado = nan(N, Nfuentes);

for i = 1:N
    for f = 1:Nfuentes
        
        num = 0; denom = 0;
        
        for j = 1:Nh
            
            num_total = sum(W3p(j,:) .* H3p(:,i)');
            
            if isnan(P_matrix(i,j)) || num_total == 0
                continue
            end
            
            w_rel = (W3p(j,f) * H3p(f,i)) / num_total;
            
            num = num + w_rel * P_matrix(i,j);
            denom = denom + w_rel;
        end
        
        if denom > 0
            delta2H_modelado(i,f) = num / denom;
        end
    end
end

[edad_sorted, idx_sort] = sort(edad_Peixao);
delta2H_sorted = delta2H_modelado(idx_sort, :);
P_matrix_sorted = P_matrix(idx_sort, :);

H_sorted = H3p(:,idx_sort);  % reordenar pesos con las edades
w12 = H_sorted(1,:) + H_sorted(2,:);  % peso total terrestre
frac1 = H_sorted(1,:) ./ (w12+eps);
frac2 = H_sorted(2,:) ./ (w12+eps);

delta2H_terr = frac1' .* delta2H_sorted(:,1) + frac2' .* delta2H_sorted(:,2);

window = 5;
terr_smooth = movmean(delta2H_terr, window, 'omitnan');

% Plot
figure; hold on;
plot(edad_sorted, terr_smooth, 's-', 'Color',[0.3 0.3 0.8], ...
    'MarkerFaceColor',[0.3 0.3 0.8],'LineWidth',2,'DisplayName','Terrestrial (EM1+EM2)');

nC = [25 27 29 31 33];
cmap = lines(size(P_matrix,2));
for j = 1:size(P_matrix,2)
    plot(edad_sorted, P_matrix_sorted(:,j), '.', 'MarkerSize', 13, ...
        'Color', cmap(j,:), 'DisplayName',['C' num2str(nC(j))]);
end

legend('Location','northeastoutside');
xlabel('Age (cal yr BP)');
ylabel('\delta^2H (‰ VSMOW)');
title('Hydrogen isotope NMF terrestrial + measured values');

set(gca,'XDir','reverse','FontSize',14);
xlim([-1000 12000]);
ylim([-210 -140]);
grid on; hold off;

%% =========================
% 10. Results δ13Cterr from NMF EM1+EM2

delta13C_modelado = nan(N, Nfuentes); % [muestras x fuentes]

% --- LOOP PRINCIPAL ---
for i = 1:N        % Loop sobre muestras
    for f = 1:Nfuentes    % Loop sobre fuentes (end-members)
        num = 0; denom = 0;
        for j = 1:Nh     % Loop sobre homólogos d13C
            Wj = W3p(j, f);       % perfil fuente f en homólogo j
            Hj = H3p(f, i);       % peso de fuente f en muestra i

            num_total = sum(W3p(j,:) .* H3p(:,i)');
            if isnan(C_matrix_Peixao(i,j)) || num_total == 0
                continue
            end
            w_rel = (Wj * Hj) / num_total;
            num = num + w_rel * C_matrix_Peixao(i,j);
            denom = denom + w_rel;
        end
        if denom > 0
            delta13C_modelado(i,f) = num / denom;
        end
    end
end

[edad_sorted, idx_sort] = sort(edad_Peixao);
delta13C_sorted = delta13C_modelado(idx_sort, :);
C_matrix_sorted = C_matrix_Peixao(idx_sort, :);

H_sorted = H3p(:, idx_sort);
w12 = H_sorted(1,:) + H_sorted(2,:);
frac1 = H_sorted(1,:) ./ (w12 + eps);
frac2 = H_sorted(2,:) ./ (w12 + eps);

delta13C_terr = frac1' .* delta13C_sorted(:,1) + ...
                frac2' .* delta13C_sorted(:,2);
)
window = 5;
terr_smooth_c = movmean(delta13C_terr, window, 'omitnan');

figure; hold on;

plot(edad_sorted, terr_smooth_c, 's-', ...
    'Color', [0.4 0 0.7], 'MarkerFaceColor', [0.4 0 0.7], 'LineWidth', 2, ...
    'DisplayName','Terrestrial (EM1+EM2)');

nC = [25 27 29 31 33];
colors_points = lines(Nh); % paleta de colores automática
for j = 1:Nh
    plot(edad_sorted, C_matrix_sorted(:,j), '.', ...
        'MarkerSize', 12, 'Color', colors_points(j,:), ...
        'DisplayName',['C' num2str(nC(j))]);
end

legend('Location','best');
xlabel('Age (years BP)');
ylabel('\delta^{13}C (‰ VPDB)');
title('delta^{13}C NMF terrestrial + measured values');

xlim([-1000 12000]);
ylim([-40 -28]); % ajusta según tu rango real

grid on;
hold off;

%% =========================
% 11. Export isotope results EM1+EM2 from NMF

U = table(edad_sorted(:), terr_smooth(:), terr_smooth_c(:),...
    'VariableNames', {'Edad_calBP', 'Deuterio Terrestre_NMF', 'Carbono_Terrestre_NMF'});

writetable(U, fullfile(pwd, 'Isotope_results_alkanes_by_NMF.xlsx'))

%% =========================
% 12. Model performance (NMF): SSE, R2 & RELATVE IMPROVEMENT BETWEEN MODELS

% Errores finales
E2 = stats2p.E(end);
E3 = stats3p.E(end);

% Etrue:
Etrue2 = stats2p.Etrue(end);
Etrue3 = stats3p.Etrue(end);

% Explained variance
EV2 = stats2p.EV;  if numel(EV2)>1, EV2 = EV2(end); end
EV3 = stats3p.EV;  if numel(EV3)>1, EV3 = EV3(end); end

fprintf('K=2: E(end)=%.6g | Etrue(end)=%.6g | EV=%.6g | stop=%s\n', E2, Etrue2, EV2, string(stats2p.stop_reason));
fprintf('K=3: E(end)=%.6g | Etrue(end)=%.6g | EV=%.6g | stop=%s\n', E3, Etrue3, EV3, string(stats3p.stop_reason));

fprintf('Relative improvement in E from K=2 to K=3: %.2f %%\n', 100*(E2 - E3)/E2);

% Gain in explained variance
fprintf('Gain in EV from K=2 to K=3: %.6g\n', EV3 - EV2);

% For K=3
Efac3 = stats3p.Efac_cum;
if numel(Efac3)>1
    fprintf('K=3 cumulative explained by factors (Efac_cum):\n');
    disp(Efac3(:)');
end

% For K=2
Efac2 = stats2p.Efac_cum;
if numel(Efac2)>1
    fprintf('K=2 cumulative explained by factors (Efac_cum):\n');
    disp(Efac2(:)');
end

% SSE
X = Matriz_Pei;
Xhat2 = W2p * H2p;
Xhat3 = W3p * H3p;

% SSE
SSE2 = norm(X - Xhat2, 'fro')^2;
SSE3 = norm(X - Xhat3, 'fro')^2;

% SST (total)
Xm = X - mean(X,2);
SST = norm(Xm, 'fro')^2;

R2_2 = 1 - SSE2/SST;
R2_3 = 1 - SSE3/SST;

fprintf('K=2: SSE=%.6g, R2=%.6f\n', SSE2, R2_2);
fprintf('K=3: SSE=%.6g, R2=%.6f\n', SSE3, R2_3);
fprintf('Gain in R2 (K=3 - K=2): %.6f\n', R2_3 - R2_2);

%% =========================
% 13. Modern vegetation comparison

% 1. Read modern vegetation dataset

Tmod = readtable(filename, 'Sheet', 'Modern_vegetation');

% Homologue labels (must match your NMF compound order!)
alkLabels = arrayfun(@(c) sprintf('C%d',c), 17:35, 'UniformOutput', false);

% Extract matrix of modern abundances (nModern x nComp)
Xmod_raw = Tmod{:, alkLabels};    % 10 x 19
groups   = string(Tmod.("EcologicalForm"));
samples  = string(Tmod.Sample);

% Normalize each modern sample to sum 1 (compositional)
Xmod = Xmod_raw ./ sum(Xmod_raw, 2);

% 2. Make sure W is normalized consistently (optional but recommended)
% W: (nComp x K) from your NMF (sediment endmembers)
W = W3p; 
H = H3p;

Wn = W ./ sum(W,1);
K = size(Wn,2);   % K = 2

% Sanity check: dimensions
[nModern, nComp] = size(Xmod);
assert(size(Wn,1) == nComp, 'W rows must match C17–C35 (19 compounds).');

% 3. Pearson and Spearman correlation modern vs endmembers

nModern = size(Xmod,1);
K = size(Wn,2);

Rpear = nan(nModern, K);
Ppear = nan(nModern, K);
Rspear = nan(nModern, K);
Pspear = nan(nModern, K);

for i = 1:nModern
    x = Xmod(i,:)';
    if any(isnan(x)), continue; end

    for k = 1:K
        y = Wn(:,k);

        % Pearson
        [r,p] = corr(x, y, 'Type','Pearson', 'Rows','complete');
        Rpear(i,k) = r; Ppear(i,k) = p;

        % Spearman (rank-based, more robust)
        [r,p] = corr(x, y, 'Type','Spearman', 'Rows','complete');
        Rspear(i,k) = r; Pspear(i,k) = p;
    end
end

% Put into tables
K = size(Rpear,2);
EMnames = arrayfun(@(k) sprintf('EM%d',k), 1:K, 'UniformOutput', false);

Tpear    = array2table(Rpear,   'VariableNames', strcat(EMnames,'_rP'));  % Pearson r
Tpear_p  = array2table(Ppear,   'VariableNames', strcat(EMnames,'_pP'));  % Pearson p
Tspear   = array2table(Rspear,  'VariableNames', strcat(EMnames,'_rS'));  % Spearman rho
Tspear_p = array2table(Pspear,  'VariableNames', strcat(EMnames,'_pS'));  % Spearman p

CorrResults = [table(samples, groups), Tpear, Tpear_p, Tspear, Tspear_p];
disp(CorrResults);

writetable(CorrResults, 'Modern_vs_Endmembers_Correlations_K3.xlsx');

% 4. NNLS

A = nan(nModern, K);      % mixing proportions per modern sample
fitErr = nan(nModern, 1); % reconstruction error

for j = 1:nModern
    x = Xmod(j,:)';                 % 19x1
    a = lsqnonneg(Wn, x);           % Kx1, a>=0
    if sum(a) > 0
        a = a ./ sum(a);            % convert to proportions that sum to 1
    end
    A(j,:) = a';
    fitErr(j) = norm(x - Wn*a);     % smaller = better match
end

% Put results in a table
A_tbl = array2table(A, 'VariableNames', arrayfun(@(k) sprintf('EM%d',k), 1:K, 'UniformOutput', false));
R = [table(samples, groups, fitErr), A_tbl];
disp(R);

% 5. Bray-Curtis similarity
% Bray-Curtis distance: d = sum(|x-y|)/sum(x+y) ; similarity = 1-d
BCsim = nan(nModern, K);
for j = 1:nModern
    x = Xmod(j,:);
    for k = 1:K
        y = Wn(:,k)';
        d = sum(abs(x - y)) / sum(x + y);
        BCsim(j,k) = 1 - d;
    end
end

BC_tbl = array2table(BCsim, 'VariableNames', A_tbl.Properties.VariableNames);
BC_res = [table(samples, groups), BC_tbl];
disp(BC_res);

%% =========================
% 14. Save workspace

save('your_workspace.mat')

