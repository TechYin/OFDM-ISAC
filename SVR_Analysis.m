%% Analysis for Singular Value Ratios (SVRs) of the optimized pilot-pattern dictionaries
% Totally 5 pilot-pattern dictionaries, each has Q=100 pilot-patterns in it.

clear; clc; close all; rng(1);
L = 128;
K = 512;
N = 64;
DFT = zeros(N,K);   % N x K DFT Matirx
for n = 0:1:N-1
    for k = 0:1:K-1
        DFT(n+1,k+1) = exp(-1j*2*pi*k*n/K);
    end
end
load PilotPatterns.mat SymbolInPilot PilotPattern_Dictionaries
%% 1 - SVR Analysis - Unoptimized pilot-patterns after primary cleaning to random pilot-patterns
PatternsNumInitial = 10000;
PatternsNumAftClean = 8000;
pop = zeros(L,PatternsNumInitial);
for nt = 1:PatternsNumInitial
    [~,PilotIndex,~] = func_RandomPilotIndex(L,K); % Random generation of pilot-patterns
    pop(:,nt) = PilotIndex;
end
SVR_rand = zeros(1,PatternsNumInitial);
SVR_randAftClean = zeros(1,PatternsNumAftClean);
for nt = 1:PatternsNumInitial
    F_P =  DFT(:,pop(:,nt));
    P = F_P*diag(SymbolInPilot);
    [U, S, V] = svd(P);
    max_singular_value = max(diag(S));  % The maximum singular value
    min_singular_value = min(diag(S));  % The minimum singular value
    SVR_rand(nt) = max_singular_value / min_singular_value;
end
figure('Name','All SVRs of random pilot-patterns'); plot(SVR_rand); grid on; hold on
% Clean up pilot-patterns with extremely-high SVRs
RandPosition = 1:PatternsNumInitial;
SVR_Threshold = 500;
HighSVRnum = find(SVR_rand > SVR_Threshold);
RandPosition(HighSVRnum) = [];
idx = randperm(length(RandPosition),PatternsNumAftClean); 
selectPosition = RandPosition(idx);                  % 屏蔽掉SVR超过300的位置后，随机抽取的位置
plot(selectPosition,SVR_rand(selectPosition), 'ro'); % 红色圆圈标记条件数最大的前popSize个位置（用于生成初始种群）
pop = pop(:,selectPosition);

for nt = 1:PatternsNumAftClean
    F_P =  DFT(:,pop(:,nt)); 
    P = F_P*diag(SymbolInPilot);
    [U, S, V] = svd(P);
    max_singular_value = max(diag(S));  % The maximum singular value
    min_singular_value = min(diag(S));  % The minimum singular value
    SVR_randAftClean(nt) = max_singular_value / min_singular_value;
end
%% 2 - SVR Analysis - The optimized pilot-patterns


% All pilot-patterns in total 5 dictionaries are getting involved in SVR distribution analysis.
PilotPatternSave = [PilotPattern_Dictionaries{1},PilotPattern_Dictionaries{2},PilotPattern_Dictionaries{3},PilotPattern_Dictionaries{4},PilotPattern_Dictionaries{5}];
ConditionNum_opt = zeros(1,size(PilotPatternSave,2));
for nt = 1:size(PilotPatternSave,2)
    F_P =  DFT(:,PilotPatternSave(:,nt)); 
    P = F_P*diag(SymbolInPilot);
    [U, S, V] = svd(P);
    max_singular_value = max(diag(S));  % The maximum singular value
    min_singular_value = min(diag(S));  % The maximum singular value
    ConditionNum_opt(nt) = max_singular_value / min_singular_value;
end



%% 3 - 曲线拟合1 - 初始随机导频图案的SVR分布
x = SVR_randAftClean;
x = x(:);
x = x(~isnan(x));   % 去掉NaN

xmin = min(x);

% 初值（比较稳）：theta 取略小于最小值，避免 x-theta 太小
theta0 = xmin - 0.05 * range(x);
z0 = x - theta0;
z0 = z0(z0 > 0);

mu0 = mean(log(z0));
sigma0 = std(log(z0), 1);
p0 = [theta0, mu0, max(sigma0, 1e-3)];

% 负对数似然函数
nll = @(p) negloglik_3plognorm(p, x);

% 优化选项
opts = optimset('Display','off', 'MaxIter', 5e4, 'MaxFunEvals', 5e4);
pHat = fminsearch(nll, p0, opts);

thetaHat = pHat(1);
muHat    = pHat(2);
sigmaHat = abs(pHat(3)); % 防止数值上出现负值

fprintf('3-parameter lognormal MLE:\n');
fprintf('  theta (loc) = %.6f\n', thetaHat);
fprintf('  mu          = %.6f\n', muHat);
fprintf('  sigma       = %.6f\n', sigmaHat);

% 绘制 PDF：三参数对数正态
xgrid_rand = linspace(min(x), 25, 2500);
pdf_rand = arrayfun(@(t) lognpdf_shift(t, thetaHat, muHat, sigmaHat), xgrid_rand);

% 核函数拟合 - 初始随机导频图案的SVR分布
[density_rand, x_rand] = ksdensity(x);  
%% 4 - 曲线拟合2 - 优化后导频图案的SVR分布
clear x
x = ConditionNum_opt;
x = x(:);
x = x(~isnan(x));   % 去掉NaN

xmin = min(x);

% 初值（比较稳）：theta 取略小于最小值，避免 x-theta 太小
theta0 = xmin - 0.05 * range(x);
z0 = x - theta0;
z0 = z0(z0 > 0);

mu0 = mean(log(z0));
sigma0 = std(log(z0), 1);
p0 = [theta0, mu0, max(sigma0, 1e-3)];

% 负对数似然函数
nll = @(p) negloglik_3plognorm(p, x);

% 优化选项
opts = optimset('Display','off', 'MaxIter', 5e4, 'MaxFunEvals', 5e4);
pHat = fminsearch(nll, p0, opts);

thetaHat = pHat(1);
muHat    = pHat(2);
sigmaHat = abs(pHat(3)); % 防止数值上出现负值

fprintf('3-parameter lognormal MLE:\n');
fprintf('  theta (loc) = %.6f\n', thetaHat);
fprintf('  mu          = %.6f\n', muHat);
fprintf('  sigma       = %.6f\n', sigmaHat);

% 绘制 PDF：三参数对数正态
xgrid_opt = linspace(min(x), 25, 2500);
pdf_opt = arrayfun(@(t) lognpdf_shift(t, thetaHat, muHat, sigmaHat), xgrid_opt);

% 核函数拟合 - 初始随机导频图案的SVR分布
[density_opt, x_opt] = ksdensity(x);  
%% 5 - 绘图：直方图 + 核函数拟合曲线

% 1. 绘制柱状图样式
figure; BinWidth = 0.5;
        histogram(SVR_randAftClean, 'BinWidth', BinWidth,'facecolor',[224 106 68]/255,'EdgeColor','none','Normalization','pdf'); 
        hold on; 
        BinWidth = 0.5;
        histogram(ConditionNum_opt, 'BinWidth', BinWidth,'facecolor',[102 188 152]/255,'EdgeColor','none','Normalization','pdf'); 

        plot(xgrid_rand, pdf_rand, 'k--','LineWidth', 0.25);% 对数正态分布拟合
        plot(xgrid_opt, pdf_opt, 'k--','LineWidth', 0.25); % 对数正态分布拟合
        xlim([0 25]);
        xlabel('Singular Value Ratios'); ylabel('Probability');
        legend('Initial','Optimized','LogNorm Fitting');
        set(gca,'FontSize',11,'Fontname', 'Times New Roman')
        set(gcf, 'Position', [100, 50, 350, 250]); 




% 2. 绘制点样式
        BinWidth = 0.5;
        [counts, edges] = histcounts(SVR_randAftClean, 'BinWidth', BinWidth);
        centers = edges(1:end-1) + diff(edges)/2;
        prob = counts / sum(counts) / diff(edges(1:2));

figure; plot(centers, prob, 'ks', 'MarkerSize', 2, 'MarkerFaceColor', 'black');
        hold on;
        BinWidth = 0.5;
        [counts, edges] = histcounts(ConditionNum_opt, 'BinWidth', BinWidth);
        centers = edges(1:end-1) + diff(edges)/2;
        prob = counts / sum(counts) / diff(edges(1:2));

        plot(centers, prob, 'ks', 'MarkerSize', 2, 'MarkerFaceColor', 'black');
        hold on;
        plot(xgrid_rand, pdf_rand, 'k-','LineWidth', 0.5); % 对数正态分布拟合
        plot(xgrid_opt, pdf_opt, 'k-','LineWidth', 0.5);   % 对数正态分布拟合
        xlim([0 25])
        xlabel('Singular Value Ratios'); ylabel('Probability');
        legend('Optimized','Initial','LogNorm Fit')
        set(gcf, 'Position', [100, 50, 350, 250]);

%% ====== Functions ======
function val = negloglik_3plognorm(p, x)
    theta = p(1);
    mu    = p(2);
    sigma = abs(p(3));

    z = x - theta;

    % 任何 z<=0 都是非法（logpdf 不存在），给大惩罚
    if sigma <= 0 || any(z <= 0)
        val = 1e12;
        return;
    end

    % 三参数对数正态：对 z 用 lognpdf，再加上 z 的雅可比(已经在lognpdf里体现)
    % 直接用 lognpdf(z, mu, sigma) 即可
    ll = sum(log(lognpdf(z, mu, sigma)));
    val = -ll;
end

function y = lognpdf_shift(x, theta, mu, sigma)
    z = x - theta;
    if z <= 0 || sigma <= 0
        y = 0;
    else
        y = lognpdf(z, mu, sigma);
    end
end