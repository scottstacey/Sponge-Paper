clear; clc;

%% Parameters
F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k_base = [0.000322195	4000	22235	7152.095014	157337.8291	1.702385143	1.835658901 ...
          2.254209123	2.023631011	1.370890126	1.606015384	1.84551856	0.623690942 ...
          0.032968342	0.030786527	0.101382613	0.074767546	0.89477407	0.004803181 ...
          8.167464174	8565.435226	2.814544497	0.055457107	5.621638845	18.85233629 ...
          1.768328496	0.004710273	5.820217585	80527.43275	2.50692042	0.028545884];

mu = 0.0003468419;

%% Settings
GFP_transcription = 1;   % fixed (nM/s)

% spRNA transcription rates (nM/s)
spRNA_t_vals = [0, 0.03, 0.08, 1, 10];

% OC6 sweep
N = 10000;
OC6_min_pos = 1e-3;   
OC6_max     = 10000;    

OC6_pos = logspace(log10(OC6_min_pos), log10(OC6_max), N-1);
OC6_vec = [0, OC6_pos]; 

tspan = [0, 1e6];

% Initial conditions
y0 = zeros(9,1);

odeOpts = odeset( ...
    'RelTol', 1e-7, ...
    'AbsTol', 1e-10, ...
    'NonNegative', 1:9 ...
);

fixedR       = false;
annihilation = false;

%% Run sweep for each spRNA_t
k = k_base;
if annihilation
    k(18) = 0;
end

Y_out_all = NaN(numel(OC6_vec), numel(spRNA_t_vals));

for m = 1:numel(spRNA_t_vals)

    spRNA_t = spRNA_t_vals(m);

    Y_out = NaN(size(OC6_vec));

    yIC = y0;
    for i = 1:numel(OC6_vec)
        OC6 = OC6_vec(i);

        rhs = @(t,y) odesys_rc_predation_OC6(t,y,k,F,OC6,spRNA_t,GFP_transcription,mu,fixedR);

        [~, Y] = ode15s(rhs, tspan, yIC, odeOpts);
        y_ss = max(Y(end,:)', 0);

        Y_out(i) = max(0, y_ss(2) + y_ss(3));
        yIC = y_ss;
    end

    Y_out_all(:,m) = Y_out;
end

%% Shared baseline = (OC6=0, spRNA_t=0)
[~, idx0] = min(abs(spRNA_t_vals - 0)); 
Y_baseline = Y_out_all(1, idx0);

eps0 = 1e-12;
if Y_baseline < eps0
    warning('Shared baseline at (OC6=0, spRNA_t=0) is ~0 (%.3g). Using eps for denominator.', Y_baseline);
    Y_baseline = eps0;
end

FC_all = Y_out_all ./ Y_baseline;


OC6_fit = OC6_vec(2:end);          % positive OC6 only (log axis)
fitResults = table('Size',[numel(spRNA_t_vals) 7], ...
    'VariableTypes', {'double','double','double','double','double','double','double'}, ...
    'VariableNames', {'spRNA_t','ymin','ymax','dynRange','K_EC50','n','R2'});

for m = 1:numel(spRNA_t_vals)

    
    xAll = OC6_fit(:);
    yAll = FC_all(2:end, m);
    yAll = yAll(:);

    
    ok = isfinite(xAll) & isfinite(yAll);

    x = xAll(ok);
    y = yAll(ok);

    
    ymin0 = min(y);
    ymax0 = max(y);
    if ymax0 < ymin0
        tmp = ymax0; ymax0 = ymin0; ymin0 = tmp;
    end
    K0 = median(x);
    n0 = 2;

    
    dy0 = max(ymax0 - ymin0, 1e-6);
    p0  = [log(max(n0,1e-3)); log(max(K0,1e-12)); ymin0; log(dy0)];

    obj = @(p) sum( (y - hill_decreasing(x, p)).^2 );

    opts = optimset('Display','off', 'MaxIter', 5e4, 'MaxFunEvals', 5e5);
    pHat = fminsearch(obj, p0, opts);

    % Unpack
    nHat    = exp(pHat(1));
    KHat    = exp(pHat(2));
    yminHat = pHat(3);
    dyHat   = exp(pHat(4));
    ymaxHat = yminHat + dyHat;

    % R^2
    yhat  = hill_decreasing(x, pHat);
    SSres = sum((y - yhat).^2);
    SStot = sum((y - mean(y)).^2);
    R2    = 1 - SSres/max(SStot, eps);

    % Store
    fitResults.spRNA_t(m)   = spRNA_t_vals(m);
    fitResults.ymin(m)      = yminHat;
    fitResults.ymax(m)      = ymaxHat;
    fitResults.dynRange(m)  = ymaxHat - yminHat;
    fitResults.K_EC50(m)    = KHat;
    fitResults.n(m)         = nHat;
    fitResults.R2(m)        = R2;
end

disp(fitResults);

writetable(fitResults, 'OC6_HillFits_by_spRNA.csv');

function yhat = hill_decreasing(x, p)
    n    = exp(p(1));
    K    = exp(p(2));
    ymin = p(3);
    ymax = ymin + exp(p(4));    % enforce ymax>ymin
    yhat = ymin + (ymax - ymin) ./ (1 + (x./K).^n);
end

figure('Color','w');
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');

xlims_log = [OC6_min_pos, OC6_max]; % MUST be positive

% Colours + legend
colors = lines(numel(spRNA_t_vals));
leg = cell(1, numel(spRNA_t_vals));
for m = 1:numel(spRNA_t_vals)
    leg{m} = sprintf('spRNA TX Rate$\\ = %.3g\\ \\mathrm{nM/s}$', spRNA_t_vals(m));
end

nexttile; hold on;
for m = 1:numel(spRNA_t_vals)
    semilogx(OC6_vec(2:end), FC_all(2:end,m), 'LineWidth', 2, 'Color', colors(m,:));
    plot(OC6_min_pos, FC_all(1,m), 'o', 'MarkerSize', 6, 'LineWidth', 1.5, 'Color', colors(m,:), 'HandleVisibility','off');
end
hold off;
ax = gca;
ax.XScale = 'log'; % force
xlim(xlims_log);

xlabel('OC6 concentration (nM)', 'Interpreter','latex');
ylabel('Relative GFPmut3', 'Interpreter','latex');
title('sRNA Induction Dose-Response Curve');
legend(leg, 'Interpreter','latex', 'Location','best');
grid off; box off;




function dy = odesys_rc_predation_OC6(~, y, k, F, OC6, spRNA_t, G, mu, fixedR)

    % sRNA transcription driven by OC6
    sRNA_t = (k(23) + (k(24)-k(23))*((OC6^k(26))/(OC6^k(26) + k(25)^k(26))));

    if fixedR
        denom = 1;
    else
        denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
        denom = max(denom, 1e-12);
    end

    term1 = (k(1)*y(1)*y(4)*k(2)) / (k(4)*denom);
    term2 = (k(13)*y(5)*y(4)*k(2)) / (k(5)*denom);

    dy = zeros(9,1);

    dy(1) = G - (mu + F(1))*y(1) - term1;
    dy(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dy(3) = F(2)*y(2) - mu*y(3);

    dy(4) = sRNA_t ...
            - (mu + F(3))*y(4) ...
            - term1 ...
            - term2 ...
            + k(18)*k(17)*y(6);

    dy(5) = spRNA_t ...
            - (mu + k(14))*y(5) ...
            - term2;

    dy(6) = term1 - (k(17)+mu)*y(6);
    dy(7) = term2 - (k(17)+mu)*y(7);

    dy(8) = k(15)*y(5) + k(15)*y(7) - mu*y(8) - F(4)*y(8);
    dy(9) = F(4)*y(8) - mu*y(9);
end

function format_axes(ax, isNormalised)
    grid(ax, 'off');
    box(ax, 'off');
    ax.YAxis.Exponent = 0;

    if isNormalised
        ylim(ax, [0 1.05]);
        yticks(ax, 0:0.25:1);
        ytickformat(ax, '%.2f');
    else
        ytickformat(ax, '%.0f');
    end
end