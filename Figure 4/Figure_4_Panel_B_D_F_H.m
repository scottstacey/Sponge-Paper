clear; clc; close all;

run_fitting = false;

Data = readtable('Data/Pooled_Chi.Bio_Experimental_ModelFitting_pp.csv');
Data = table2array(Data);

F = [0.001699, 0.002818, 0.0003610, 0.0004495];

k_fixed = [0.000322195 4000 22235 7152.095014 157337.8291 1.702385143 1.835658901 ...
           2.254209123 2.023631011 1.370890126 1.606015384 1.84551856 0.623690942 ...
           0.032968342 0.030786527 0.101382613 0.074767546 0.89477407 0.004803181 ...
           8.167464174 8565.435226 2.814544497 0.055457107 5.621638845 18.85233629 ...
           1.768328496 0.004710273 5.820217585 80527.43275 2.50692042 0.028545884];

D.pSS_02_005_1 = get_array(Data, 83:88);
D.pSS_02_005_2 = get_array(Data, 89:94);

D.pSS_02_006_1 = get_array(Data, 95:101);
D.pSS_02_006_2 = get_array(Data, 102:108);
D.pSS_02_006_3 = get_array(Data, 109:115);
D.pSS_02_006_4 = get_array(Data, 116:122);

D.pSS_02_001 = get_array(Data, 73:77);

D.PVanCC_2 = get_array(Data, 29:32);
D.PVanCC_3 = get_array(Data, 33:36);
D.PVanCC_4 = get_array(Data, 45:48);

D.PLuxB_1 = get_array(Data, 9:12);
D.PLuxB_2 = get_array(Data, 13:16);
D.PLuxB_3 = get_array(Data, 17:20);
D.PLuxB_4 = get_array(Data, 21:24);
D.PLuxB_5 = get_array(Data, 53:56);

D.PTac_2 = get_array(Data, 5:8);
D.PTac_3 = get_array(Data, 49:52);

k0 = k_fixed;
lb = k_fixed;
ub = k_fixed;

options = optimoptions( ...
    'fmincon', ...
    'FunctionTolerance', 1e-11, ...
    'StepTolerance', 1e-11, ...
    'ConstraintTolerance', 1e-11, ...
    'MaxFunctionEvaluations', 1e11, ...
    'MaxIterations', 1e11 ...
);

fixed_cost = combined_cost(k_fixed, F, D);
fprintf('Cost using fixed parameter set: %.6g\n', fixed_cost);

if run_fitting
    
    import globaloptim.*
    
    num_starts = 12;
    k0_multistart = lb + rand(num_starts, numel(lb)) .* (ub - lb);
    k0_multistart(1,:) = k0;
    
    gs = GlobalSearch( ...
        'FunctionTolerance', 1e-7, ...
        'NumTrialPoints', 11000, ...
        'NumStageOnePoints', 7000 ...
    );
    
    results = cell(num_starts, 1);
    costs = zeros(num_starts, 1);
    
    parfor i = 1:num_starts
        
        problem = createOptimProblem( ...
            'fmincon', ...
            'objective', @(k) combined_cost(k, F, D), ...
            'x0', k0_multistart(i,:), ...
            'lb', lb, ...
            'ub', ub, ...
            'nonlcon', @(k) nonlcon(k), ...
            'options', options ...
        );
        
        [k_i, cost_i] = run(gs, problem);
        
        results{i} = k_i;
        costs(i) = cost_i;
    end
    
    [cost, best_idx] = min(costs);
    k = results{best_idx};
    
    fprintf('Best fitted cost: %.6g\n', cost);
    
else
    
    k = k_fixed;
    
end

[t_005_1, y_005_1] = simulate_pSS_02_005(D.pSS_02_005_1, k, F, 10);
[t_005_2, y_005_2] = simulate_pSS_02_005(D.pSS_02_005_2, k, F, 11);

[t_006_1, y_006_1] = simulate_pSS_02_006(D.pSS_02_006_1, k, F, 6);
[t_006_2, y_006_2] = simulate_pSS_02_006(D.pSS_02_006_2, k, F, 7);
[t_006_3, y_006_3] = simulate_pSS_02_006(D.pSS_02_006_3, k, F, 8);
[t_006_4, y_006_4] = simulate_pSS_02_006(D.pSS_02_006_4, k, F, 9);

[t_001, y_001] = simulate_pSS_02_001(D.pSS_02_001, k, F, 12);

[t_PVanCC_2, y_PVanCC_2] = simulate_promoter(D.PVanCC_2, k, F, "Van");
[t_PVanCC_3, y_PVanCC_3] = simulate_promoter(D.PVanCC_3, k, F, "Van");
[t_PVanCC_4, y_PVanCC_4] = simulate_promoter(D.PVanCC_4, k, F, "Van");

[t_PLuxB_1, y_PLuxB_1] = simulate_promoter(D.PLuxB_1, k, F, "Lux");
[t_PLuxB_2, y_PLuxB_2] = simulate_promoter(D.PLuxB_2, k, F, "Lux");
[t_PLuxB_3, y_PLuxB_3] = simulate_promoter(D.PLuxB_3, k, F, "Lux");
[t_PLuxB_4, y_PLuxB_4] = simulate_promoter(D.PLuxB_4, k, F, "Lux");
[t_PLuxB_5, y_PLuxB_5] = simulate_promoter(D.PLuxB_5, k, F, "Lux");

[t_PTac_2, y_PTac_2] = simulate_promoter(D.PTac_2, k, F, "Tac");
[t_PTac_3, y_PTac_3] = simulate_promoter(D.PTac_3, k, F, "Tac");

plot_fit(1,  t_005_1, y_005_1(:,3)*k(10), D.pSS_02_005_1(:,1), D.pSS_02_005_1(:,5), 'g', ...
    'Mature GFPmut3 Protein (nM)', 'pSS-02-005 Replicate 1');

plot_fit(2,  t_005_2, y_005_2(:,3)*k(11), D.pSS_02_005_2(:,1), D.pSS_02_005_2(:,5), 'g', ...
    'Mature GFPmut3 Protein (nM)', 'pSS-02-005 Replicate 2');

plot_fit(3,  t_006_1, y_006_1(:,3)*k(6), D.pSS_02_006_1(:,1), D.pSS_02_006_1(:,5), 'g', ...
    'Mature GFPmut3 Protein (nM)', 'pSS-02-006 GFPmut3 Replicate 1');

plot_fit(4,  t_006_2, y_006_2(:,3)*k(7), D.pSS_02_006_2(:,1), D.pSS_02_006_2(:,5), 'g', ...
    'Mature GFPmut3 Protein (nM)', 'pSS-02-006 GFPmut3 Replicate 2');

plot_fit(5,  t_006_3, y_006_3(:,3)*k(8), D.pSS_02_006_3(:,1), D.pSS_02_006_3(:,5), 'g', ...
    'Mature GFPmut3 Protein (nM)', 'pSS-02-006 GFPmut3 Replicate 3');

plot_fit(6,  t_006_4, y_006_4(:,3)*k(9), D.pSS_02_006_4(:,1), D.pSS_02_006_4(:,5), 'g', ...
    'Mature GFPmut3 Protein (nM)', 'pSS-02-006 GFPmut3 Replicate 4');

plot_fit(7,  t_001, y_001(:,3)*k(12), D.pSS_02_001(:,1), D.pSS_02_001(:,4), 'g', ...
    'Mature GFPmut3 Protein (nM)', 'pSS-02-001');

plot_fit(8,  t_006_1, y_006_1(:,9)*k(6), D.pSS_02_006_1(:,1), D.pSS_02_006_1(:,6), 'r', ...
    'Mature mScarlet-I Protein (nM)', 'pSS-02-006 mScarlet-I Replicate 1');

plot_fit(9,  t_006_2, y_006_2(:,9)*k(7), D.pSS_02_006_2(:,1), D.pSS_02_006_2(:,6), 'r', ...
    'Mature mScarlet-I Protein (nM)', 'pSS-02-006 mScarlet-I Replicate 2');

plot_fit(10, t_006_3, y_006_3(:,9)*k(8), D.pSS_02_006_3(:,1), D.pSS_02_006_3(:,6), 'r', ...
    'Mature mScarlet-I Protein (nM)', 'pSS-02-006 mScarlet-I Replicate 3');

plot_fit(11, t_006_4, y_006_4(:,9)*k(9), D.pSS_02_006_4(:,1), D.pSS_02_006_4(:,6), 'r', ...
    'Mature mScarlet-I Protein (nM)', 'pSS-02-006 mScarlet-I Replicate 4');

plot_multi_fit(12, ...
    {t_PVanCC_2, t_PVanCC_3, t_PVanCC_4}, ...
    {y_PVanCC_2(:,3), y_PVanCC_3(:,3), y_PVanCC_4(:,3)}, ...
    {D.PVanCC_2, D.PVanCC_3, D.PVanCC_4}, ...
    {'k', 'b', 'g'}, ...
    'Mature GFPmut3 Protein (nM)', ...
    'pSS-01-001');

plot_multi_fit(13, ...
    {t_PLuxB_1, t_PLuxB_2, t_PLuxB_3, t_PLuxB_4, t_PLuxB_5}, ...
    {y_PLuxB_1(:,3), y_PLuxB_2(:,3), y_PLuxB_3(:,3), y_PLuxB_4(:,3), y_PLuxB_5(:,3)}, ...
    {D.PLuxB_1, D.PLuxB_2, D.PLuxB_3, D.PLuxB_4, D.PLuxB_5}, ...
    {'r', 'k', 'b', 'g', 'c'}, ...
    'Mature GFPmut3 Protein (nM)', ...
    'pSS-01-002');

plot_multi_fit(14, ...
    {t_PTac_2, t_PTac_3}, ...
    {y_PTac_2(:,3), y_PTac_3(:,3)}, ...
    {D.PTac_2, D.PTac_3}, ...
    {'k', 'b'}, ...
    'Mature GFPmut3 Protein (nM)', ...
    'pSS-01-003');

plot_mean_fit(15, ...
    {t_006_1, t_006_2, t_006_3, t_006_4}, ...
    {y_006_1(:,3)*k(6), y_006_2(:,3)*k(7), y_006_3(:,3)*k(8), y_006_4(:,3)*k(9)}, ...
    {D.pSS_02_006_1(:,1), D.pSS_02_006_2(:,1), D.pSS_02_006_3(:,1), D.pSS_02_006_4(:,1)}, ...
    {D.pSS_02_006_1(:,5), D.pSS_02_006_2(:,5), D.pSS_02_006_3(:,5), D.pSS_02_006_4(:,5)}, ...
    'g', 'Mature GFPmut3 (nM)', 'pSS-02-006 GFPmut3 Mean and SD', true);

plot_mean_fit(16, ...
    {t_006_1, t_006_2, t_006_3, t_006_4}, ...
    {y_006_1(:,9)*k(6), y_006_2(:,9)*k(7), y_006_3(:,9)*k(8), y_006_4(:,9)*k(9)}, ...
    {D.pSS_02_006_1(:,1), D.pSS_02_006_2(:,1), D.pSS_02_006_3(:,1), D.pSS_02_006_4(:,1)}, ...
    {D.pSS_02_006_1(:,6), D.pSS_02_006_2(:,6), D.pSS_02_006_3(:,6), D.pSS_02_006_4(:,6)}, ...
    'r', 'Mature mScarlet-I (nM)', 'pSS-02-006 mScarlet-I Mean and SD', true);

plot_mean_fit(17, ...
    {t_005_1, t_005_2}, ...
    {y_005_1(:,3)*k(10), y_005_2(:,3)*k(11)}, ...
    {D.pSS_02_005_1(:,1), D.pSS_02_005_2(:,1)}, ...
    {D.pSS_02_005_1(:,5), D.pSS_02_005_2(:,5)}, ...
    'g', 'Mature GFPmut3 (nM)', 'pSS-02-005 GFPmut3 Mean and SD', true);

plot_mean_fit(18, ...
    {t_PLuxB_1, t_PLuxB_2, t_PLuxB_3, t_PLuxB_4}, ...
    {y_PLuxB_1(:,3), y_PLuxB_2(:,3), y_PLuxB_3(:,3), y_PLuxB_4(:,3)}, ...
    {D.PLuxB_1(:,1), D.PLuxB_2(:,1), D.PLuxB_3(:,1), D.PLuxB_4(:,1)}, ...
    {D.PLuxB_1(:,3), D.PLuxB_2(:,3), D.PLuxB_3(:,3), D.PLuxB_4(:,3)}, ...
    'g', 'Mature GFPmut3 (nM)', 'pSS-01-002 Mean and SD (excluding rep 5)', true);


function A = get_array(Data, cols)

    A = Data(:, cols);
    A = A(~any(isnan(A), 2), :);
end


function cost = combined_cost(k, F, D)

    costs = [
        pSS_02_006_cost(k, F, D.pSS_02_006_1, 6), ...
        pSS_02_006_cost(k, F, D.pSS_02_006_2, 7), ...
        pSS_02_006_cost(k, F, D.pSS_02_006_3, 8), ...
        pSS_02_006_cost(k, F, D.pSS_02_006_4, 9), ...
        pSS_02_005_cost(k, F, D.pSS_02_005_1, 10), ...
        pSS_02_005_cost(k, F, D.pSS_02_005_2, 11), ...
        pSS_02_001_cost(k, F, D.pSS_02_001, 12), ...
        promoter_cost(k, F, D.PVanCC_2, "Van"), ...
        promoter_cost(k, F, D.PVanCC_3, "Van"), ...
        promoter_cost(k, F, D.PVanCC_4, "Van"), ...
        promoter_cost(k, F, D.PLuxB_1, "Lux"), ...
        promoter_cost(k, F, D.PLuxB_2, "Lux"), ...
        promoter_cost(k, F, D.PLuxB_3, "Lux"), ...
        promoter_cost(k, F, D.PLuxB_4, "Lux"), ...
        promoter_cost(k, F, D.PLuxB_5, "Lux"), ...
        promoter_cost(k, F, D.PTac_2, "Tac"), ...
        promoter_cost(k, F, D.PTac_3, "Tac")
    ];

    weights = [
        4/8, 4/8, 4/8, 4/8, ...
        2/2, 2/2, ...
        1, ...
        3/3, 3/3, 3/3, ...
        5/5, 5/5, 5/5, 5/5, 5/5, ...
        2/2, 2/2
    ];

    cost = sum(weights .* costs, 'omitnan') / sum(weights(~isnan(costs)));
end


function cost = pSS_02_006_cost(k, F, A, scale_idx)

    [t, y] = simulate_pSS_02_006(A, k, F, scale_idx);

    GFP_model = y(:,3) * k(scale_idx);
    GFP_exp = interp1(A(:,1), A(:,5), t, 'linear', 'extrap');

    mSc_model = y(:,9) * k(scale_idx);
    mSc_exp = interp1(A(:,1), A(:,6), t, 'linear', 'extrap');

    cost = mean([
        normal_mse(GFP_exp, GFP_model), ...
        normal_mse(mSc_exp, mSc_model)
    ], 'omitnan');
end


function cost = pSS_02_005_cost(k, F, A, scale_idx)

    [t, y] = simulate_pSS_02_005(A, k, F, scale_idx);

    model = y(:,3) * k(scale_idx);
    exp = interp1(A(:,1), A(:,5), t, 'linear', 'extrap');

    cost = normal_mse(exp, model);
end


function cost = pSS_02_001_cost(k, F, A, scale_idx)

    [t, y] = simulate_pSS_02_001(A, k, F, scale_idx);

    model = y(:,3) * k(scale_idx);
    exp = interp1(A(:,1), A(:,4), t, 'linear', 'extrap');

    cost = normal_mse(exp, model);
end


function cost = promoter_cost(k, F, A, promoter)

    [t, y] = simulate_promoter(A, k, F, promoter);

    model = y(:,3);
    exp = interp1(A(:,1), A(:,3), t, 'linear', 'extrap');

    cost = normal_mse(exp, model);
end


function error = normal_mse(exp, model)

    eps_floor = max(1, 0.01*prctile(abs(exp), 95));
    denom = max((abs(exp) + abs(model))/2, eps_floor);

    error = mean(((exp - model)./denom).^2, 'omitnan');
end


function [c, ceq] = nonlcon(k)

    c = [
        k(20) - 2000*k(19)
        k(19) - 0.02*k(20)
        k(24) - 2000*k(23)
        k(23) - 0.02*k(24)
        k(28) - 2000*k(27)
        k(27) - 0.02*k(28)
        k(14) + 1e-9 - k(16)
    ];

    ceq = [];
end


function [t, y] = simulate_pSS_02_006(A, k, F, scale_idx)

    y0 = [0, 0, A(1,5)/k(scale_idx), 0, 0, 0, 0, 0, A(1,6)/k(scale_idx)];
    y0(y0 < 0) = 0;

    [t, y] = ode15s(@(t,y) pSS_02_006_ODEs(t, y, k, F, A), [0 max(A(:,1))], y0);
end


function [t, y] = simulate_pSS_02_005(A, k, F, scale_idx)

    y0 = [0, 0, A(1,5)/k(scale_idx), 0, 0, 0, 0];
    y0(y0 < 0) = 0;

    [t, y] = ode15s(@(t,y) pSS_02_005_ODEs(t, y, k, F, A), [0 max(A(:,1))], y0);
end


function [t, y] = simulate_pSS_02_001(A, k, F, scale_idx)

    y0 = [0, 0, A(1,4)/k(scale_idx), 0, 0];
    y0(y0 < 0) = 0;

    [t, y] = ode15s(@(t,y) pSS_02_001_ODEs(t, y, k, F, A), [0 max(A(:,1))], y0);
end


function [t, y] = simulate_promoter(A, k, F, promoter)

    y0 = [
        (A(1,4)*((A(1,4)*A(1,3))/F(2)) + F(2)*((A(1,4)*A(1,3))/F(2)))/k(31), ...
        ((A(1,4)*A(1,3))/F(2)), ...
        A(1,3)
    ];

    y0(y0 < 0) = 0;

    [t, y] = ode15s(@(t,y) promoter_ODEs(t, y, k, F, A, promoter), [0 max(A(:,1))], y0);
end


function plot_fit(fig_num, t_model, y_model, t_exp, y_exp, colour, ylab, fig_title)

    figure(fig_num); clf;
    plot(t_model, y_model, [colour ':'], t_exp, y_exp, colour, 'LineWidth', 2);
    xlabel('Time (seconds)');
    ylabel(ylab);
    title(fig_title);
    legend('Model', 'Experimental Data', 'Location', 'best');
    grid off;
end


function plot_multi_fit(fig_num, t_models, y_models, exp_arrays, colours, ylab, fig_title)

    figure(fig_num); clf; hold on;

    leg = cell(1, 2*numel(t_models));

    for i = 1:numel(t_models)
        plot(t_models{i}, y_models{i}, [colours{i} ':'], 'LineWidth', 2);
        plot(exp_arrays{i}(:,1), exp_arrays{i}(:,3), colours{i}, 'LineWidth', 2);

        leg{2*i - 1} = sprintf('Model %d', i);
        leg{2*i} = sprintf('Experimental Data %d', i);
    end

    xlabel('Time (seconds)');
    ylabel(ylab);
    title(fig_title);
    legend(leg, 'Location', 'best');
    grid off;
    hold off;
end


function plot_mean_fit(fig_num, t_models, y_models, t_exps, y_exps, colour, ylab, fig_title, show_sd)

    t_end = min(cellfun(@(x) x(end), t_models));
    t_common = linspace(0, t_end, 200);

    model_all = nan(numel(t_models), numel(t_common));
    exp_all = nan(numel(t_exps), numel(t_common));

    for i = 1:numel(t_models)
        model_all(i,:) = interp1(t_models{i}, y_models{i}, t_common, 'pchip', 'extrap');
        exp_all(i,:) = interp1(t_exps{i}, y_exps{i}, t_common, 'pchip', 'extrap');
    end

    model_mean = mean(model_all, 1, 'omitnan');
    exp_mean = mean(exp_all, 1, 'omitnan');

    figure(fig_num); clf; hold on;

    if show_sd
        model_std = std(model_all, 0, 1, 'omitnan');
        exp_std = std(exp_all, 0, 1, 'omitnan');

        fill([t_common, fliplr(t_common)], ...
             [max(exp_mean + exp_std, 0), fliplr(max(exp_mean - exp_std, 0))], ...
             colour, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

        fill([t_common, fliplr(t_common)], ...
             [max(model_mean + model_std, 0), fliplr(max(model_mean - model_std, 0))], ...
             'k', 'FaceAlpha', 0.10, 'EdgeColor', 'none');
    end

    plot(t_common, exp_mean, colour, 'LineWidth', 2);
    plot(t_common, model_mean, 'k--', 'LineWidth', 2);

    xlabel('Time (seconds)');
    ylabel(ylab);
    title(fig_title);

    if show_sd
        legend('Exp \pm SD', 'Model \pm SD', 'Exp mean', 'Model mean', 'Location', 'best');
    else
        legend('Experimental mean', 'Model mean', 'Location', 'best');
    end

    grid off;
    hold off;
end


function dy = pSS_02_006_ODEs(t, y, k, F, A)

    mu = get_input(A(:,1), A(:,7), t, 'nearest');
    Van = get_input(A(:,1), A(:,2), t, 'previous');
    OC6 = get_input(A(:,1), A(:,3), t, 'previous');
    IPTG = get_input(A(:,1), A(:,4), t, 'previous');

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
    denom = max(denom, 1e-12);

    GFP_t = k(19) + (k(20)-k(19))*Van^k(22)/(Van^k(22)+k(21)^k(22));
    sRNA_t = k(23) + (k(24)-k(23))*OC6^k(26)/(OC6^k(26)+k(25)^k(26));
    spRNA_t = k(27) + (k(28)-k(27))*IPTG^k(30)/(IPTG^k(30)+k(29)^k(30));

    mRNA_term = k(1)*y(1)*y(4)*k(2)/(k(4)*denom);
    spRNA_term = k(13)*y(5)*y(4)*k(2)/(k(5)*denom);

    dy = zeros(9,1);

    dy(1) = GFP_t - (mu + F(1))*y(1) - mRNA_term;
    dy(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dy(3) = F(2)*y(2) - mu*y(3);

    dy(4) = sRNA_t ...
          - (mu + F(3))*y(4) ...
          - mRNA_term ...
          - spRNA_term ...
          + k(18)*k(17)*y(6);

    dy(5) = spRNA_t - (mu + k(14))*y(5) - spRNA_term;

    dy(6) = mRNA_term - (k(17)+mu)*y(6);
    dy(7) = spRNA_term - (k(17)+mu)*y(7);

    dy(8) = k(15)*y(5) + k(15)*y(7) - mu*y(8) - F(4)*y(8);
    dy(9) = F(4)*y(8) - mu*y(9);
end


function dy = pSS_02_005_ODEs(t, y, k, F, A)

    mu = get_input(A(:,1), A(:,6), t, 'nearest');
    Van = get_input(A(:,1), A(:,2), t, 'previous');
    OC6 = get_input(A(:,1), A(:,3), t, 'previous');
    IPTG = get_input(A(:,1), A(:,4), t, 'previous');

    denom = 1 + (y(4)*y(1))/k(4) + (y(4)*y(5))/k(5);
    denom = max(denom, 1e-12);

    GFP_t = k(19) + (k(20)-k(19))*Van^k(22)/(Van^k(22)+k(21)^k(22));
    sRNA_t = k(23) + (k(24)-k(23))*OC6^k(26)/(OC6^k(26)+k(25)^k(26));
    spRNA_t = k(27) + (k(28)-k(27))*IPTG^k(30)/(IPTG^k(30)+k(29)^k(30));

    mRNA_term = k(1)*y(1)*y(4)*k(2)/(k(4)*denom);
    spRNA_term = k(13)*y(5)*y(4)*k(2)/(k(5)*denom);

    dy = zeros(7,1);

    dy(1) = GFP_t - (mu + F(1))*y(1) - mRNA_term;
    dy(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dy(3) = F(2)*y(2) - mu*y(3);

    dy(4) = sRNA_t ...
          - (mu + F(3))*y(4) ...
          - mRNA_term ...
          - spRNA_term ...
          + k(18)*k(17)*y(6);

    dy(5) = spRNA_t - (mu + k(16))*y(5) - spRNA_term;

    dy(6) = mRNA_term - (k(17)+mu)*y(6);
    dy(7) = spRNA_term - (k(17)+mu)*y(7);
end


function dy = pSS_02_001_ODEs(t, y, k, F, A)

    mu = get_input(A(:,1), A(:,5), t, 'nearest');
    Van = get_input(A(:,1), A(:,2), t, 'previous');
    OC6 = get_input(A(:,1), A(:,3), t, 'previous');

    denom = 1 + (y(4)*y(1))/k(4);
    denom = max(denom, 1e-12);

    GFP_t = k(19) + (k(20)-k(19))*Van^k(22)/(Van^k(22)+k(21)^k(22));
    sRNA_t = k(23) + (k(24)-k(23))*OC6^k(26)/(OC6^k(26)+k(25)^k(26));
    mRNA_term = k(1)*y(1)*y(4)*k(2)/(k(4)*denom);

    dy = zeros(5,1);

    dy(1) = GFP_t - (mu + F(1))*y(1) - mRNA_term;
    dy(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dy(3) = F(2)*y(2) - mu*y(3);
    dy(4) = sRNA_t - (mu + F(3))*y(4) - mRNA_term + k(18)*k(17)*y(5);
    dy(5) = mRNA_term - (mu + k(17))*y(5);
end


function dy = promoter_ODEs(t, y, k, F, A, promoter)

    mu = get_input(A(:,1), A(:,4), t, 'nearest');
    inducer = get_input(A(:,1), A(:,2), t, 'previous');

    switch promoter
        case "Van"
            tx = k(19) + (k(20)-k(19))*inducer^k(22)/(inducer^k(22)+k(21)^k(22));
        case "Lux"
            tx = k(23) + (k(24)-k(23))*inducer^k(26)/(inducer^k(26)+k(25)^k(26));
        case "Tac"
            tx = k(27) + (k(28)-k(27))*inducer^k(30)/(inducer^k(30)+k(29)^k(30));
    end

    dy = zeros(3,1);

    dy(1) = tx - (mu + F(1))*y(1);
    dy(2) = k(31)*y(1) - (mu + F(2))*y(2);
    dy(3) = F(2)*y(2) - mu*y(3);
end


function value = get_input(t_exp, x_exp, t_model, method)

    value = interp1(t_exp, x_exp, t_model, method, 'extrap');
end