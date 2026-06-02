%% 
clear; clc; close all

% start figure counter - at zero as will first be increased
fig_cntr = 0;

%% Importing global optimum toolbox
import globaloptim.*;

%% Specify which kind of fitting we are running

% select the indicator of which growth rate calulcation method we run the fitting for
% options:
%   diff mums - different maximum growth rates for each experiment
%   same mum - same maximum growth rate for all experiments
%   fixed mus - different but FIXED (i.e. no growth feedback and mu=mum always) growth rates for each experiment
mu_calc_method='diff mums';

% set cost weights for protein level and growth rate fitting
% cost_weights = [(0)=protein_cost_weight, (1)=growth_cost_weight];
cost_weights = [1, 1];

%% Input Data
Data = readtable('Data/Figure_4_Panels_KL/Pooled_Chi.Bio_Experimental_ModelFitting_pp.csv');
Data_array = table2array(Data);

% Sponge Circuit 1.3.1
SpongeCircuit_1_3_1 = double.empty;
SpongeCircuit_1_3_1(:,1) = Data_array(:,83); % Time
SpongeCircuit_1_3_1(:,2) = Data_array(:,84); % Van
SpongeCircuit_1_3_1(:,3) = Data_array(:,85); % OC6
SpongeCircuit_1_3_1(:,4) = Data_array(:,86); % IPTG
SpongeCircuit_1_3_1(:,5) = Data_array(:,87); % EGFP
SpongeCircuit_1_3_1(:,6) = Data_array(:,88); % Growth Rate
SpongeCircuit_1_3_1 = SpongeCircuit_1_3_1(~any(isnan(SpongeCircuit_1_3_1), 2), :);

% Sponge Circuit 1.3.2
SpongeCircuit_1_3_2 = double.empty;
SpongeCircuit_1_3_2(:,1) = Data_array(:,89); % Time
SpongeCircuit_1_3_2(:,2) = Data_array(:,90); % Van
SpongeCircuit_1_3_2(:,3) = Data_array(:,91); % OC6
SpongeCircuit_1_3_2(:,4) = Data_array(:,92); % IPTG
SpongeCircuit_1_3_2(:,5) = Data_array(:,93); % GFPmut3
SpongeCircuit_1_3_2(:,6) = Data_array(:,94); % Growth Rate
SpongeCircuit_1_3_2 = SpongeCircuit_1_3_2(~any(isnan(SpongeCircuit_1_3_2), 2), :);

% Sponge Circuit 1.4.1
SpongeCircuit_1_4_1 = double.empty;
SpongeCircuit_1_4_1(:,1) = Data_array(:,95); % Time1
SpongeCircuit_1_4_1(:,2) = Data_array(:,96); % Van
SpongeCircuit_1_4_1(:,3) = Data_array(:,97); % OC6
SpongeCircuit_1_4_1(:,4) = Data_array(:,98); % IPTG
SpongeCircuit_1_4_1(:,5) = Data_array(:,99); % GFPmut3
SpongeCircuit_1_4_1(:,6) = Data_array(:,100); % mScarlet-I
SpongeCircuit_1_4_1(:,7) = Data_array(:,101); % Growth Rate
SpongeCircuit_1_4_1 = SpongeCircuit_1_4_1(~any(isnan(SpongeCircuit_1_4_1), 2), :);

% Sponge Circuit 1.4.2
SpongeCircuit_1_4_2 = double.empty;
SpongeCircuit_1_4_2(:,1) = Data_array(:,102); % Time
SpongeCircuit_1_4_2(:,2) = Data_array(:,103); % Van
SpongeCircuit_1_4_2(:,3) = Data_array(:,104); % OC6
SpongeCircuit_1_4_2(:,4) = Data_array(:,105); % IPTG
SpongeCircuit_1_4_2(:,5) = Data_array(:,106); % GFPmut3
SpongeCircuit_1_4_2(:,6) = Data_array(:,107); % mScarlet-I
SpongeCircuit_1_4_2(:,7) = Data_array(:,108); % Growth Rate
SpongeCircuit_1_4_2 = SpongeCircuit_1_4_2(~any(isnan(SpongeCircuit_1_4_2), 2), :);

% Sponge Circuit 1.4.3
SpongeCircuit_1_4_3 = double.empty;
SpongeCircuit_1_4_3(:,1) = Data_array(:,109); % Time
SpongeCircuit_1_4_3(:,2) = Data_array(:,110); % Van
SpongeCircuit_1_4_3(:,3) = Data_array(:,111); % OC6
SpongeCircuit_1_4_3(:,4) = Data_array(:,112); % IPTG
SpongeCircuit_1_4_3(:,5) = Data_array(:,113); % GFPmut3
SpongeCircuit_1_4_3(:,6) = Data_array(:,114); % mScarlet-I
SpongeCircuit_1_4_3(:,7) = Data_array(:,115); % Growth Rate
SpongeCircuit_1_4_3 = SpongeCircuit_1_4_3(~any(isnan(SpongeCircuit_1_4_3), 2), :);

% Sponge Circuit 1.4.4
SpongeCircuit_1_4_4 = double.empty;
SpongeCircuit_1_4_4(:,1) = Data_array(:,116); % Time
SpongeCircuit_1_4_4(:,2) = Data_array(:,117); % Van
SpongeCircuit_1_4_4(:,3) = Data_array(:,118); % OC6
SpongeCircuit_1_4_4(:,4) = Data_array(:,119); % IPTG
SpongeCircuit_1_4_4(:,5) = Data_array(:,120); % GFPmut3
SpongeCircuit_1_4_4(:,6) = Data_array(:,121); % mScarlet-I
SpongeCircuit_1_4_4(:,7) = Data_array(:,122); % Growth Rate
SpongeCircuit_1_4_4 = SpongeCircuit_1_4_4(~any(isnan(SpongeCircuit_1_4_4), 2), :);

% sRNA  only circuit
sRNACircuit_Array_1 = double.empty;
sRNACircuit_Array_1(:,1) = Data_array(:,73); % Time
sRNACircuit_Array_1(:,2) = Data_array(:,74); % Van 
sRNACircuit_Array_1(:,3) = Data_array(:,75); % OC6 
sRNACircuit_Array_1(:,4) = Data_array(:,76); % GFPmut3
sRNACircuit_Array_1(:,5) = Data_array(:,77); % Growth Rate

% PVanCC-GFPmut3 Replicate 2 Array
PVanCC_EGFP_Array2 = double.empty;
PVanCC_EGFP_Array2(:,1) = Data_array(:,29); % Time (seconds)
PVanCC_EGFP_Array2(:,2) = Data_array(:,30); % Vanillic Acid (nM)
PVanCC_EGFP_Array2(:,3) = Data_array(:,31); % EGFP (nM)
PVanCC_EGFP_Array2(:,4) = Data_array(:,32); % Growth rate (s^-1)
PVanCC_EGFP_Array2      = rmmissing(PVanCC_EGFP_Array2);

% PVanCC-EGFP Replicate 3 Array
PVanCC_EGFP_Array3 = double.empty;
PVanCC_EGFP_Array3(:,1) = Data_array(:,33); % Time (seconds)
PVanCC_EGFP_Array3(:,2) = Data_array(:,34); % Vanillic Acid (nM)
PVanCC_EGFP_Array3(:,3) = Data_array(:,35); % EGFP (nM)
PVanCC_EGFP_Array3(:,4) = Data_array(:,36); % Growth rate (s^-1)
PVanCC_EGFP_Array3      = rmmissing(PVanCC_EGFP_Array3);

% PVanCC-EGFP Replicate 4 Array 
PVanCC_EGFP_Array4 = double.empty;
PVanCC_EGFP_Array4(:,1) = Data_array(:,45); % Time (seconds)
PVanCC_EGFP_Array4(:,2) = Data_array(:,46); % Vanillic Acid (nM)
PVanCC_EGFP_Array4(:,3) = Data_array(:,47); % EGFP (nM)
PVanCC_EGFP_Array4(:,4) = Data_array(:,48); % Growth rate (s^-1)
PVanCC_EGFP_Array4      = rmmissing(PVanCC_EGFP_Array4);

% PLuxB-EGFP Replicate 1 Array
PLuxB_EGFP_Array1 = double.empty;
PLuxB_EGFP_Array1(:,1) = Data_array(:,9); % Time (seconds)
PLuxB_EGFP_Array1(:,2) = Data_array(:,10); % OC6 (nM)
PLuxB_EGFP_Array1(:,3) = Data_array(:,11); % EGFP (nM)
PLuxB_EGFP_Array1(:,4) = Data_array(:,12); % Growth rate (s^-1)
PLuxB_EGFP_Array1      = rmmissing(PLuxB_EGFP_Array1);

% PLuxB-EGFP Replicate 2 Array
PLuxB_EGFP_Array2 = double.empty;
PLuxB_EGFP_Array2(:,1) = Data_array(:,13); % Time (seconds)
PLuxB_EGFP_Array2(:,2) = Data_array(:,14); % OC6 (nM)
PLuxB_EGFP_Array2(:,3) = Data_array(:,15); % EGFP (nM)
PLuxB_EGFP_Array2(:,4) = Data_array(:,16); % Growth rate (s^-1)
PLuxB_EGFP_Array2      = rmmissing(PLuxB_EGFP_Array2);

% PLuxB-EGFP Replicate 3 Array
PLuxB_EGFP_Array3 = double.empty;
PLuxB_EGFP_Array3(:,1) = Data_array(:,17); % Time (seconds)
PLuxB_EGFP_Array3(:,2) = Data_array(:,18); % OC6 (nM)
PLuxB_EGFP_Array3(:,3) = Data_array(:,19); % EGFP (nM)
PLuxB_EGFP_Array3(:,4) = Data_array(:,20); % Growth rate (s^-1)
PLuxB_EGFP_Array3      = rmmissing(PLuxB_EGFP_Array3);


% PLuxB-EGFP Replicate 4 Array
PLuxB_EGFP_Array4 = double.empty;
PLuxB_EGFP_Array4(:,1) = Data_array(:,21); % Time (seconds)
PLuxB_EGFP_Array4(:,2) = Data_array(:,22); % OC6 (nM)
PLuxB_EGFP_Array4(:,3) = Data_array(:,23); % EGFP (nM)
PLuxB_EGFP_Array4(:,4) = Data_array(:,24); % Growth rate (s^-1)
PLuxB_EGFP_Array4      = rmmissing(PLuxB_EGFP_Array4);

% PLuxB-EGFP Replicate 5 Array
PLuxB_EGFP_Array5 = double.empty;
PLuxB_EGFP_Array5(:,1) = Data_array(:,53); % Time (seconds)
PLuxB_EGFP_Array5(:,2) = Data_array(:,54); % OC6 (nM)
PLuxB_EGFP_Array5(:,3) = Data_array(:,55); % EGFP (nM)
PLuxB_EGFP_Array5(:,4) = Data_array(:,56); % Growth rate (s^-1)
PLuxB_EGFP_Array5      = rmmissing(PLuxB_EGFP_Array5);

% PTac-EGFP Replicate 2 Array
PTac_EGFP_Array2 = double.empty;
PTac_EGFP_Array2(:,1) = Data_array(:,5); % Time (seconds)
PTac_EGFP_Array2(:,2) = Data_array(:,6); % IPTG (nM)
PTac_EGFP_Array2(:,3) = Data_array(:,7); % EGFP (nM)
PTac_EGFP_Array2(:,4) = Data_array(:,8); % Growth rate (s^-1)
PTac_EGFP_Array2      = rmmissing(PTac_EGFP_Array2);

% PTac-EGFP Replicate 3 Array 
PTac_EGFP_Array3 = double.empty;
PTac_EGFP_Array3(:,1) = Data_array(:,49); % Time (seconds)
PTac_EGFP_Array3(:,2) = Data_array(:,50); % IPTG (nM)
PTac_EGFP_Array3(:,3) = Data_array(:,51); % EGFP (nM)
PTac_EGFP_Array3(:,4) = Data_array(:,52); % Growth rate (s^-1)
PTac_EGFP_Array3      = rmmissing(PTac_EGFP_Array3);

%% Fixed parameters from prior fits

F_fit = [0.001699,    0.002818,    0.0003610,   0.0004495];

% Parameters set 1
k_fit = [0.000322195, 4000, ...
    1, ... % irrelevant, void k_fit(3) value)
    7152.095014,157337.8291,1.702385143,1.835658901,2.254209123,2.023631011,1.370890126,1.606015384,1.84551856,0.623690942,0.032968342,0.030786527,0.101382613,0.074767546,0.89477407,0.004803181,8.167464174,8565.435226,2.814544497,0.055457107,5.621638845,18.85233629,1.768328496,0.004710273,5.820217585,80527.43275,2.50692042,0.028545884...
    ];


%% Parameters to fit (burden thresholds)

% diff mums, same mum:
% J = [J(1)=J_GFPmut3; J(2)=J_mScarlet];
% fixed mus:
% burden thresholds set to negative values and not optimised => no entries

% units of mM - and starting at around 1e5 nM
if (~strcmp(mu_calc_method,'fixed mus'))
    J0 = [1, 1];
    J_ub = J0.*100;
    J_lb = J0./100;
else
    J0=[]; J_ub=[]; J_lb=[];
end

%% Parameters to fit (maximum possible growth rates for each trajectory)

% diff mums, fixed mus:
% mums = [  mums(1)=mum_sp141;  mums(2)=mum_sp142; mums(3)=mum_sp143;  mums(4)=mum_sp144; 
%           mums(5)=mum_sp131;  mums(6)=mum_sp132; 
%           mums(7)=mum_srna; 
%           mums(8)=mum_pvancc2; mums(9)=mum_pvancc3; mums(10)=mum_pvancc4;
%           mums(11)=mum_pluxb1; mums(12)=mum_pluxb2; mums(13)=mum_pluxb3; mums(14)=mum_pluxb4; mums(15)=mum_pluxb5;
%           mums(16)=mum_ptac2; mums(17)=mum_ptac3;
% ];
% same mum
% mums = [mum];

% default initial estimate = first measured growth rate, when there's no induction yet
mums0_default = [ ...
    SpongeCircuit_1_4_1(1,7), SpongeCircuit_1_4_2(1,7),  SpongeCircuit_1_4_3(1,7), SpongeCircuit_1_4_4(1,7), ...
    SpongeCircuit_1_3_1(1,6), SpongeCircuit_1_3_2(1,6), ...
    sRNACircuit_Array_1(1,5), ...
    PVanCC_EGFP_Array2(1,4), PVanCC_EGFP_Array3(1,4), PVanCC_EGFP_Array4(1,4), ...
    PLuxB_EGFP_Array1(1,4), PLuxB_EGFP_Array2(1,4), PLuxB_EGFP_Array3(1,4), PLuxB_EGFP_Array4(1,4), PLuxB_EGFP_Array5(1,4), ...
    PTac_EGFP_Array2(1,4), PTac_EGFP_Array3(1,4) ...
    ] .* 3600; % multiply by 3600 as mum units are 1/h

if(~strcmp(mu_calc_method, 'same mum'))
    mums0 = mums0_default; % if mus different, initial estimate particular to each experiments
else
    mums0=max(mums0_default); % if mus same, take the maximum of all initial estimates as the common one
end

% same upper and lower max growth rate bounds for all reactors
mum_ub=2.5; % upper bound 2.5/h - VERY fast-growing cells, see e.g. Chure and Cremer, eLife, 2023
mum_lb=0.75; % lower bound 0.75/h - cells growing THIS slowly actually never observed
mums_ub = mum_ub.*ones(size(mums0)); % upper bound same for all reactors
mums_lb = mum_lb.*ones(size(mums0)); % lower bound same for all reactors

%% Combine parameters to fit

Jmums0 = [J0, mums0];
ub = [J_ub, mums_ub];
lb = [J_lb, mums_lb];

% until we run the fitting, treat the starting point as the estimate
Jmums_est = Jmums0;
fit_cost=0;

%% Setting up the Fmincon an ODE integration options:

options = optimoptions('fmincon', FunctionTolerance = 1e-5, TolX = 1e-5, ConstraintTolerance = 1e-5, MaxFunctionEvaluations = 1e5, MaxIterations = 1e5);
% Scott's mega-powerful options
% options = optimoptions('fmincon', FunctionTolerance = 1e-11, TolX = 1e-11, ConstraintTolerance = 1e-11, MaxFunctionEvaluations = 1e11, MaxIterations = 1e11);
ode_options = odeset(RelTol=1e-5,AbsTol=1e-5);

%% Setting up ODE time vectors:
tspan_SpongeCircuit_1_4_1 = [0 max(SpongeCircuit_1_4_1(:,1))];
tspan_SpongeCircuit_1_4_2 = [0 max(SpongeCircuit_1_4_2(:,1))];
tspan_SpongeCircuit_1_4_3 = [0 max(SpongeCircuit_1_4_3(:,1))];
tspan_SpongeCircuit_1_4_4 = [0 max(SpongeCircuit_1_4_4(:,1))];

tspan_SpongeCircuit_1_3_1 = [0 max(SpongeCircuit_1_3_1(:,1))];
tspan_SpongeCircuit_1_3_2 = [0 max(SpongeCircuit_1_3_2(:,1))];

tspan_sRNA_circuit_1 = [0 max(sRNACircuit_Array_1(:,1))];

tspan_PVanCC_EGFP_2 = [0 max(PVanCC_EGFP_Array2(:,1))];
tspan_PVanCC_EGFP_3 = [0 max(PVanCC_EGFP_Array3(:,1))];
tspan_PVanCC_EGFP_4 = [0 max(PVanCC_EGFP_Array4(:,1))];

tspan_PLuxB_EGFP_1 = [0 max(PLuxB_EGFP_Array1(:,1))];
tspan_PLuxB_EGFP_2 = [0 max(PLuxB_EGFP_Array2(:,1))];
tspan_PLuxB_EGFP_3 = [0 max(PLuxB_EGFP_Array3(:,1))];
tspan_PLuxB_EGFP_4 = [0 max(PLuxB_EGFP_Array4(:,1))];
tspan_PLuxB_EGFP_5 = [0 max(PLuxB_EGFP_Array5(:,1))];

tspan_PTac_EGFP_2 = [0 max(PTac_EGFP_Array2(:,1))];
tspan_PTac_EGFP_3 = [0 max(PTac_EGFP_Array3(:,1))];

%% Load old results

load("Data/Figure_4_Panels_KL/Qmums_est_Muf_GrPr.mat")


%% OR, Run GlobalSearch with multiple starting points

% % Generate multiple starting points within the bounds
% num_starts = 6; % Number of starting points
% Jmums0_multistart = lb + rand(num_starts, length(lb)) .* (ub - lb);
% Jmums0_multistart(1,:) = Jmums0;
% 
% % Run GlobalSearch with multiple starting points
% gs = GlobalSearch('FunctionTolerance', 1e-6, 'NumTrialPoints', 3000, 'NumStageOnePoints', 1000);
% % Scott's mega-powerful options
% % gs = GlobalSearch('FunctionTolerance', 1e-7, 'NumTrialPoints', 11000, 'NumStageOnePoints', 7000);
% 
% % Create and run optimization problems with different starting points
% results = cell(num_starts, 1);
% costs = zeros(num_starts, 1);
% parfor i = 1:num_starts
%     problem = createOptimProblem('fmincon', 'objective', ...
%         @(Jmums) Combined_CostFunction( ...
%             k_fit, F_fit, ...
%             SpongeCircuit_1_3_1, tspan_SpongeCircuit_1_3_1, ...
%             SpongeCircuit_1_3_2, tspan_SpongeCircuit_1_3_2, ...
%             SpongeCircuit_1_4_1, tspan_SpongeCircuit_1_4_1, ...
%             SpongeCircuit_1_4_2, tspan_SpongeCircuit_1_4_2, ...
%             SpongeCircuit_1_4_3, tspan_SpongeCircuit_1_4_3, ...
%             SpongeCircuit_1_4_4, tspan_SpongeCircuit_1_4_4, ...
%             sRNACircuit_Array_1, tspan_sRNA_circuit_1, ...
%             PVanCC_EGFP_Array2, tspan_PVanCC_EGFP_2, ...
%             PVanCC_EGFP_Array3, tspan_PVanCC_EGFP_3, ...
%             PVanCC_EGFP_Array4, tspan_PVanCC_EGFP_4, ...
%             PLuxB_EGFP_Array1, tspan_PLuxB_EGFP_1, ...
%             PLuxB_EGFP_Array2, tspan_PLuxB_EGFP_2, ...
%             PLuxB_EGFP_Array3, tspan_PLuxB_EGFP_3, ...
%             PLuxB_EGFP_Array4, tspan_PLuxB_EGFP_4, ...
%             PLuxB_EGFP_Array5, tspan_PLuxB_EGFP_5, ...
%             PTac_EGFP_Array2, tspan_PTac_EGFP_2, ...
%             PTac_EGFP_Array3, tspan_PTac_EGFP_3, ... % up to now, original cost function arguments
%             Jmums, ...          % parameters being fitted: GFPmut3 and mScarlet burden thresholds J; maximum growth rates for different experiments
%             mu_calc_method, ...  % growth rate calclation method for which we run the fitting
%             cost_weights, ... % weights of protein and growth rate fitting costs
%             ode_options ... % ODE integration options
%         ), ...
%         'x0', Jmums0_multistart(i, :), 'lb', lb, 'ub', ub, 'options', options);
%     [Jmums_est, cost] = run(gs, problem);
%     results{i} = Jmums_est;
%     costs(i) = cost;
% end
% 
% % Extract the best result
% [minValue, minIndex] = min(costs);
% Jmums_est = results{minIndex};
% cost = minValue;
% results = cell2table(results);
% results = table2array(results);
% 
% % save
% save('Jmums_est.mat', "Jmums_est", "costs", ...
%     'mu_calc_method','cost_weights' ... % for record-keeping, also save the growth calculation method and the weights associated with the cost
%     )

%% Find pure protein and pure growth costs for the fit

% pure protein costs
pure_protein_cost = Combined_CostFunction( ...
            k_fit, F_fit, ...
            SpongeCircuit_1_3_1, tspan_SpongeCircuit_1_3_1, ...
            SpongeCircuit_1_3_2, tspan_SpongeCircuit_1_3_2, ...
            SpongeCircuit_1_4_1, tspan_SpongeCircuit_1_4_1, ...
            SpongeCircuit_1_4_2, tspan_SpongeCircuit_1_4_2, ...
            SpongeCircuit_1_4_3, tspan_SpongeCircuit_1_4_3, ...
            SpongeCircuit_1_4_4, tspan_SpongeCircuit_1_4_4, ...
            sRNACircuit_Array_1, tspan_sRNA_circuit_1, ...
            PVanCC_EGFP_Array2, tspan_PVanCC_EGFP_2, ...
            PVanCC_EGFP_Array3, tspan_PVanCC_EGFP_3, ...
            PVanCC_EGFP_Array4, tspan_PVanCC_EGFP_4, ...
            PLuxB_EGFP_Array1, tspan_PLuxB_EGFP_1, ...
            PLuxB_EGFP_Array2, tspan_PLuxB_EGFP_2, ...
            PLuxB_EGFP_Array3, tspan_PLuxB_EGFP_3, ...
            PLuxB_EGFP_Array4, tspan_PLuxB_EGFP_4, ...
            PLuxB_EGFP_Array5, tspan_PLuxB_EGFP_5, ...
            PTac_EGFP_Array2, tspan_PTac_EGFP_2, ...
            PTac_EGFP_Array3, tspan_PTac_EGFP_3, ... % up to now, original cost function arguments
            Jmums_est, ...          % parameters being fitted: GFPmut3 and mScarlet burden thresholds J; maximum growth rates for different experiments
            mu_calc_method, ...  % growth rate calclation method for which we run the fitting
            [1, 0], ... % weights of protein and growth rate fitting costs
            ode_options ... % ODE integration options
        );

% pure growth costs
pure_growth_cost = Combined_CostFunction( ...
            k_fit, F_fit, ...
            SpongeCircuit_1_3_1, tspan_SpongeCircuit_1_3_1, ...
            SpongeCircuit_1_3_2, tspan_SpongeCircuit_1_3_2, ...
            SpongeCircuit_1_4_1, tspan_SpongeCircuit_1_4_1, ...
            SpongeCircuit_1_4_2, tspan_SpongeCircuit_1_4_2, ...
            SpongeCircuit_1_4_3, tspan_SpongeCircuit_1_4_3, ...
            SpongeCircuit_1_4_4, tspan_SpongeCircuit_1_4_4, ...
            sRNACircuit_Array_1, tspan_sRNA_circuit_1, ...
            PVanCC_EGFP_Array2, tspan_PVanCC_EGFP_2, ...
            PVanCC_EGFP_Array3, tspan_PVanCC_EGFP_3, ...
            PVanCC_EGFP_Array4, tspan_PVanCC_EGFP_4, ...
            PLuxB_EGFP_Array1, tspan_PLuxB_EGFP_1, ...
            PLuxB_EGFP_Array2, tspan_PLuxB_EGFP_2, ...
            PLuxB_EGFP_Array3, tspan_PLuxB_EGFP_3, ...
            PLuxB_EGFP_Array4, tspan_PLuxB_EGFP_4, ...
            PLuxB_EGFP_Array5, tspan_PLuxB_EGFP_5, ...
            PTac_EGFP_Array2, tspan_PTac_EGFP_2, ...
            PTac_EGFP_Array3, tspan_PTac_EGFP_3, ... % up to now, original cost function arguments
            Jmums_est, ...          % parameters being fitted: GFPmut3 and mScarlet burden thresholds J; maximum growth rates for different experiments
            mu_calc_method, ...  % growth rate calclation method for which we run the fitting
            [0, 1], ... % weights of protein and growth rate fitting costs
            ode_options ... % ODE integration options
        );

%% Simulating ODEs after optimisation - setting initial conditions
y0_SpongeCircuit_1_4_1 = [0, 0, SpongeCircuit_1_4_1(1,5)/k_fit(6), 0, 0, 0, 0, 0, SpongeCircuit_1_4_1(1,6)/k_fit(6)];
y0_SpongeCircuit_1_4_1(y0_SpongeCircuit_1_4_1 < 0) = 0;

y0_SpongeCircuit_1_4_2 = [0, 0, SpongeCircuit_1_4_2(1,5)/k_fit(7), 0, 0, 0, 0, 0, SpongeCircuit_1_4_2(1,6)/k_fit(7)];
y0_SpongeCircuit_1_4_2(y0_SpongeCircuit_1_4_2 < 0) = 0;

y0_SpongeCircuit_1_4_3 = [0, 0, SpongeCircuit_1_4_3(1,5)/k_fit(8), 0, 0, 0, 0, 0, SpongeCircuit_1_4_3(1,6)/k_fit(8)];
y0_SpongeCircuit_1_4_3(y0_SpongeCircuit_1_4_3 < 0) = 0;

y0_SpongeCircuit_1_4_4 = [0, 0, SpongeCircuit_1_4_4(1,5)/k_fit(9), 0, 0, 0, 0, 0, SpongeCircuit_1_4_4(1,6)/k_fit(9)];
y0_SpongeCircuit_1_4_4(y0_SpongeCircuit_1_4_4 < 0) = 0;

y0_SpongeCircuit_1_3_1 = [0, 0, SpongeCircuit_1_3_1(1,5)/k_fit(10), 0, 0, 0, 0];
y0_SpongeCircuit_1_3_1(y0_SpongeCircuit_1_3_1 < 0) = 0;

y0_SpongeCircuit_1_3_2 = [0, 0, SpongeCircuit_1_3_2(1,5)/k_fit(11), 0, 0, 0, 0];
y0_SpongeCircuit_1_3_2(y0_SpongeCircuit_1_3_2 < 0) = 0;

y0_sRNACircuit = [0, 0, sRNACircuit_Array_1(1,4)/k_fit(12), 0, 0];
y0_sRNACircuit(y0_sRNACircuit<0) = 0;

y0_PVanCC_EGFP_2 = [(PVanCC_EGFP_Array2(1,4)*((PVanCC_EGFP_Array2(1,4)*PVanCC_EGFP_Array2(1,3))/F_fit(2))+ F_fit(2)*((PVanCC_EGFP_Array2(1,4)*PVanCC_EGFP_Array2(1,3))/F_fit(2)))/k_fit(31), ((PVanCC_EGFP_Array2(1,4)*PVanCC_EGFP_Array2(1,3))/F_fit(2)), PVanCC_EGFP_Array2(1,3)];
y0_PVanCC_EGFP_3 = [(PVanCC_EGFP_Array3(1,4)*((PVanCC_EGFP_Array3(1,4)*PVanCC_EGFP_Array3(1,3))/F_fit(2))+ F_fit(2)*((PVanCC_EGFP_Array3(1,4)*PVanCC_EGFP_Array3(1,3))/F_fit(2)))/k_fit(31), ((PVanCC_EGFP_Array3(1,4)*PVanCC_EGFP_Array3(1,3))/F_fit(2)), PVanCC_EGFP_Array3(1,3)];
y0_PVanCC_EGFP_4 = [(PVanCC_EGFP_Array4(1,4)*((PVanCC_EGFP_Array4(1,4)*PVanCC_EGFP_Array4(1,3))/F_fit(2))+ F_fit(2)*((PVanCC_EGFP_Array4(1,4)*PVanCC_EGFP_Array4(1,3))/F_fit(2)))/k_fit(31), ((PVanCC_EGFP_Array4(1,4)*PVanCC_EGFP_Array4(1,3))/F_fit(2)), PVanCC_EGFP_Array4(1,3)];
y0_PLuxB_EGFP_1 = [(PLuxB_EGFP_Array1(1,4)*((PLuxB_EGFP_Array1(1,4)*PLuxB_EGFP_Array1(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array1(1,4)*PLuxB_EGFP_Array1(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array1(1,4)*PLuxB_EGFP_Array1(1,3))/F_fit(2)), PLuxB_EGFP_Array1(1,3)];
y0_PLuxB_EGFP_2 = [(PLuxB_EGFP_Array2(1,4)*((PLuxB_EGFP_Array2(1,4)*PLuxB_EGFP_Array2(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array2(1,4)*PLuxB_EGFP_Array2(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array2(1,4)*PLuxB_EGFP_Array2(1,3))/F_fit(2)), PLuxB_EGFP_Array2(1,3)];
y0_PLuxB_EGFP_3 = [(PLuxB_EGFP_Array3(1,4)*((PLuxB_EGFP_Array3(1,4)*PLuxB_EGFP_Array3(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array3(1,4)*PLuxB_EGFP_Array3(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array3(1,4)*PLuxB_EGFP_Array3(1,3))/F_fit(2)), PLuxB_EGFP_Array3(1,3)];
y0_PLuxB_EGFP_4 = [(PLuxB_EGFP_Array4(1,4)*((PLuxB_EGFP_Array4(1,4)*PLuxB_EGFP_Array4(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array4(1,4)*PLuxB_EGFP_Array4(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array4(1,4)*PLuxB_EGFP_Array4(1,3))/F_fit(2)), PLuxB_EGFP_Array4(1,3)];
y0_PLuxB_EGFP_5 = [(PLuxB_EGFP_Array5(1,4)*((PLuxB_EGFP_Array5(1,4)*PLuxB_EGFP_Array5(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array5(1,4)*PLuxB_EGFP_Array5(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array5(1,4)*PLuxB_EGFP_Array5(1,3))/F_fit(2)), PLuxB_EGFP_Array5(1,3)];

y0_PTac_EGFP_2 = [(PTac_EGFP_Array2(1,4)*((PTac_EGFP_Array2(1,4)*PTac_EGFP_Array2(1,3))/F_fit(2))+ F_fit(2)*((PTac_EGFP_Array2(1,4)*PTac_EGFP_Array2(1,3))/F_fit(2)))/k_fit(31), ((PTac_EGFP_Array2(1,4)*PTac_EGFP_Array2(1,3))/F_fit(2)), PTac_EGFP_Array2(1,3)];
y0_PTac_EGFP_3 = [(PTac_EGFP_Array3(1,4)*((PTac_EGFP_Array3(1,4)*PTac_EGFP_Array3(1,3))/F_fit(2))+ F_fit(2)*((PTac_EGFP_Array3(1,4)*PTac_EGFP_Array3(1,3))/F_fit(2)))/k_fit(31), ((PTac_EGFP_Array3(1,4)*PTac_EGFP_Array3(1,3))/F_fit(2)), PTac_EGFP_Array3(1,3)];

y0_PVanCC_EGFP_2(y0_PVanCC_EGFP_2 < 0) = 0;
y0_PVanCC_EGFP_3(y0_PVanCC_EGFP_3 < 0) = 0;
y0_PVanCC_EGFP_4(y0_PVanCC_EGFP_4 < 0) = 0;
y0_PLuxB_EGFP_1(y0_PLuxB_EGFP_1 < 0) = 0;
y0_PLuxB_EGFP_2(y0_PLuxB_EGFP_2 < 0) = 0;
y0_PLuxB_EGFP_3(y0_PLuxB_EGFP_3 < 0) = 0;
y0_PLuxB_EGFP_4(y0_PLuxB_EGFP_4 < 0) = 0;

y0_PTac_EGFP_2(y0_PTac_EGFP_2 < 0) = 0;
y0_PTac_EGFP_3(y0_PTac_EGFP_3 < 0) = 0;

%% Simulating ODEs after optimisation - run simulations

% make sure the right J and mum values are used for the calcualtion
if(strcmp(mu_calc_method, 'same mum'))
    Jmums_forsim(1:2)=Jmums_est(1:2);
    Jmums_forsim(3:19)=Jmums_est(3);
elseif(strcmp(mu_calc_method, 'fixed mus'))
    Jmums_forsim = [-1, -1, Jmums_est];
else
    Jmums_forsim=Jmums_est;
end

% Sponge 1.4
[t_model_SpongeCircuit_1_4_1, y_model_SpongeCircuit_1_4_1] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_1, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(3) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_1, y0_SpongeCircuit_1_4_1, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_SpongeCircuit_1_4_2, y_model_SpongeCircuit_1_4_2] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_2, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(4) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_2, y0_SpongeCircuit_1_4_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_SpongeCircuit_1_4_3, y_model_SpongeCircuit_1_4_3] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_3, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(5) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_3, y0_SpongeCircuit_1_4_3, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_SpongeCircuit_1_4_4, y_model_SpongeCircuit_1_4_4] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_4, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(6) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_4, y0_SpongeCircuit_1_4_4, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% Sponge 1.3
[t_model_SpongeCircuit_1_3_1, y_model_SpongeCircuit_1_3_1] = ode15s( ...
    @(t, y) Sponge_1_3_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_3_1, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(7) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_3_1, y0_SpongeCircuit_1_3_1, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_SpongeCircuit_1_3_2, y_model_SpongeCircuit_1_3_2] = ode15s( ...
    @(t, y) Sponge_1_3_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_3_2, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(8) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_3_2, y0_SpongeCircuit_1_3_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% sRNA only
[t_model_sRNACircuit, y_model_sRNACircuit] = ode15s( ...
    @(t, y) sRNA_ODEs( ...
    t, y, k_fit, F_fit, sRNACircuit_Array_1, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(9) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_sRNA_circuit_1, y0_sRNACircuit, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% PVanCC characterisation
[t_model_PVanCC_EGFP2, y_model_PVanCC_EGFP2] = ode15s( ...
    @(t, y) PVanCC_EGFP_ODEs(...
        t, y, k_fit, F_fit, PVanCC_EGFP_Array2, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(10) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PVanCC_EGFP_2, y0_PVanCC_EGFP_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_PVanCC_EGFP3, y_model_PVanCC_EGFP3] = ode15s( ...
    @(t, y) PVanCC_EGFP_ODEs(...
        t, y, k_fit, F_fit, PVanCC_EGFP_Array3, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(11) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PVanCC_EGFP_3, y0_PVanCC_EGFP_3, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_PVanCC_EGFP4, y_model_PVanCC_EGFP4] = ode15s( ...
    @(t, y) PVanCC_EGFP_ODEs(...
        t, y, k_fit, F_fit, PVanCC_EGFP_Array4, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(12) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PVanCC_EGFP_4, y0_PVanCC_EGFP_4, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% PLuxB characterisation
[t_model_PLuxB_EGFP1, y_model_PLuxB_EGFP1] = ode15s( ...
    @(t, y) PLuxB_EGFP_ODEs(...
        t, y, k_fit, F_fit, PLuxB_EGFP_Array1, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(13) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PLuxB_EGFP_1, y0_PLuxB_EGFP_1, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_PLuxB_EGFP2, y_model_PLuxB_EGFP2] = ode15s( ...
    @(t, y) PLuxB_EGFP_ODEs(...
        t, y, k_fit, F_fit, PLuxB_EGFP_Array2, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(14) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PLuxB_EGFP_2, y0_PLuxB_EGFP_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_PLuxB_EGFP3, y_model_PLuxB_EGFP3] = ode15s( ...
    @(t, y) PLuxB_EGFP_ODEs(...
        t, y, k_fit, F_fit, PLuxB_EGFP_Array3, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(15) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PLuxB_EGFP_3, y0_PLuxB_EGFP_3, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_PLuxB_EGFP4, y_model_PLuxB_EGFP4] = ode15s( ...
    @(t, y) PLuxB_EGFP_ODEs(...
        t, y, k_fit, F_fit, PLuxB_EGFP_Array4, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(16) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PLuxB_EGFP_4, y0_PLuxB_EGFP_4, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_PLuxB_EGFP5, y_model_PLuxB_EGFP5] = ode15s( ...
    @(t, y) PLuxB_EGFP_ODEs(...
        t, y, k_fit, F_fit, PLuxB_EGFP_Array5, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(17) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PLuxB_EGFP_5, y0_PLuxB_EGFP_5, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% PTac characterisation
[t_model_PTac_EGFP2, y_model_PTac_EGFP2] = ode15s( ...
    @(t, y) PTac_EGFP_ODEs(...
        t, y, k_fit, F_fit, PTac_EGFP_Array2, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(18) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PTac_EGFP_2, y0_PTac_EGFP_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_PTac_EGFP3, y_model_PTac_EGFP3] = ode15s( ...
    @(t, y) PTac_EGFP_ODEs(...
        t, y, k_fit, F_fit, PTac_EGFP_Array3, ... % original arguments for the model ODE
        Jmums_forsim(1:2), Jmums_forsim(19) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_PTac_EGFP_3, y0_PTac_EGFP_3, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

%% Get growth rates predicted from simulations

mu_model_SpongeCircuit_1_4_1 = mu_from_growth_fbck(...
        y_model_SpongeCircuit_1_4_1(:,2) + y_model_SpongeCircuit_1_4_1(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        y_model_SpongeCircuit_1_4_1(:,8) + y_model_SpongeCircuit_1_4_1(:,9), ... % amount of mScarlet proteins in the cell, nascent AND mature
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(3)); % maximum possible growth rate
mu_model_SpongeCircuit_1_4_2 = mu_from_growth_fbck(...
        y_model_SpongeCircuit_1_4_2(:,2) + y_model_SpongeCircuit_1_4_2(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        y_model_SpongeCircuit_1_4_2(:,8) + y_model_SpongeCircuit_1_4_2(:,9), ... % amount of mScarlet proteins in the cell, nascent AND mature
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(4)); % maximum possible growth rate
mu_model_SpongeCircuit_1_4_3 = mu_from_growth_fbck(...
        y_model_SpongeCircuit_1_4_3(:,2) + y_model_SpongeCircuit_1_4_3(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        y_model_SpongeCircuit_1_4_3(:,8) + y_model_SpongeCircuit_1_4_3(:,9), ... % amount of mScarlet proteins in the cell, nascent AND mature
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(5)); % maximum possible growth rate
mu_model_SpongeCircuit_1_4_4 = mu_from_growth_fbck(...
        y_model_SpongeCircuit_1_4_4(:,2) + y_model_SpongeCircuit_1_4_4(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        y_model_SpongeCircuit_1_4_4(:,8) + y_model_SpongeCircuit_1_4_4(:,9), ... % amount of mScarlet proteins in the cell, nascent AND mature
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(6)); % maximum possible growth rate

mu_model_SpongeCircuit_1_3_1 = mu_from_growth_fbck(...
        y_model_SpongeCircuit_1_3_1(:,2) + y_model_SpongeCircuit_1_3_1(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none in Circuit 1.3)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(7)); % maximum possible growth rate
mu_model_SpongeCircuit_1_3_2 = mu_from_growth_fbck(...
        y_model_SpongeCircuit_1_3_2(:,2) + y_model_SpongeCircuit_1_3_2(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none in Circuit 1.3)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(8)); % maximum possible growth rate

mu_model_sRNACircuit = mu_from_growth_fbck(...
        y_model_sRNACircuit(:,2) + y_model_sRNACircuit(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for sRNA circuit)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(9)); % maximum possible growth rate

mu_model_PVanCC_EGFP2 = mu_from_growth_fbck(...
        y_model_PVanCC_EGFP2(:,2) + y_model_PVanCC_EGFP2(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(10)); % maximum possible growth rate
mu_model_PVanCC_EGFP3 = mu_from_growth_fbck(...
        y_model_PVanCC_EGFP3(:,2) + y_model_PVanCC_EGFP3(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(11)); % maximum possible growth rate
mu_model_PVanCC_EGFP4 = mu_from_growth_fbck(...
        y_model_PVanCC_EGFP4(:,2) + y_model_PVanCC_EGFP4(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(12)); % maximum possible growth rate

mu_model_PLuxB_EGFP1 = mu_from_growth_fbck(...
        y_model_PLuxB_EGFP1(:,2) + y_model_PLuxB_EGFP1(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(13)); % maximum possible growth rate
mu_model_PLuxB_EGFP2 = mu_from_growth_fbck(...
        y_model_PLuxB_EGFP2(:,2) + y_model_PLuxB_EGFP2(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(14)); % maximum possible growth rate
mu_model_PLuxB_EGFP3 = mu_from_growth_fbck(...
        y_model_PLuxB_EGFP3(:,2) + y_model_PLuxB_EGFP3(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(15)); % maximum possible growth rate
mu_model_PLuxB_EGFP4 = mu_from_growth_fbck(...
        y_model_PLuxB_EGFP4(:,2) + y_model_PLuxB_EGFP4(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(16)); % maximum possible growth rate
mu_model_PLuxB_EGFP5 = mu_from_growth_fbck(...
        y_model_PLuxB_EGFP5(:,2) + y_model_PLuxB_EGFP5(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(17)); % maximum possible growth rate


mu_model_PTac_EGFP2 = mu_from_growth_fbck(...
        y_model_PTac_EGFP2(:,2) + y_model_PTac_EGFP2(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(18)); % maximum possible growth rate
mu_model_PTac_EGFP3 = mu_from_growth_fbck(...
        y_model_PTac_EGFP3(:,2) + y_model_PTac_EGFP3(:,3), ... % amount of GFPmut3 proteins in the cell, nascent AND  mature
        0, ... % amount of mScarlet proteins in the cell (none for promoter characterisation)
        Jmums_forsim(1:2), ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        Jmums_forsim(19)); % maximum possible growth rate


%% Mean growth rates
% Common time base per circuit type
n_common = 2600;

t_common_14 = linspace(0, min([ ...
    max(SpongeCircuit_1_4_1(:,1)), max(SpongeCircuit_1_4_2(:,1)), max(SpongeCircuit_1_4_3(:,1)), max(SpongeCircuit_1_4_4(:,1)) ...
    ]), n_common);

t_common_13 = linspace(0, min([ ...
    max(SpongeCircuit_1_3_1(:,1)), max(SpongeCircuit_1_3_2(:,1)) ...
    ]), n_common);

t_common_s  = linspace(0, max(sRNACircuit_Array_1(:,1)), n_common);

% Interpolate exp growth rates onto common grids (1/h)
mu14_exp_all = nan(4, numel(t_common_14));
mu14_exp_all(1,:) = interp1(SpongeCircuit_1_4_1(:,1), SpongeCircuit_1_4_1(:,7).*3600, t_common_14, 'pchip', 'extrap');
mu14_exp_all(2,:) = interp1(SpongeCircuit_1_4_2(:,1), SpongeCircuit_1_4_2(:,7).*3600, t_common_14, 'pchip', 'extrap');
mu14_exp_all(3,:) = interp1(SpongeCircuit_1_4_3(:,1), SpongeCircuit_1_4_3(:,7).*3600, t_common_14, 'pchip', 'extrap');
mu14_exp_all(4,:) = interp1(SpongeCircuit_1_4_4(:,1), SpongeCircuit_1_4_4(:,7).*3600, t_common_14, 'pchip', 'extrap');

mu13_exp_all = nan(2, numel(t_common_13));
mu13_exp_all(1,:) = interp1(SpongeCircuit_1_3_1(:,1), SpongeCircuit_1_3_1(:,6).*3600, t_common_13, 'pchip', 'extrap');
mu13_exp_all(2,:) = interp1(SpongeCircuit_1_3_2(:,1), SpongeCircuit_1_3_2(:,6).*3600, t_common_13, 'pchip', 'extrap');

muS_exp = interp1(sRNACircuit_Array_1(:,1), sRNACircuit_Array_1(:,5).*3600, t_common_s, 'pchip', 'extrap');

% Interpolate model growth rates onto common grids (1/h)
mu14_mod_all = nan(4, numel(t_common_14));
mu14_mod_all(1,:) = interp1(t_model_SpongeCircuit_1_4_1, mu_model_SpongeCircuit_1_4_1.*3600, t_common_14, 'pchip', 'extrap');
mu14_mod_all(2,:) = interp1(t_model_SpongeCircuit_1_4_2, mu_model_SpongeCircuit_1_4_2.*3600, t_common_14, 'pchip', 'extrap');
mu14_mod_all(3,:) = interp1(t_model_SpongeCircuit_1_4_3, mu_model_SpongeCircuit_1_4_3.*3600, t_common_14, 'pchip', 'extrap');
mu14_mod_all(4,:) = interp1(t_model_SpongeCircuit_1_4_4, mu_model_SpongeCircuit_1_4_4.*3600, t_common_14, 'pchip', 'extrap');

mu13_mod_all = nan(2, numel(t_common_13));
mu13_mod_all(1,:) = interp1(t_model_SpongeCircuit_1_3_1, mu_model_SpongeCircuit_1_3_1.*3600, t_common_13, 'pchip', 'extrap');
mu13_mod_all(2,:) = interp1(t_model_SpongeCircuit_1_3_2, mu_model_SpongeCircuit_1_3_2.*3600, t_common_13, 'pchip', 'extrap');

muS_mod = interp1(t_model_sRNACircuit, mu_model_sRNACircuit.*3600, t_common_s, 'pchip', 'extrap');

% Means and SDs
mu14_exp_mean = mean(mu14_exp_all, 1, 'omitnan');
mu14_exp_std  = std( mu14_exp_all, 0, 1, 'omitnan');
mu14_mod_mean = mean(mu14_mod_all, 1, 'omitnan');

mu13_exp_mean = mean(mu13_exp_all, 1, 'omitnan');
mu13_exp_std  = std( mu13_exp_all, 0, 1, 'omitnan');
mu13_mod_mean = mean(mu13_mod_all, 1, 'omitnan');

muS_exp_mean  = muS_exp;                 % No mean (n=1)
muS_exp_std   = zeros(size(muS_exp));    % SD = 0 for n=1
muS_mod_mean  = muS_mod;

% Plot
fig_cntr = fig_cntr + 1; figure(fig_cntr); clf; hold on;

% pSS-02-006
fill([t_common_14, fliplr(t_common_14)], ...
     [mu14_exp_mean + mu14_exp_std, fliplr((mu14_exp_mean - mu14_exp_std))], ...
     'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'pSS-02-006 exp \pm SD');
plot(t_common_14, mu14_exp_mean, 'g',  'LineWidth', 2, 'DisplayName', 'pSS-02-006 exp mean');
plot(t_common_14, mu14_mod_mean, 'g--','LineWidth', 2, 'DisplayName', 'pSS-02-006 model mean');

% pSS-02-005
fill([t_common_13, fliplr(t_common_13)], ...
     [mu13_exp_mean + mu13_exp_std, fliplr((mu13_exp_mean - mu13_exp_std))], ...
     'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'pSS-02-005 exp \pm SD');
plot(t_common_13, mu13_exp_mean, 'b',  'LineWidth', 2, 'DisplayName', 'pSS-02-005 exp mean');
plot(t_common_13, mu13_mod_mean, 'b--','LineWidth', 2, 'DisplayName', 'pSS-02-005 model mean');

% pSS-02-001
plot(t_common_s, muS_exp_mean, 'm',  'LineWidth', 2, 'DisplayName', 'pSS-02-001 exp');
plot(t_common_s, muS_mod_mean, 'm--','LineWidth', 2, 'DisplayName', 'sRNA model');

xlabel('Time (seconds)');
ylabel('Growth rate (1/h)');
title('Mean growth rates by circuit type (exp mean \pm SD vs model mean)');
legend('Location', 'best');
grid off;
hold off;

%% FUNCTIONS

%% Define interpolation for inducers

function inducer = inducer_function(t_exp, I, t_model)
    inducer = interp1(t_exp, I, t_model, "previous");
end

%% Define growth rate calculation functions

% interpolate mu from experimental data
function mu_value = mu_from_exp(...
    t_exp, mu_exp, ... % times of experimental measurements and measured growth rates
    t_model...  % time for which model simulation requires finding the growth rate
    )

    mu_value = interp1(t_exp, mu_exp, t_model,"nearest");
end

% calculate mu according to LINEAR growth feedback relations
function mu_value = mu_from_growth_fbck(...
    GFPmut3_conc, mScarlet_conc, ... % amounts of GFPmut3 and mScarlet proteins in the cell
    J, ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
    mum ... % maximum possible growth rate
    )

    % calculate the burdens from each protein - zero if J set to negative value
    GFPmut3_burden = (GFPmut3_conc./J(1)) .* (J(1)>=0) .* 1e-6; % extra 1e-6 factor because protein concs are nM, J values are in mM
    mScarlet_burden = (mScarlet_conc./J(2)) .* (J(2)>=0) .*1e-6; % extra 1e-6 factor because protein concs are nM, J values are in mM
    
    % calculate the mu value according to the growth feedback model
    % mu_value_calc = mum./(1+GFPmut3_burden+mScarlet_burden);
    mu_value_calc = mum.*(1-GFPmut3_burden-mScarlet_burden)./3600; % extra 3600 division because mum is in 1/h but growth rate is in 1/s

    % clip the mu value from below at 0.5/h=1.39e-4/s to avoid unrealistic simulations
    mu_value = max(mu_value_calc, 1.39e-4);
end

%% ODE Model for Sponge_1.4 Circuit

function Sponge_1_4 = Sponge_1_4_ODEs(t, y, k_fit, F_fit, Sponge_1_4_Array,... % original arguments for the model ODE
    J, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(Sponge_1_4_Array(:,1), Sponge_1_4_Array(:,6), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = y(8) + y(9); % total (nascent and mature)  mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, J, mum);
    end
    
    % find inducer levels
    Van  = inducer_function(Sponge_1_4_Array(:,1), Sponge_1_4_Array(:,2), t);
    OC6  = inducer_function(Sponge_1_4_Array(:,1), Sponge_1_4_Array(:,3), t);
    IPTG = inducer_function(Sponge_1_4_Array(:,1), Sponge_1_4_Array(:,4), t);
    
    % find the derivative
    Sponge_1_4 = [((k_fit(19) + (k_fit(20)-k_fit(19))* ((Van^k_fit(22))/(Van^k_fit(22)+k_fit(21)^k_fit(22))))) - (mu + F_fit(1))*y(1) - ((k_fit(1) * y(1)*y(4)*k_fit(2)) / (k_fit(4) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5))));
                  k_fit(31)*y(1) - (mu + F_fit(2))*y(2);
                  F_fit(2)*y(2) - mu * y(3);
                  ((k_fit(23) + (k_fit(24)-k_fit(23))*((OC6^k_fit(26))/(OC6^k_fit(26)+k_fit(25)^k_fit(26)))))  - (mu + F_fit(3))*y(4) - ((k_fit(1) * y(1)*y(4)*k_fit(2)) / (k_fit(4) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5)))) - ((k_fit(13) * y(5)*y(4)*k_fit(2)) / (k_fit(5) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5)))) + k_fit(18)*k_fit(17)*y(6);
                  ((k_fit(27) + (k_fit(28)-k_fit(27))*((IPTG^k_fit(30))/(IPTG^k_fit(30)+k_fit(29)^k_fit(30))))) - (mu + k_fit(14))*y(5) - ((k_fit(13) * y(5)*y(4)*k_fit(2)) / (k_fit(5) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5))));
                  ((k_fit(1)   * y(1)*y(4)*k_fit(2)) / (k_fit(4) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5)))) - (k_fit(17)+mu)*y(6);
                  ((k_fit(13)  * y(5)*y(4)*k_fit(2)) / (k_fit(5) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5)))) - (k_fit(17)+mu)*y(7);
                  k_fit(15)*y(5) + k_fit(15)*y(7) - mu*y(8) - F_fit(4)*y(8);
                  F_fit(4)*y(8) - mu*y(9)];
end

%% ODE Model for Sponge_1.3 Circuit

function Sponge_1_3 = Sponge_1_3_ODEs(t, y, k_fit, F_fit, Sponge_1_3_Array,... % original arguments for the model ODE
    J, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(Sponge_1_3_Array(:,1), Sponge_1_3_Array(:,6), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = 0; % circuit 1.3 has no mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, J, mum);
    end
    
    % find inducer levels
    Van  = inducer_function(Sponge_1_3_Array(:,1), Sponge_1_3_Array(:,2), t);
    OC6  = inducer_function(Sponge_1_3_Array(:,1), Sponge_1_3_Array(:,3), t);
    IPTG = inducer_function(Sponge_1_3_Array(:,1), Sponge_1_3_Array(:,4), t);
    
    % find the derivative
    Sponge_1_3 = [((k_fit(19) + (k_fit(20)-k_fit(19))* ((Van^k_fit(22))/(Van^k_fit(22)+k_fit(21)^k_fit(22))))) - (mu + F_fit(1))*y(1) - ((k_fit(1) * y(1)*y(4)*k_fit(2)) / (k_fit(4) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5))));
                  k_fit(31)*y(1) - (mu + F_fit(2))*y(2);
                  F_fit(2)*y(2) - mu * y(3);
                  ((k_fit(23) + (k_fit(24)-k_fit(23))*((OC6^k_fit(26))/(OC6^k_fit(26)+k_fit(25)^k_fit(26)))))  - (mu + F_fit(3))*y(4) - ((k_fit(1) * y(1)*y(4)*k_fit(2)) / (k_fit(4) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5)))) - ((k_fit(13) * y(5)*y(4)*k_fit(2)) / (k_fit(5) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5)))) + k_fit(18)*k_fit(17)*y(6);
                  ((k_fit(27) + (k_fit(28)-k_fit(27))*((IPTG^k_fit(30))/(IPTG^k_fit(30)+k_fit(29)^k_fit(30))))) - (mu + k_fit(16))*y(5) - ((k_fit(13) * y(5)*y(4)*k_fit(2)) / (k_fit(5) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5))));
                  ((k_fit(1)   * y(1)*y(4)*k_fit(2)) / (k_fit(4) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5)))) - (k_fit(17)+mu)*y(6);
                  ((k_fit(13)  * y(5)*y(4)*k_fit(2)) / (k_fit(5) * (1 + (y(4)*y(1))/k_fit(4) + (y(4)*y(5))/k_fit(5)))) - (k_fit(17)+mu)*y(7)];
end

%% ODE Model for sRNA Circuit

% y = [ y(1) = GFPmut3 mRNA;   y(2) = nascent GFPmut3;    y(3) = mature GFPmut3;
%       y(4) = sRNA;        y(5) = sRNA-mRNA];

function sRNA = sRNA_ODEs(t, y, k_fit, F_fit, sRNACircuit_Array_1,... % original arguments for the model ODE
    J, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(sRNACircuit_Array(:,1), sRNACircuit_Array(:,5), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = 0; % sRNA-only circuit has no mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, J, mum);
    end

    % find inductions
    Van  = inducer_function(sRNACircuit_Array_1(:,1), sRNACircuit_Array_1(:,2), t);
    OC6  = inducer_function(sRNACircuit_Array_1(:,1), sRNACircuit_Array_1(:,3), t);

    % find the derivative
    sRNA = [((k_fit(19) + (k_fit(20)-k_fit(19))* ((Van^k_fit(22))/(Van^k_fit(22)+k_fit(21)^k_fit(22))))) - (mu + F_fit(1))*y(1) - (k_fit(1)/k_fit(4))*y(1)*y(4)*k_fit(2)/(1 + (y(4)*y(1))/k_fit(4));
            k_fit(31)*y(1) - (mu + F_fit(2))*y(2);
            F_fit(2)*y(2) - mu * y(3);
            ((k_fit(23) + (k_fit(24)-k_fit(23))*((OC6^k_fit(26))/(OC6^k_fit(26)+k_fit(25)^k_fit(26)))))  - (mu + F_fit(3))*y(4) - (k_fit(1)/k_fit(4))*y(1)*y(4)*k_fit(2)/(1 + (y(4)*y(1))/k_fit(4)) + k_fit(18)*k_fit(17)*y(5);
            (k_fit(1)/k_fit(4))*y(1)*y(4)*k_fit(2)/(1 + (y(4)*y(1))/k_fit(4)) - (mu+k_fit(17))*y(5)];
end

%% ODE Model for PVanCC-EGFP

function PVanCC_EGFP = PVanCC_EGFP_ODEs(t, y, k_fit, F_fit, PVanCC_Array,... % original arguments for the model ODE
    J, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(PVanCC_Array(:,1), PVanCC_Array(:,5), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = 0; % sRNA-only circuit has no mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, J, mum);
    end

    % find inductions
    Van = inducer_function(PVanCC_Array(:,1), PVanCC_Array(:,2), t);

    % find the derivative
    PVanCC_EGFP = [((k_fit(19) + (k_fit(20)-k_fit(19))* ((Van^k_fit(22))/(Van^k_fit(22)+k_fit(21)^k_fit(22))))) - (mu + F_fit(1))*y(1);
                   (k_fit(31)*y(1)) - mu*y(2)                                        - F_fit(2)*y(2);
                   (F_fit(2)*y(2)) - mu*y(3)];
end

%% ODE Model for PLuxB-EGFP
function PLuxB_EGFP = PLuxB_EGFP_ODEs(t, y, k_fit, F_fit, PLuxB_Array,... % original arguments for the model ODE
    J, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(PLuxB_Array(:,1), PLuxB_Array(:,5), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = 0; % sRNA-only circuit has no mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, J, mum);
    end

    % find inductions
    OC6 = inducer_function(PLuxB_Array(:,1), PLuxB_Array(:,2), t);

    % find the derivative
    PLuxB_EGFP = [((k_fit(23) + (k_fit(24)-k_fit(23))*((OC6^k_fit(26))/(OC6^k_fit(26)+k_fit(25)^k_fit(26))))) - (mu + F_fit(1))*y(1);
                   (k_fit(31)*y(1)) - mu*y(2)                                        - F_fit(2)*y(2);
                   (F_fit(2)*y(2)) - mu*y(3)];
end

%% ODE Model for PLuxB-EGFP
function PTac_EGFP = PTac_EGFP_ODEs(t, y, k_fit, F_fit, PTac_Array,... % original arguments for the model ODE
    J, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(PTac_Array(:,1), PTac_Array(:,5), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = 0; % sRNA-only circuit has no mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, J, mum);
    end

    % find inductions
    IPTG = inducer_function(PTac_Array(:,1), PTac_Array(:,2), t);

    % find the derivative
    PTac_EGFP =   [((k_fit(27) + (k_fit(28)-k_fit(27))*((IPTG^k_fit(30))/(IPTG^k_fit(30)+k_fit(29)^k_fit(30))))) - (mu + F_fit(1))*y(1);
                   (k_fit(31)*y(1)) - mu*y(2)                                        - F_fit(2)*y(2);
                   (F_fit(2)*y(2)) - mu*y(3)];
end

%% COST CALCULATIONS

%% Scott's Normal MSE (not how I would define it, but consistency matters more)
function Normal_MSE = Normal_MSE_function(exp, model)
    eps_floor = max(1, 0.01*prctile(abs(exp),95));
    denom     = max((exp + model)./2, eps_floor); 
    Normal_MSE = mean( ((exp - model)./denom).^2, 'omitnan');
end

%% Cost functions for different sRNA circuits

%  Cost function for Sponge Circuits 1.4.i
function Sponge_1_4_i_Cost = Sponge_1_4_i_Cost_Function(...
    k_fit, F_fit, SpongeCircuit_1_4_i, tspan_SpongeCircuit_1_4_i, y0_SpongeCircuit_1_4_i,... % original arguments for the model ODE
    J, mum_SpongeCircuit_1_4_i, ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    SF_SpongeCircuit_1_4_i, ... % SF factor for the particular experiment (from k_fit)
    ode_options ... % ODE integration otpions
    )

    % getting a simulated trajectory 
    [t_model, y_model] = ode15s(@(t, y) ...
        Sponge_1_4_ODEs(...
            t, y, k_fit, F_fit, SpongeCircuit_1_4_i, ... % original arguments for the model ODE
            J, mum_SpongeCircuit_1_4_i... % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
            ), ...  
        tspan_SpongeCircuit_1_4_i, y0_SpongeCircuit_1_4_i, ...  % simulation timespan and initial condition
        ode_options ... % ODE integration options
        );

    % getting GFPmut3 cost
    y_model_output = y_model(:,3).*SF_SpongeCircuit_1_4_i;
    tsin = timeseries(SpongeCircuit_1_4_i(:,5), SpongeCircuit_1_4_i(:,1));
    tsout = resample(tsin, t_model);
    y_exp_output = tsout.data;
    GFPmut3_cost = Normal_MSE_function(y_exp_output, y_model_output);

    % getting mScarlet cost
    y_model_output = y_model(:,9).*SF_SpongeCircuit_1_4_i;
    tsin = timeseries(SpongeCircuit_1_4_i(:,6), SpongeCircuit_1_4_i(:,1));
    tsout = resample(tsin, t_model);
    y_exp_output = tsout.data;
    mScarlet_cost = Normal_MSE_function(y_exp_output, y_model_output);

    % geting growth cost (from now till the end of the function)
    % get modelled growth rates
    GFPmut3_concs = y_model(:,2) + y_model(:,3);% total (nascent and mature) GFPmut3 level
    mScarlet_concs = y_model(:,8) + y_model(:,9);% total (nascent and mature) mScarlet level
    mu_model=mu_from_growth_fbck(...
        GFPmut3_concs, mScarlet_concs, ... % amounts of GFPmut3 and mScarlet proteins in the cell
        J, ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        mum_SpongeCircuit_1_4_i); % maximum possible growth rate

    % get measured growth rates (resampled at model simulation timepoints)
    tsin = timeseries(SpongeCircuit_1_4_i(:,7), SpongeCircuit_1_4_i(:,1));
    tsout = resample(tsin, t_model);
    mu_exp = tsout.data;
    
    % get a normalised MSE cost - growth rate units are 1/h
    growth_cost = Normal_MSE_function(mu_exp.*3600, mu_model.*3600);

    % gather all costs in an array
    Sponge_1_4_i_Cost = [GFPmut3_cost, mScarlet_cost, growth_cost];
end

% Cost function for Sponge Circuits 1.3.i - for observed GFPmut3, mScarlet and growth rate
function Sponge_1_3_i_Cost = Sponge_1_3_i_Cost_Function(...
    k_fit, F_fit, SpongeCircuit_1_3_i, tspan_SpongeCircuit_1_3_i, y0_SpongeCircuit_1_3_i,... % original arguments for the model ODE
    J, mum_SpongeCircuit_1_3_i, ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    SF_SpongeCircuit_1_3_i, ... % SF factor for the particular experiment (from k_fit)
    ode_options ... % ODE integration otpions
    )

    % getting a simulated trajectory 
    [t_model, y_model] = ode15s(@(t, y) ...
        Sponge_1_3_ODEs(...
            t, y, k_fit, F_fit, SpongeCircuit_1_3_i, ... % original arguments for the model ODE
            J, mum_SpongeCircuit_1_3_i ...  % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        ), ...   
        tspan_SpongeCircuit_1_3_i, y0_SpongeCircuit_1_3_i, ...  % simulation timespan and initial condition
        ode_options ... % ODE integration options
        );

    % getting GFPmut3 cost
    y_model_output = y_model(:,3).*SF_SpongeCircuit_1_3_i;
    tsin = timeseries(SpongeCircuit_1_3_i(:,5), SpongeCircuit_1_3_i(:,1));
    tsout = resample(tsin, t_model);
    y_exp_output = tsout.data;
    GFPmut3_cost = Normal_MSE_function(y_exp_output, y_model_output);

    % getting mScarlet cost
    mScarlet_cost = 0; % mScarlet not present

    % geting growth cost (from now till the end of the function)
    % get modelled growth rates
    GFPmut3_concs = y_model(:,2) + y_model(:,3);% total (nascent and mature) GFPmut3 level
    mScarlet_concs = 0; % circuit 1.3 has no mScarlet
    mu_model=mu_from_growth_fbck(...
        GFPmut3_concs, mScarlet_concs, ... % amounts of GFPmut3 and mScarlet proteins in the cell
        J, ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        mum_SpongeCircuit_1_3_i); % maximum possible growth rate

    % get measured growth rates (resampled at model simulation timepoints)
    tsin = timeseries(SpongeCircuit_1_3_i(:,6), SpongeCircuit_1_3_i(:,1));
    tsout = resample(tsin, t_model);
    mu_exp = tsout.data;
    
    % get a normalised MSE GROWTH cost - growth rate units are 1/h
    growth_cost = Normal_MSE_function(mu_exp.*3600, mu_model.*3600);

    % gather all costs in an array
    Sponge_1_3_i_Cost = [GFPmut3_cost, mScarlet_cost, growth_cost];
end

% Growth cost for sRNA-only circuits
function sRNA_i_Cost = sRNA_i_Cost_Function( ...
    k_fit, F_fit, sRNACircuit_Array_i, tspan_sRNACircuit_i, y0_sRNACircuit_i,... % original arguments for the model ODE
    J, mum_sRNACircuit_i, ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    SF_sRNACircuit_i, ... % SF factor for the particular experiment (from k_fit)
    ode_options ... % ODE integration otpions
    )
    
    % getting a simulated trajectory 
    [t_model, y_model] = ode15s(@(t, y) ...
        sRNA_ODEs(...
            t, y, k_fit, F_fit, sRNACircuit_Array_i, ... % original arguments for the model ODE
            J, mum_sRNACircuit_i... % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
            ), ...  
        tspan_sRNACircuit_i, y0_sRNACircuit_i, ...  % simulation timespan and initial condition
        ode_options ... % ODE integration options
        );
    
    % getting GFPmut3 cost
    y_model_output = y_model(:,3).*SF_sRNACircuit_i;
    tsin = timeseries(sRNACircuit_Array_i(:,4), sRNACircuit_Array_i(:,1));
    tsout = resample(tsin, t_model);
    y_exp_output = tsout.data;
    GFPmut3_cost = Normal_MSE_function(y_exp_output, y_model_output);

    % getting mScarlet cost
    mScarlet_cost = 0; % mScarlet not present

    % geting growth cost (from now till the end of the function)
    % get modelled growth rates
    GFPmut3_concs = y_model(:,2) + y_model(:,3);% total (nascent and mature) GFPmut3 level
    mScarlet_concs = 0; % sRNA circuit has no mScarlet
    mu_model=mu_from_growth_fbck(...
        GFPmut3_concs, mScarlet_concs, ... % amounts of GFPmut3 and mScarlet proteins in the cell
        J, ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        mum_sRNACircuit_i); % maximum possible growth rate

    % get measured growth rates (resampled at model simulation timepoints)
    tsin = timeseries(sRNACircuit_Array_i(:,5), sRNACircuit_Array_i(:,1));
    tsout = resample(tsin, t_model);
    mu_exp = tsout.data;

    % get a normalised MSE cost - growth rate units are 1/h
    growth_cost = Normal_MSE_function(mu_exp.*3600, mu_model.*3600);

    % gather all costs in an array
    sRNA_i_Cost = [GFPmut3_cost, mScarlet_cost, growth_cost];
end

%% Cost functions for promoter characterisation

% PVanCC
function PVanCC_EGFP_i_Cost = PVanCC_i_Cost_Function( ...
    k_fit, F_fit, PVanCC_EGFP_Array_i, tspan_PVanCC_EGFP_i, y0_PVanCC_EGFP_i,... % original arguments for the model ODE
    J, mum_PVanCC_EGFP_i, ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    SF_PVanCC_EGFP_i, ... % SF factor for the particular experiment (1 for promoter characterisation)
    ode_options ... % ODE integration otpions
    )
   
    % getting a simulated trajectory 
    [t_model, y_model] = ode15s(@(t, y) ...
        PVanCC_EGFP_ODEs(...
            t, y, k_fit, F_fit, PVanCC_EGFP_Array_i, ... % original arguments for the model ODE
            J, mum_PVanCC_EGFP_i... % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
            ), ...  
        tspan_PVanCC_EGFP_i, y0_PVanCC_EGFP_i, ...  % simulation timespan and initial condition
        ode_options ... % ODE integration options
        );
    
    % getting GFPmut3 cost
    y_model_output = y_model(:,3).*SF_PVanCC_EGFP_i;
    tsin = timeseries(PVanCC_EGFP_Array_i(:,3), PVanCC_EGFP_Array_i(:,1));
    tsout = resample(tsin, t_model);
    y_exp_output = tsout.data;
    GFPmut3_cost = Normal_MSE_function(y_exp_output, y_model_output);

    % getting mScarlet cost
    mScarlet_cost = 0; % mScarlet not present

    % geting growth cost (from now till the end of the function)
    % get modelled growth rates
    GFPmut3_concs = y_model(:,2) + y_model(:,3);% total (nascent and mature) GFPmut3 level
    mScarlet_concs = 0; % sRNA circuit has no mScarlet
    mu_model=mu_from_growth_fbck(...
        GFPmut3_concs, mScarlet_concs, ... % amounts of GFPmut3 and mScarlet proteins in the cell
        J, ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        mum_PVanCC_EGFP_i); % maximum possible growth rate

    % get measured growth rates (resampled at model simulation timepoints)
    tsin = timeseries(PVanCC_EGFP_Array_i(:,4), PVanCC_EGFP_Array_i(:,1));
    tsout = resample(tsin, t_model);
    mu_exp = tsout.data;

    % get a normalised MSE cost - growth rate units are 1/h
    growth_cost = Normal_MSE_function(mu_exp.*3600, mu_model.*3600);

    % gather all costs in an array
    PVanCC_EGFP_i_Cost = [GFPmut3_cost, mScarlet_cost, growth_cost];
end

% PLuxB
function PLuxB_EGFP_i_Cost = PLuxB_i_Cost_Function( ...
    k_fit, F_fit, PLuxB_EGFP_Array_i, tspan_PLuxB_EGFP_i, y0_PLuxB_EGFP_i,... % original arguments for the model ODE
    J, mum_PLuxB_EGFP_i, ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    SF_PLuxB_EGFP_i, ... % SF factor for the particular experiment (1 for promoter characterisation)
    ode_options ... % ODE integration otpions
    )
   
    % getting a simulated trajectory 
    [t_model, y_model] = ode15s(@(t, y) ...
        PLuxB_EGFP_ODEs(...
            t, y, k_fit, F_fit, PLuxB_EGFP_Array_i, ... % original arguments for the model ODE
            J, mum_PLuxB_EGFP_i... % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
            ), ...  
        tspan_PLuxB_EGFP_i, y0_PLuxB_EGFP_i, ...  % simulation timespan and initial condition
        ode_options ... % ODE integration options
        );
    
    % getting GFPmut3 cost
    y_model_output = y_model(:,3).*SF_PLuxB_EGFP_i;
    tsin = timeseries(PLuxB_EGFP_Array_i(:,3), PLuxB_EGFP_Array_i(:,1));
    tsout = resample(tsin, t_model);
    y_exp_output = tsout.data;
    GFPmut3_cost = Normal_MSE_function(y_exp_output, y_model_output);

    % getting mScarlet cost
    mScarlet_cost = 0; % mScarlet not present

    % geting growth cost (from now till the end of the function)
    % get modelled growth rates
    GFPmut3_concs = y_model(:,2) + y_model(:,3);% total (nascent and mature) GFPmut3 level
    mScarlet_concs = 0; % sRNA circuit has no mScarlet
    mu_model=mu_from_growth_fbck(...
        GFPmut3_concs, mScarlet_concs, ... % amounts of GFPmut3 and mScarlet proteins in the cell
        J, ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        mum_PLuxB_EGFP_i); % maximum possible growth rate

    % get measured growth rates (resampled at model simulation timepoints)
    tsin = timeseries(PLuxB_EGFP_Array_i(:,4), PLuxB_EGFP_Array_i(:,1));
    tsout = resample(tsin, t_model);
    mu_exp = tsout.data;

    % get a normalised MSE cost - growth rate units are 1/h
    growth_cost = Normal_MSE_function(mu_exp.*3600, mu_model.*3600);

    % gather all costs in an array
    PLuxB_EGFP_i_Cost = [GFPmut3_cost, mScarlet_cost, growth_cost];
end

% PTac
function PTac_EGFP_i_Cost = PTac_i_Cost_Function( ...
    k_fit, F_fit, PTac_EGFP_Array_i, tspan_PTac_EGFP_i, y0_PTac_EGFP_i,... % original arguments for the model ODE
    J, mum_PTac_EGFP_i, ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    SF_PTac_EGFP_i, ... % SF factor for the particular experiment (1 for promoter characterisation)
    ode_options ... % ODE integration otpions
    )
   
    % getting a simulated trajectory 
    [t_model, y_model] = ode15s(@(t, y) ...
        PTac_EGFP_ODEs(...
            t, y, k_fit, F_fit, PTac_EGFP_Array_i, ... % original arguments for the model ODE
            J, mum_PTac_EGFP_i... % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
            ), ...  
        tspan_PTac_EGFP_i, y0_PTac_EGFP_i, ...  % simulation timespan and initial condition
        ode_options ... % ODE integration options
        );
    
    % getting GFPmut3 cost
    y_model_output = y_model(:,3).*SF_PTac_EGFP_i;
    tsin = timeseries(PTac_EGFP_Array_i(:,3), PTac_EGFP_Array_i(:,1));
    tsout = resample(tsin, t_model);
    y_exp_output = tsout.data;
    GFPmut3_cost = Normal_MSE_function(y_exp_output, y_model_output);

    % getting mScarlet cost
    mScarlet_cost = 0; % mScarlet not present

    % geting growth cost (from now till the end of the function)
    % get modelled growth rates
    GFPmut3_concs = y_model(:,2) + y_model(:,3);% total (nascent and mature) GFPmut3 level
    mScarlet_concs = 0; % sRNA circuit has no mScarlet
    mu_model=mu_from_growth_fbck(...
        GFPmut3_concs, mScarlet_concs, ... % amounts of GFPmut3 and mScarlet proteins in the cell
        J, ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
        mum_PTac_EGFP_i); % maximum possible growth rate

    % get measured growth rates (resampled at model simulation timepoints)
    tsin = timeseries(PTac_EGFP_Array_i(:,4), PTac_EGFP_Array_i(:,1));
    tsout = resample(tsin, t_model);
    mu_exp = tsout.data;

    % get a normalised MSE cost - growth rate units are 1/h
    growth_cost = Normal_MSE_function(mu_exp.*3600, mu_model.*3600);

    % gather all costs in an array
    PTac_EGFP_i_Cost = [GFPmut3_cost, mScarlet_cost, growth_cost];
end

%% Combined cost function

function Combined_Cost = Combined_CostFunction(k_fit, F_fit, ...
    SpongeCircuit_1_3_1, tspan_SpongeCircuit_1_3_1, ...
    SpongeCircuit_1_3_2, tspan_SpongeCircuit_1_3_2, ...
    SpongeCircuit_1_4_1, tspan_SpongeCircuit_1_4_1, ...
    SpongeCircuit_1_4_2, tspan_SpongeCircuit_1_4_2, ...
    SpongeCircuit_1_4_3, tspan_SpongeCircuit_1_4_3, ...
    SpongeCircuit_1_4_4, tspan_SpongeCircuit_1_4_4, ...
    sRNACircuit_Array_1, tspan_sRNA_Circuit_1, ...
    PVanCC_EGFP_Array2, tspan_PVanCC_EGFP_2, ...
    PVanCC_EGFP_Array3, tspan_PVanCC_EGFP_3, ...
    PVanCC_EGFP_Array4, tspan_PVanCC_EGFP_4, ...
    PLuxB_EGFP_Array1, tspan_PLuxB_EGFP_1, ...
    PLuxB_EGFP_Array2, tspan_PLuxB_EGFP_2, ...
    PLuxB_EGFP_Array3, tspan_PLuxB_EGFP_3, ...
    PLuxB_EGFP_Array4, tspan_PLuxB_EGFP_4, ...
    PLuxB_EGFP_Array5, tspan_PLuxB_EGFP_5, ...
    PTac_EGFP_Array2, tspan_PTac_EGFP_2, ...
    PTac_EGFP_Array3, tspan_PTac_EGFP_3, ...  % up to now, original cost function arguments
    Jmums, ...          % parameters being fitted: GFPmut3 and mScarlet burden thresholds J; maximum growth rates for different experiments
    mu_calc_method, ...  % growth rate calclation method for which we run the fitting
    cost_weights, ... % weights of protein and growth rate fitting costs
    ode_options ... % ODE integration options
    )

    % set parameters according to the growth rate caclulation method
    if(strcmp(mu_calc_method,'same mums'))
        Js=Jmums(1:2);
        mums=Jmums(3)*ones(1,size(Jmums,2)-2);
    elseif(strcmp(mu_calc_method,'fixed mus'))
        Js=[-1, -1]; % negative burden thresholds will cause zero growth fbck
        mums=Jmums;
    else
        Js=Jmums(1:2);
        mums=Jmums(3:end);
    end
    
    % Initial conditions for RNA circuits
    y0_SpongeCircuit_1_3_1 = [0, 0, SpongeCircuit_1_3_1(1,5)/k_fit(10), 0, 0, 0, 0];
    y0_SpongeCircuit_1_3_1(y0_SpongeCircuit_1_3_1 < 0) = 0;
    
    y0_SpongeCircuit_1_3_2 = [0, 0, SpongeCircuit_1_3_2(1,5)/k_fit(11), 0, 0, 0, 0];
    y0_SpongeCircuit_1_3_2(y0_SpongeCircuit_1_3_2 < 0) = 0;
    
    y0_SpongeCircuit_1_4_1 = [0, 0, SpongeCircuit_1_4_1(1,5)/k_fit(6), 0, 0, 0, 0, 0, SpongeCircuit_1_4_1(1,6)/k_fit(6)];
    y0_SpongeCircuit_1_4_1(y0_SpongeCircuit_1_4_1 < 0) = 0;
    
    y0_SpongeCircuit_1_4_2 = [0, 0, SpongeCircuit_1_4_2(1,5)/k_fit(7), 0, 0, 0, 0, 0, SpongeCircuit_1_4_2(1,6)/k_fit(7)];
    y0_SpongeCircuit_1_4_2(y0_SpongeCircuit_1_4_2 < 0) = 0;
    
    y0_SpongeCircuit_1_4_3 = [0, 0, SpongeCircuit_1_4_3(1,5)/k_fit(8), 0, 0, 0, 0, 0, SpongeCircuit_1_4_3(1,6)/k_fit(8)];
    y0_SpongeCircuit_1_4_3(y0_SpongeCircuit_1_4_3 < 0) = 0;
    
    y0_SpongeCircuit_1_4_4 = [0, 0, SpongeCircuit_1_4_4(1,5)/k_fit(9), 0, 0, 0, 0, 0, SpongeCircuit_1_4_4(1,6)/k_fit(9)];
    y0_SpongeCircuit_1_4_4(y0_SpongeCircuit_1_4_4 < 0) = 0;
    
    y0_sRNACircuit = [0, 0, sRNACircuit_Array_1(1,4)/k_fit(12), 0, 0];
    y0_sRNACircuit(y0_sRNACircuit<0) = 0;

    % Circuit 1.4 costs
    Sponge_1_4_1_Cost = Sponge_1_4_i_Cost_Function(k_fit, F_fit, SpongeCircuit_1_4_1, tspan_SpongeCircuit_1_4_1, y0_SpongeCircuit_1_4_1, ... % original arguments for the model ODE
        Js(1:2), mums(1), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        k_fit(6), ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    Sponge_1_4_2_Cost = Sponge_1_4_i_Cost_Function(k_fit, F_fit, SpongeCircuit_1_4_2, tspan_SpongeCircuit_1_4_2, y0_SpongeCircuit_1_4_2, ... % original arguments for the model ODE
        Js(1:2), mums(2), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        k_fit(7), ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    Sponge_1_4_3_Cost = Sponge_1_4_i_Cost_Function(k_fit, F_fit, SpongeCircuit_1_4_3, tspan_SpongeCircuit_1_4_3, y0_SpongeCircuit_1_4_3, ... % original arguments for the model ODE
        Js(1:2), mums(3), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        k_fit(8), ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    Sponge_1_4_4_Cost = Sponge_1_4_i_Cost_Function(k_fit, F_fit, SpongeCircuit_1_4_4, tspan_SpongeCircuit_1_4_4, y0_SpongeCircuit_1_4_4, ... % original arguments for the model ODE
        Js(1:2), mums(4), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        k_fit(9), ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    Sponge_1_4_Costs = [Sponge_1_4_1_Cost(1), Sponge_1_4_2_Cost(1), Sponge_1_4_3_Cost(1), Sponge_1_4_4_Cost(1)];
    mean_Sponge_1_4_Cost = mean(Sponge_1_4_Costs, 'omitnan');
    Sponge_1_4_mScarlet_Costs = [Sponge_1_4_1_Cost(2), Sponge_1_4_2_Cost(2), Sponge_1_4_3_Cost(2), Sponge_1_4_4_Cost(2)];
    mean_Sponge_1_4_mScarlet_Cost = mean(Sponge_1_4_mScarlet_Costs, 'omitnan');
    Sponge_1_4_growth_Costs = [Sponge_1_4_1_Cost(3), Sponge_1_4_2_Cost(3), Sponge_1_4_3_Cost(3), Sponge_1_4_4_Cost(3)];
    mean_Sponge_1_4_growth_Cost = mean(Sponge_1_4_growth_Costs, 'omitnan');
    
    % Circuit 1.3 costs
    Sponge_1_3_1_Cost = Sponge_1_3_i_Cost_Function(k_fit, F_fit, SpongeCircuit_1_3_1, tspan_SpongeCircuit_1_3_1, y0_SpongeCircuit_1_3_1, ... % original arguments for the model ODE
        Js(1:2), mums(5), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        k_fit(10), ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    Sponge_1_3_2_Cost = Sponge_1_3_i_Cost_Function(k_fit, F_fit, SpongeCircuit_1_3_2, tspan_SpongeCircuit_1_3_2, y0_SpongeCircuit_1_3_2, ... % original arguments for the model ODE
        Js(1:2), mums(6), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        k_fit(11), ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    Sponge_1_3_Costs = [Sponge_1_3_1_Cost(1), Sponge_1_3_2_Cost(1)];
    mean_sponge_1_3_cost = mean(Sponge_1_3_Costs, 'omitnan');
    Sponge_1_3_growth_Costs = [Sponge_1_3_1_Cost(3), Sponge_1_3_2_Cost(3)];
    mean_Sponge_1_3_growth_Cost = mean(Sponge_1_3_growth_Costs, 'omitnan');
    
    % sRNA-only circuit costs
    sRNA_1_Cost = sRNA_i_Cost_Function(k_fit, F_fit, sRNACircuit_Array_1, tspan_sRNA_Circuit_1, y0_sRNACircuit, ... % original arguments for the model ODE
        Js(1:2), mums(7), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        k_fit(6), ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    sRNA_Cost = sRNA_1_Cost(1);
    sRNA_growth_Cost = sRNA_1_Cost(3);
    
    % Initial conditions for promoter characetrisation
    y0_PVanCC_2 = [(PVanCC_EGFP_Array2(1,4)*((PVanCC_EGFP_Array2(1,4)*PVanCC_EGFP_Array2(1,3))/F_fit(2))+ F_fit(2)*((PVanCC_EGFP_Array2(1,4)*PVanCC_EGFP_Array2(1,3))/F_fit(2)))/k_fit(31), ((PVanCC_EGFP_Array2(1,4)*PVanCC_EGFP_Array2(1,3))/F_fit(2)), PVanCC_EGFP_Array2(1,3)];
    y0_PVanCC_3 = [(PVanCC_EGFP_Array3(1,4)*((PVanCC_EGFP_Array3(1,4)*PVanCC_EGFP_Array3(1,3))/F_fit(2))+ F_fit(2)*((PVanCC_EGFP_Array3(1,4)*PVanCC_EGFP_Array3(1,3))/F_fit(2)))/k_fit(31), ((PVanCC_EGFP_Array3(1,4)*PVanCC_EGFP_Array3(1,3))/F_fit(2)), PVanCC_EGFP_Array3(1,3)];
    y0_PVanCC_4 = [(PVanCC_EGFP_Array4(1,4)*((PVanCC_EGFP_Array4(1,4)*PVanCC_EGFP_Array4(1,3))/F_fit(2))+ F_fit(2)*((PVanCC_EGFP_Array4(1,4)*PVanCC_EGFP_Array4(1,3))/F_fit(2)))/k_fit(31), ((PVanCC_EGFP_Array4(1,4)*PVanCC_EGFP_Array4(1,3))/F_fit(2)), PVanCC_EGFP_Array4(1,3)];
    
    y0_PLuxB_1 = [(PLuxB_EGFP_Array1(1,4)*((PLuxB_EGFP_Array1(1,4)*PLuxB_EGFP_Array1(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array1(1,4)*PLuxB_EGFP_Array1(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array1(1,4)*PLuxB_EGFP_Array1(1,3))/F_fit(2)), PLuxB_EGFP_Array1(1,3)];
    y0_PLuxB_2 = [(PLuxB_EGFP_Array2(1,4)*((PLuxB_EGFP_Array2(1,4)*PLuxB_EGFP_Array2(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array2(1,4)*PLuxB_EGFP_Array2(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array2(1,4)*PLuxB_EGFP_Array2(1,3))/F_fit(2)), PLuxB_EGFP_Array2(1,3)];
    y0_PLuxB_3 = [(PLuxB_EGFP_Array3(1,4)*((PLuxB_EGFP_Array3(1,4)*PLuxB_EGFP_Array3(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array3(1,4)*PLuxB_EGFP_Array3(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array3(1,4)*PLuxB_EGFP_Array3(1,3))/F_fit(2)), PLuxB_EGFP_Array3(1,3)];
    y0_PLuxB_4 = [(PLuxB_EGFP_Array4(1,4)*((PLuxB_EGFP_Array4(1,4)*PLuxB_EGFP_Array4(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array4(1,4)*PLuxB_EGFP_Array4(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array4(1,4)*PLuxB_EGFP_Array4(1,3))/F_fit(2)), PLuxB_EGFP_Array4(1,3)];
    y0_PLuxB_5 = [(PLuxB_EGFP_Array5(1,4)*((PLuxB_EGFP_Array5(1,4)*PLuxB_EGFP_Array5(1,3))/F_fit(2))+ F_fit(2)*((PLuxB_EGFP_Array5(1,4)*PLuxB_EGFP_Array5(1,3))/F_fit(2)))/k_fit(31), ((PLuxB_EGFP_Array5(1,4)*PLuxB_EGFP_Array5(1,3))/F_fit(2)), PLuxB_EGFP_Array5(1,3)];
    
    
    y0_PTac_2 = [(PTac_EGFP_Array2(1,4)*((PTac_EGFP_Array2(1,4)*PTac_EGFP_Array2(1,3))/F_fit(2))+ F_fit(2)*((PTac_EGFP_Array2(1,4)*PTac_EGFP_Array2(1,3))/F_fit(2)))/k_fit(31), ((PTac_EGFP_Array2(1,4)*PTac_EGFP_Array2(1,3))/F_fit(2)), PTac_EGFP_Array2(1,3)];
    y0_PTac_3 = [(PTac_EGFP_Array3(1,4)*((PTac_EGFP_Array3(1,4)*PTac_EGFP_Array3(1,3))/F_fit(2))+ F_fit(2)*((PTac_EGFP_Array3(1,4)*PTac_EGFP_Array3(1,3))/F_fit(2)))/k_fit(31), ((PTac_EGFP_Array3(1,4)*PTac_EGFP_Array3(1,3))/F_fit(2)), PTac_EGFP_Array3(1,3)];
    
    y0_PVanCC_2(y0_PVanCC_2 < 0) = 0;
    y0_PVanCC_3(y0_PVanCC_3 < 0) = 0;
    y0_PVanCC_4(y0_PVanCC_4 < 0) = 0;
    y0_PLuxB_1(y0_PLuxB_1 < 0) = 0;
    y0_PLuxB_2(y0_PLuxB_2 < 0) = 0;
    y0_PLuxB_3(y0_PLuxB_3 < 0) = 0;
    y0_PLuxB_4(y0_PLuxB_4 < 0) = 0;
    y0_PTac_2(y0_PTac_2 < 0) = 0;
    y0_PTac_3(y0_PTac_3 < 0) = 0;

    % PVanCC characterisation cost
    PVanCC_Cost_2 = PVanCC_i_Cost_Function(k_fit, F_fit, PVanCC_EGFP_Array2, tspan_PVanCC_EGFP_2, y0_PVanCC_2, ... % original arguments for the model ODE
        Js(1:2), mums(8), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PVanCC_Cost_3 = PVanCC_i_Cost_Function(k_fit, F_fit, PVanCC_EGFP_Array3, tspan_PVanCC_EGFP_3, y0_PVanCC_3, ... % original arguments for the model ODE
        Js(1:2), mums(9), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PVanCC_Cost_4 = PVanCC_i_Cost_Function(k_fit, F_fit, PVanCC_EGFP_Array4, tspan_PVanCC_EGFP_4, y0_PVanCC_4, ... % original arguments for the model ODE
        Js(1:2), mums(10), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PVanCC_Costs = [PVanCC_Cost_2(1), PVanCC_Cost_3(1), PVanCC_Cost_4(1)];
    PVanCC_Mean_Cost = mean(PVanCC_Costs, 'omitnan');
    PVanCC_growth_Costs = [PVanCC_Cost_2(3), PVanCC_Cost_3(3), PVanCC_Cost_4(3)];
    PVanCC_growth_Mean_Cost = mean(PVanCC_growth_Costs, 'omitnan');

    
    PLuxB_Cost_1 = PLuxB_i_Cost_Function(k_fit, F_fit, PLuxB_EGFP_Array1, tspan_PLuxB_EGFP_1, y0_PLuxB_1, ... % original arguments for the model ODE
        Js(1:2), mums(11), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PLuxB_Cost_2 = PLuxB_i_Cost_Function(k_fit, F_fit, PLuxB_EGFP_Array2, tspan_PLuxB_EGFP_2, y0_PLuxB_2, ... % original arguments for the model ODE
        Js(1:2), mums(12), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PLuxB_Cost_3 = PLuxB_i_Cost_Function(k_fit, F_fit, PLuxB_EGFP_Array3, tspan_PLuxB_EGFP_3, y0_PLuxB_3, ... % original arguments for the model ODE
        Js(1:2), mums(13), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PLuxB_Cost_4 = PLuxB_i_Cost_Function(k_fit, F_fit, PLuxB_EGFP_Array4, tspan_PLuxB_EGFP_4, y0_PLuxB_4, ... % original arguments for the model ODE
        Js(1:2), mums(14), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PLuxB_Cost_5 = PLuxB_i_Cost_Function(k_fit, F_fit, PLuxB_EGFP_Array5, tspan_PLuxB_EGFP_5, y0_PLuxB_5, ... % original arguments for the model ODE
        Js(1:2), mums(15), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PLuxB_Costs = [PLuxB_Cost_1(1), PLuxB_Cost_2(1), PLuxB_Cost_3(1), PLuxB_Cost_4(1), PLuxB_Cost_5(1)];
    PLuxB_Mean_Cost = mean(PLuxB_Costs, 'omitnan');
    PLuxB_growth_Costs = [PLuxB_Cost_1(3), PLuxB_Cost_2(3), PLuxB_Cost_3(3), PLuxB_Cost_4(3), PLuxB_Cost_5(3)];
    PLuxB_growth_Mean_Cost = mean(PLuxB_growth_Costs, 'omitnan');
    
    PTac_Cost_2 = PTac_i_Cost_Function(k_fit, F_fit, PTac_EGFP_Array2, tspan_PTac_EGFP_2, y0_PTac_2, ... % original arguments for the model ODE
        Js(1:2), mums(16), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PTac_Cost_3 = PTac_i_Cost_Function(k_fit, F_fit, PTac_EGFP_Array3, tspan_PTac_EGFP_3, y0_PTac_3, ... % original arguments for the model ODE
        Js(1:2), mums(17), ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
        1, ... % SF factor for the particular experiment (1 for promoter characterisation)
        ode_options); % ode integration options
    PTac_Costs = [PTac_Cost_2(1), PTac_Cost_3(1)];
    PTac_Mean_Cost = mean(PTac_Costs, 'omitnan');
    PTac_growth_Costs = [PTac_Cost_2(3), PTac_Cost_3(3)];
    PTac_growth_Mean_Cost = mean(PTac_growth_Costs, 'omitnan');
    
    % Find protein costs
    costs   = [mean_Sponge_1_4_Cost, mean_Sponge_1_4_mScarlet_Cost, mean_sponge_1_3_cost, sRNA_Cost, PVanCC_Mean_Cost, PLuxB_Mean_Cost, PTac_Mean_Cost];
    weights = [4,                    4,                             2,                    1,         3,                5,               2];
    Combined_protein_Cost = sum(weights .* costs) / sum(weights); 

    % Find growth costs
    growth_costs   = [mean_Sponge_1_4_growth_Cost,  mean_Sponge_1_3_growth_Cost, sRNA_growth_Cost, PVanCC_growth_Mean_Cost, PLuxB_growth_Mean_Cost, PTac_growth_Mean_Cost];
    growth_weights = [4,                            2,                           1,                3,                       5,                      2];
    Combined_growth_Cost = sum(growth_weights .* growth_costs) / sum(growth_weights); 

    % finmd the weighted cost
    Combined_Cost = (Combined_protein_Cost*cost_weights(1) + Combined_growth_Cost*cost_weights(2)) / sum(cost_weights);
end