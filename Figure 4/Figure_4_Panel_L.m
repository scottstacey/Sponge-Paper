%% 
clear; clc; close all

% start figure counter - at zero as will first be increased
fig_cntr = 0;

%% Importing global optimum toolbox
import globaloptim.*;

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

%% Fixed parameters from prior fits

F_fit = [0.001699,    0.002818,    0.0003610,   0.0004495];

% % Parameters set 1
k_fit = [0.000322195, 4000, ...
    1, ... % irrelevant, void k_fit(3) value)
    7152.095014,157337.8291,1.702385143,1.835658901,2.254209123,2.023631011,1.370890126,1.606015384,1.84551856,0.623690942,0.032968342,0.030786527,0.101382613,0.074767546,0.89477407,0.004803181,8.167464174,8565.435226,2.814544497,0.055457107,5.621638845,18.85233629,1.768328496,0.004710273,5.820217585,80527.43275,2.50692042,0.028545884...
    ];

% Parameter set 4
% k_fit = [0.000350397,4000, ...
%     1, ... % irrelevant, void k_fit(3) value)
%     16856.44865,17835.42047,1.66591737,1.790938682,2.298912681,1.936038641,1.327266016,1.502474243,1.692113553,0.055210966,0.067807682,0.046628111,0.146280536,0.107563206,0.857089746,0.004833399,7.619956479,7941.904818,2.936469528,0.051440964,5.906894109,18.07644568,1.867842555,0.004370751,5.817963469,77652.58211,2.667920529,0.027048179...
% ];

%% Load the GROWTH FEEDBACK model parameters

load("Data/Figure_4_Panels_KL/Qmums_est_Muf_GrPr.mat")
Qmums_fbck=Qmums_est;

%% Load the FIXED GROWTH RATE model parameters

load("Data/Figure_4_Panels_KL/Qmums_est_Muf_GrPr_fix.mat")
Qmums_fix = [-1, -1, Qmums_est];

%% Get GFPmut3 upregulations for all experimental data

% Sponge Circuits 1.3
[SpongeCircuit_1_3_1_GFPmut3_upreg_abs, SpongeCircuit_1_3_1_GFPmut3_upreg_rel] = SpongeCircuit_1_3_i_GFPmut3_upreg_exp(SpongeCircuit_1_3_1);
[SpongeCircuit_1_3_2_GFPmut3_upreg_abs, SpongeCircuit_1_3_2_GFPmut3_upreg_rel] = SpongeCircuit_1_3_i_GFPmut3_upreg_exp(SpongeCircuit_1_3_2);
% Sponge Circuits 1.4
[SpongeCircuit_1_4_1_GFPmut3_upreg_abs, SpongeCircuit_1_4_1_GFPmut3_upreg_rel] = SpongeCircuit_1_4_i_GFPmut3_upreg_exp(SpongeCircuit_1_4_1);
[SpongeCircuit_1_4_2_GFPmut3_upreg_abs, SpongeCircuit_1_4_2_GFPmut3_upreg_rel] = SpongeCircuit_1_4_i_GFPmut3_upreg_exp(SpongeCircuit_1_4_2);
[SpongeCircuit_1_4_3_GFPmut3_upreg_abs, SpongeCircuit_1_4_3_GFPmut3_upreg_rel] = SpongeCircuit_1_4_i_GFPmut3_upreg_exp(SpongeCircuit_1_4_3);
[SpongeCircuit_1_4_4_GFPmut3_upreg_abs, SpongeCircuit_1_4_4_GFPmut3_upreg_rel] = SpongeCircuit_1_4_i_GFPmut3_upreg_exp(SpongeCircuit_1_4_4);

% gather the relative upregulations
SpongeCircuit_1_3_GFPmut3_upreg_rel_gathered = [SpongeCircuit_1_3_1_GFPmut3_upreg_rel, SpongeCircuit_1_3_2_GFPmut3_upreg_rel];
SpongeCircuit_1_4_GFPmut3_upreg_rel_gathered = [SpongeCircuit_1_4_1_GFPmut3_upreg_rel, SpongeCircuit_1_4_2_GFPmut3_upreg_rel, SpongeCircuit_1_4_1_GFPmut3_upreg_rel, SpongeCircuit_1_4_2_GFPmut3_upreg_rel];

% find mean relative upregulations
SpongeCircuits_1_3_GFPmut3_upreg_rel_mean=mean(SpongeCircuit_1_3_GFPmut3_upreg_rel_gathered);
SpongeCircuits_1_4_GFPmut3_upreg_rel_mean=mean(SpongeCircuit_1_4_GFPmut3_upreg_rel_gathered);

% gather the absolute upregulations
SpongeCircuit_1_3_GFPmut3_upreg_abs_gathered = [SpongeCircuit_1_3_1_GFPmut3_upreg_abs, SpongeCircuit_1_3_2_GFPmut3_upreg_abs];
SpongeCircuit_1_4_GFPmut3_upreg_abs_gathered = [SpongeCircuit_1_4_1_GFPmut3_upreg_abs, SpongeCircuit_1_4_2_GFPmut3_upreg_abs, SpongeCircuit_1_4_1_GFPmut3_upreg_abs, SpongeCircuit_1_4_2_GFPmut3_upreg_abs];

% find mean absolute upregulations
SpongeCircuits_1_3_GFPmut3_upreg_abs_mean=mean(SpongeCircuit_1_3_GFPmut3_upreg_abs_gathered);
SpongeCircuits_1_4_GFPmut3_upreg_abs_mean=mean(SpongeCircuit_1_4_GFPmut3_upreg_abs_gathered);

%% Get induction event times for all experiments

% Sponge Circuits 1.3
[SpongeCircuit_1_3_1_which_ind, SpongeCircuit_1_2_1_all_ind_events] = SpongeCircuit_1_3_i_get_all_ind_events(SpongeCircuit_1_3_1);
[SpongeCircuit_1_3_2_which_ind, SpongeCircuit_1_3_2_all_ind_events] = SpongeCircuit_1_3_i_get_all_ind_events(SpongeCircuit_1_3_2);
% Sponge Circuits 1.4
[SpongeCircuit_1_4_1_which_ind, SpongeCircuit_1_4_1_all_ind_events] = SpongeCircuit_1_4_i_get_all_ind_events(SpongeCircuit_1_4_1);
[SpongeCircuit_1_4_2_which_ind, SpongeCircuit_1_4_2_all_ind_events] = SpongeCircuit_1_4_i_get_all_ind_events(SpongeCircuit_1_4_2);
[SpongeCircuit_1_4_3_which_ind, SpongeCircuit_1_4_3_all_ind_events] = SpongeCircuit_1_4_i_get_all_ind_events(SpongeCircuit_1_4_3);
[SpongeCircuit_1_4_4_which_ind, SpongeCircuit_1_4_4_all_ind_events] = SpongeCircuit_1_4_i_get_all_ind_events(SpongeCircuit_1_4_4);
% sRNA-only circuits
[sRNACircuit_1_which_ind, sRNACircuit_1_all_ind_events] = sRNACircuit_i_get_all_ind_events(sRNACircuit_Array_1);

%% Set options, initial conditions and time spans for model simulations

% set ODE integration options
ode_options = odeset(RelTol=1e-5,AbsTol=1e-5);

% set initial conditions for observed variables (mature FPconc) to recorded values, zeros otherwise
y0_SpongeCircuit_1_3_1 = [0, 0, SpongeCircuit_1_3_1(1,5), 0, 0, 0, 0];
y0_SpongeCircuit_1_3_1(y0_SpongeCircuit_1_3_1 < 0) = 0;
y0_SpongeCircuit_1_3_2 = [0, 0, SpongeCircuit_1_3_2(1,5), 0, 0, 0, 0];
y0_SpongeCircuit_1_3_2(y0_SpongeCircuit_1_3_2 < 0) = 0;

y0_SpongeCircuit_1_4_1 = [0, 0, SpongeCircuit_1_4_1(1,5), 0, 0, 0, 0, 0, SpongeCircuit_1_4_1(1,6)];
y0_SpongeCircuit_1_4_1(y0_SpongeCircuit_1_4_1 < 0) = 0;
y0_SpongeCircuit_1_4_2 = [0, 0, SpongeCircuit_1_4_2(1,5), 0, 0, 0, 0, 0, SpongeCircuit_1_4_2(1,6)];
y0_SpongeCircuit_1_4_2(y0_SpongeCircuit_1_4_2 < 0) = 0;
y0_SpongeCircuit_1_4_3 = [0, 0, SpongeCircuit_1_4_3(1,5), 0, 0, 0, 0, 0, SpongeCircuit_1_4_3(1,6)];
y0_SpongeCircuit_1_4_3(y0_SpongeCircuit_1_4_3 < 0) = 0;
y0_SpongeCircuit_1_4_4 = [0, 0, SpongeCircuit_1_4_4(1,5), 0, 0, 0, 0, 0, SpongeCircuit_1_4_4(1,6)];
y0_SpongeCircuit_1_4_4(y0_SpongeCircuit_1_4_4 < 0) = 0;

y0_sRNACircuit = [0, 0, sRNACircuit_Array_1(1,4), 0, 0];
y0_sRNACircuit(y0_sRNACircuit<0) = 0;

% set time spans
tspan_SpongeCircuit_1_3_1 = [0 max(SpongeCircuit_1_3_1(:,1))];
tspan_SpongeCircuit_1_3_2 = [0 max(SpongeCircuit_1_3_2(:,1))];

tspan_SpongeCircuit_1_4_1 = [0 max(SpongeCircuit_1_4_1(:,1))];
tspan_SpongeCircuit_1_4_2 = [0 max(SpongeCircuit_1_4_2(:,1))];
tspan_SpongeCircuit_1_4_3 = [0 max(SpongeCircuit_1_4_3(:,1))];
tspan_SpongeCircuit_1_4_4 = [0 max(SpongeCircuit_1_4_4(:,1))];

tspan_sRNA_circuit_1 = [0 max(sRNACircuit_Array_1(:,1))];

%% Simulate with GROWTH FEEDBACK

% Sponge 1.4
[t_model_fbck_SpongeCircuit_1_4_1, y_model_fbck_SpongeCircuit_1_4_1] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_1, ... % original arguments for the model ODE
        Qmums_fbck(1:2), Qmums_fbck(3) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_1, y0_SpongeCircuit_1_4_1, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_fbck_SpongeCircuit_1_4_2, y_model_fbck_SpongeCircuit_1_4_2] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_2, ... % original arguments for the model ODE
        Qmums_fbck(1:2), Qmums_fbck(4) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_2, y0_SpongeCircuit_1_4_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_fbck_SpongeCircuit_1_4_3, y_model_fbck_SpongeCircuit_1_4_3] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_3, ... % original arguments for the model ODE
        Qmums_fbck(1:2), Qmums_fbck(5) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_3, y0_SpongeCircuit_1_4_3, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_fbck_SpongeCircuit_1_4_4, y_model_fbck_SpongeCircuit_1_4_4] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_4, ... % original arguments for the model ODE
        Qmums_fbck(1:2), Qmums_fbck(6) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_4, y0_SpongeCircuit_1_4_4, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% Sponge 1.3
[t_model_fbck_SpongeCircuit_1_3_1, y_model_fbck_SpongeCircuit_1_3_1] = ode15s( ...
    @(t, y) Sponge_1_3_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_3_1, ... % original arguments for the model ODE
        Qmums_fbck(1:2), Qmums_fbck(7) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_3_1, y0_SpongeCircuit_1_3_1, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_fbck_SpongeCircuit_1_3_2, y_model_fbck_SpongeCircuit_1_3_2] = ode15s( ...
    @(t, y) Sponge_1_3_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_3_2, ... % original arguments for the model ODE
        Qmums_fbck(1:2), Qmums_fbck(8) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_3_2, y0_SpongeCircuit_1_3_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% sRNA only
[t_model_fbck_sRNACircuit, y_model_fbck_sRNACircuit] = ode15s( ...
    @(t, y) sRNA_ODEs( ...
    t, y, k_fit, F_fit, sRNACircuit_Array_1, ... % original arguments for the model ODE
        Qmums_fbck(1:2), Qmums_fbck(9) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_sRNA_circuit_1, y0_sRNACircuit, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

%% Simulate with FIXED GROWTH RATE

% Sponge 1.4
[t_model_fix_SpongeCircuit_1_4_1, y_model_fix_SpongeCircuit_1_4_1] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_1, ... % original arguments for the model ODE
        Qmums_fix(1:2), Qmums_fix(3) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_1, y0_SpongeCircuit_1_4_1, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_fix_SpongeCircuit_1_4_2, y_model_fix_SpongeCircuit_1_4_2] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_2, ... % original arguments for the model ODE
        Qmums_fix(1:2), Qmums_fix(4) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_2, y0_SpongeCircuit_1_4_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_fix_SpongeCircuit_1_4_3, y_model_fix_SpongeCircuit_1_4_3] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_3, ... % original arguments for the model ODE
        Qmums_fix(1:2), Qmums_fix(5) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_3, y0_SpongeCircuit_1_4_3, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_fix_SpongeCircuit_1_4_4, y_model_fix_SpongeCircuit_1_4_4] = ode15s( ...
    @(t, y) Sponge_1_4_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_4_4, ... % original arguments for the model ODE
        Qmums_fix(1:2), Qmums_fix(6) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_4_4, y0_SpongeCircuit_1_4_4, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% Sponge 1.3
[t_model_fix_SpongeCircuit_1_3_1, y_model_fix_SpongeCircuit_1_3_1] = ode15s( ...
    @(t, y) Sponge_1_3_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_3_1, ... % original arguments for the model ODE
        Qmums_fix(1:2), Qmums_fix(7) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_3_1, y0_SpongeCircuit_1_3_1, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );
[t_model_fix_SpongeCircuit_1_3_2, y_model_fix_SpongeCircuit_1_3_2] = ode15s( ...
    @(t, y) Sponge_1_3_ODEs(...
        t, y, k_fit, F_fit, SpongeCircuit_1_3_2, ... % original arguments for the model ODE
        Qmums_fix(1:2), Qmums_fix(8) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_SpongeCircuit_1_3_2, y0_SpongeCircuit_1_3_2, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

% sRNA only
[t_model_fix_sRNACircuit, y_model_fix_sRNACircuit] = ode15s( ...
    @(t, y) sRNA_ODEs( ...
    t, y, k_fit, F_fit, sRNACircuit_Array_1, ... % original arguments for the model ODE
        Qmums_fix(1:2), Qmums_fix(9) ...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    ), ...
    tspan_sRNA_circuit_1, y0_sRNACircuit, ...  % simulation timespan and initial condition
    ode_options ... % ODE integration options
    );

%% Get GFPmut3 upregulations for the GROWTH FEEDBACK model

% Sponge Circuits 1.3
[SpongeCircuit_1_3_1_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_3_1_GFPmut3_upreg_rel_fbck] = SpongeCircuit_1_3_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_3_1, ...
    t_model_fbck_SpongeCircuit_1_3_1, y_model_fbck_SpongeCircuit_1_3_1, ...
    k_fit(10));
[SpongeCircuit_1_3_2_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_3_2_GFPmut3_upreg_rel_fbck] = SpongeCircuit_1_3_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_3_2, ...
    t_model_fbck_SpongeCircuit_1_3_2, y_model_fbck_SpongeCircuit_1_3_2, ...
    k_fit(11));

% Sponge Circuits 1.4
[SpongeCircuit_1_4_1_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_4_1_GFPmut3_upreg_rel_fbck] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_1, ...
    t_model_fbck_SpongeCircuit_1_4_1, y_model_fbck_SpongeCircuit_1_4_1, ...
    k_fit(6));
[SpongeCircuit_1_4_2_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_4_2_GFPmut3_upreg_rel_fbck] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_2, ...
    t_model_fbck_SpongeCircuit_1_4_2, y_model_fbck_SpongeCircuit_1_4_2, ...
    k_fit(7));
[SpongeCircuit_1_4_3_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_4_3_GFPmut3_upreg_rel_fbck] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_3, ...
    t_model_fbck_SpongeCircuit_1_4_3, y_model_fbck_SpongeCircuit_1_4_3, ...
    k_fit(8));
[SpongeCircuit_1_4_4_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_4_4_GFPmut3_upreg_rel_fbck] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_4, ...
    t_model_fbck_SpongeCircuit_1_4_4, y_model_fbck_SpongeCircuit_1_4_4, ...
    k_fit(9));

% gather the relative upregulations
SpongeCircuit_1_3_GFPmut3_upreg_rel_gathered_fbck = [SpongeCircuit_1_3_1_GFPmut3_upreg_rel_fbck, SpongeCircuit_1_3_2_GFPmut3_upreg_rel_fbck];
SpongeCircuit_1_4_GFPmut3_upreg_rel_gathered_fbck = [SpongeCircuit_1_4_1_GFPmut3_upreg_rel_fbck, SpongeCircuit_1_4_2_GFPmut3_upreg_rel_fbck, SpongeCircuit_1_4_1_GFPmut3_upreg_rel_fbck, SpongeCircuit_1_4_2_GFPmut3_upreg_rel_fbck];

% find mean relative upregulations
SpongeCircuits_1_3_GFPmut3_upreg_rel_mean_fbck=mean(SpongeCircuit_1_3_GFPmut3_upreg_rel_gathered_fbck);
SpongeCircuits_1_4_GFPmut3_upreg_rel_mean_fbck=mean(SpongeCircuit_1_4_GFPmut3_upreg_rel_gathered_fbck);

% gather the absolute upregulations
SpongeCircuit_1_3_GFPmut3_upreg_abs_gathered_fbck = [SpongeCircuit_1_3_1_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_3_2_GFPmut3_upreg_abs_fbck];
SpongeCircuit_1_4_GFPmut3_upreg_abs_gathered_fbck = [SpongeCircuit_1_4_1_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_4_2_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_4_1_GFPmut3_upreg_abs_fbck, SpongeCircuit_1_4_2_GFPmut3_upreg_abs_fbck];

% find mean absolute upregulations
SpongeCircuits_1_3_GFPmut3_upreg_abs_mean_fbck=mean(SpongeCircuit_1_3_GFPmut3_upreg_abs_gathered_fbck);
SpongeCircuits_1_4_GFPmut3_upreg_abs_mean_fbck=mean(SpongeCircuit_1_4_GFPmut3_upreg_abs_gathered_fbck);

%% Get GFPmut3 upregulations for the NO GROWTH FEEDBACK model

% Sponge Circuits 1.3
[SpongeCircuit_1_3_1_GFPmut3_upreg_abs_fix, SpongeCircuit_1_3_1_GFPmut3_upreg_rel_fix] = SpongeCircuit_1_3_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_3_1, ...
    t_model_fix_SpongeCircuit_1_3_1, y_model_fix_SpongeCircuit_1_3_1, ...
    k_fit(10));
[SpongeCircuit_1_3_2_GFPmut3_upreg_abs_fix, SpongeCircuit_1_3_2_GFPmut3_upreg_rel_fix] = SpongeCircuit_1_3_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_3_2, ...
    t_model_fix_SpongeCircuit_1_3_2, y_model_fix_SpongeCircuit_1_3_2, ...
    k_fit(11));

% Sponge Circuits 1.4
[SpongeCircuit_1_4_1_GFPmut3_upreg_abs_fix, SpongeCircuit_1_4_1_GFPmut3_upreg_rel_fix] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_1, ...
    t_model_fix_SpongeCircuit_1_4_1, y_model_fix_SpongeCircuit_1_4_1, ...
    k_fit(6));
[SpongeCircuit_1_4_2_GFPmut3_upreg_abs_fix, SpongeCircuit_1_4_2_GFPmut3_upreg_rel_fix] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_2, ...
    t_model_fix_SpongeCircuit_1_4_2, y_model_fix_SpongeCircuit_1_4_2, ...
    k_fit(7));
[SpongeCircuit_1_4_3_GFPmut3_upreg_abs_fix, SpongeCircuit_1_4_3_GFPmut3_upreg_rel_fix] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_3, ...
    t_model_fix_SpongeCircuit_1_4_3, y_model_fix_SpongeCircuit_1_4_3, ...
    k_fit(8));
[SpongeCircuit_1_4_4_GFPmut3_upreg_abs_fix, SpongeCircuit_1_4_4_GFPmut3_upreg_rel_fix] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_4, ...
    t_model_fix_SpongeCircuit_1_4_4, y_model_fix_SpongeCircuit_1_4_4, ...
    k_fit(9));

% gather the relative upregulations
SpongeCircuit_1_3_GFPmut3_upreg_rel_gathered_fix = [SpongeCircuit_1_3_1_GFPmut3_upreg_rel_fix, SpongeCircuit_1_3_2_GFPmut3_upreg_rel_fix];
SpongeCircuit_1_4_GFPmut3_upreg_rel_gathered_fix = [SpongeCircuit_1_4_1_GFPmut3_upreg_rel_fix, SpongeCircuit_1_4_2_GFPmut3_upreg_rel_fix, SpongeCircuit_1_4_1_GFPmut3_upreg_rel_fix, SpongeCircuit_1_4_2_GFPmut3_upreg_rel_fix];

% find mean relative upregulations
SpongeCircuits_1_3_GFPmut3_upreg_rel_mean_fix=mean(SpongeCircuit_1_3_GFPmut3_upreg_rel_gathered_fix);
SpongeCircuits_1_4_GFPmut3_upreg_rel_mean_fix=mean(SpongeCircuit_1_4_GFPmut3_upreg_rel_gathered_fix);

% gather the absolute upregulations
SpongeCircuit_1_3_GFPmut3_upreg_abs_gathered_fix = [SpongeCircuit_1_3_1_GFPmut3_upreg_abs_fix, SpongeCircuit_1_3_2_GFPmut3_upreg_abs_fix];
SpongeCircuit_1_4_GFPmut3_upreg_abs_gathered_fix = [SpongeCircuit_1_4_1_GFPmut3_upreg_abs_fix, SpongeCircuit_1_4_2_GFPmut3_upreg_abs_fix, SpongeCircuit_1_4_1_GFPmut3_upreg_abs_fix, SpongeCircuit_1_4_2_GFPmut3_upreg_abs_fix];

% find mean absolute upregulations
SpongeCircuits_1_3_GFPmut3_upreg_abs_mean_fix=mean(SpongeCircuit_1_3_GFPmut3_upreg_abs_gathered_fix);
SpongeCircuits_1_4_GFPmut3_upreg_abs_mean_fix=mean(SpongeCircuit_1_4_GFPmut3_upreg_abs_gathered_fix);


%% Plot relative GFPmut3 upregulations

% plot mean relative GFPmut3 upregulation
fig_cntr = fig_cntr+1; figure(fig_cntr)
hb = bar([ ...
    [SpongeCircuits_1_3_GFPmut3_upreg_rel_mean, SpongeCircuits_1_3_GFPmut3_upreg_rel_mean_fbck, SpongeCircuits_1_3_GFPmut3_upreg_rel_mean_fix]; ...
    [SpongeCircuits_1_4_GFPmut3_upreg_rel_mean, SpongeCircuits_1_4_GFPmut3_upreg_rel_mean_fbck,  SpongeCircuits_1_4_GFPmut3_upreg_rel_mean_fix] ...
    ]);
hb(1).BarWidth = 1;
hb(2).BarWidth = 1;
hold on

% add legend
hb(1).DisplayName = 'Exp.';
hb(2).DisplayName = 'Growth fbck';
hb(3).DisplayName = 'No growth fbck';
legend(Location='northwest')
% fromat the plot
ax = gca;
ax.XTickLabel = ['Sponge Circuit 1.3';'Sponge Circuit 1.4'];
title('Relative increase in GFPmut3 from sponge induction','Models fit to protein levels')
hold off



%% FUNCTIONS - get induction events for different circuits

% Sponge Circuit 1.3
function [SpongeCircuit_1_3_i_which_ind, SpongeCircuit_1_3_i_all_ind_events] = SpongeCircuit_1_3_i_get_all_ind_events( ...
    SpongeCircuit_1_3_i ... % experimental data array
    )

    % initialise
    SpongeCircuit_1_3_i_which_ind={'GFPmut3','sRNA','sponge'};
    SpongeCircuit_1_3_all_ind_events={};

    % GFPmut3 induction by Van
    SpongeCircuit_1_3_i_all_ind_events{1} = detect_ind_events(SpongeCircuit_1_3_i(:,1), SpongeCircuit_1_3_i(:,2));
    % sRNA indcution by OC6
    SpongeCircuit_1_3_i_all_ind_events{2} = detect_ind_events(SpongeCircuit_1_3_i(:,1), SpongeCircuit_1_3_i(:,3));
    % sponge induction by IPTG
    SpongeCircuit_1_3_i_all_ind_events{3} = detect_ind_events(SpongeCircuit_1_3_i(:,1), SpongeCircuit_1_3_i(:,4));
end

% Sponge Circuit 1.4
function [SpongeCircuit_1_4_i_which_ind, SpongeCircuit_1_4_i_all_ind_events] = SpongeCircuit_1_4_i_get_all_ind_events( ...
    SpongeCircuit_1_4_i ... % experimental data array
    )

    % initialise
    SpongeCircuit_1_4_i_which_ind={'GFPmut3','sRNA','sponge'};
    SpongeCircuit_1_3_all_ind_events={};

    % GFPmut3 induction by Van
    SpongeCircuit_1_4_i_all_ind_events{1} = detect_ind_events(SpongeCircuit_1_4_i(:,1), SpongeCircuit_1_4_i(:,2));
    % sRNA indcution by OC6
    SpongeCircuit_1_4_i_all_ind_events{2} = detect_ind_events(SpongeCircuit_1_4_i(:,1), SpongeCircuit_1_4_i(:,3));
    % sponge induction by IPTG
    SpongeCircuit_1_4_i_all_ind_events{3} = detect_ind_events(SpongeCircuit_1_4_i(:,1), SpongeCircuit_1_4_i(:,4));
end

% sRNA - only circuit
function [sRNACircuit_i_which_ind, sRNACircuit_i_all_ind_events] = sRNACircuit_i_get_all_ind_events( ...
    sRNACircuit_i ... % experimental data array
    )

    % initialise
    sRNACircuit_i_which_ind={'GFPmut3','sRNA'};
    sRNACircuit_all_ind_events={};

    % GFPmut3 induction by Van
    sRNACircuit_i_all_ind_events{1} = detect_ind_events(sRNACircuit_i(:,1), sRNACircuit_i(:,2));
    % sRNA indcution by OC6
    sRNACircuit_i_all_ind_events{2} = detect_ind_events(sRNACircuit_i(:,1), sRNACircuit_i(:,3));
end


%% FUNCTION - detect induction events

function ind_events = detect_ind_events(t_exp, ind_exp)
    % initialise
    ind_events = double.empty;

    % set induction event counter
    ind_event_cntr=0;

    % read through inducer levels
    for i=2:size(t_exp,1)
        % induction event = change in inducer level from previous point
        if(ind_exp(i-1)~=ind_exp(i))
            ind_event_cntr = ind_event_cntr+1;  % update induction event counter
            ind_events(ind_event_cntr)=t_exp(i);    % record the induction event
        end
    end
end

%% FUNCTION - get datapoints in the last T seconds before a given event

function [t_T_upto_T_event, y_T_upto_T_event] = get_T_upto_T_event( ...
    t, y, ...   % time and trajectory
    T_event, ...  % time of the event
    T ...       % time window before the induction event to consider
    )
    % get indices of suitable datapoints
    indices_h_before_ind = and(T_event-T<=t,t<T_event);

    % get time points
    t_T_upto_T_event=t(indices_h_before_ind);

    % get trajectory points
    y_T_upto_T_event=y(indices_h_before_ind);
end

%% FUNCTION - get the difference in average value of variable in T seconds up to T_before and up to T_after
% (as well as the before and afetr values for reference)

function [y_change, y_before, y_after] = before_and_after( ...
    t, y, ...   % time and trajectory
    T_before, T_after, ...  % before and after times
    T ...       % time window to consider
    )
    % get mean variable value before
    [t_T_upto_T_before, y_T_upto_T_before] = get_T_upto_T_event(t, y, T_before, T);
    y_before=mean(y_T_upto_T_before);

    % get mean variable value before
    [t_T_upto_T_after, y_T_upto_T_after] = get_T_upto_T_event(t, y, T_after, T);
    y_after=mean(y_T_upto_T_after);

    % get the change
    y_change = y_after - y_before;
end

%% FUNCTIONS - get GFPmut3 upregulations by the sponge - for experimental data

% Sponge Circuit 1.3
function [SpongeCircuit_1_3_i_GFPmut3_upreg_abs, SpongeCircuit_1_3_i_GFPmut3_upreg_rel] = SpongeCircuit_1_3_i_GFPmut3_upreg_exp( ...
    SpongeCircuit_1_3_i...
    )
    % get induction events
    [SpongeCircuit_1_3_i_which_ind, SpongeCircuit_1_3_i_all_ind_events] = SpongeCircuit_1_3_i_get_all_ind_events(SpongeCircuit_1_3_i);

    % get GFPmut3 levels before any sponge induction and (the evntual level) after the second event
    [SpongeCircuit_1_4_i_GFPmut3_change, SpongeCircuit_1_4_i_GFPmut3_before, SpongeCircuit_1_3_i_GFPmut3_after] = before_and_after( ...
        SpongeCircuit_1_3_i(:,1), SpongeCircuit_1_3_i(:,5), ...   % time and trajectory
        SpongeCircuit_1_3_i_all_ind_events{3}(1), ... % before time - up to first sponge induction
        SpongeCircuit_1_3_i_all_ind_events{2}(2), ...  % after time - up to second sRNA induction
        1800 ...       % time window to consider - 30 min = 1800s
        );

    % calculate the upregulation
    SpongeCircuit_1_3_i_GFPmut3_upreg_abs = SpongeCircuit_1_4_i_GFPmut3_change;
    SpongeCircuit_1_3_i_GFPmut3_upreg_rel = SpongeCircuit_1_4_i_GFPmut3_change./SpongeCircuit_1_4_i_GFPmut3_before;
end

% Sponge Circuit 1.4
function [SpongeCircuit_1_4_i_GFPmut3_upreg_abs, SpongeCircuit_1_4_i_GFPmut3_upreg_rel] = SpongeCircuit_1_4_i_GFPmut3_upreg_exp( ...
    SpongeCircuit_1_4_i...
    )
    % get induction events
    [SpongeCircuit_1_4_i_which_ind, SpongeCircuit_1_4_i_all_ind_events] = SpongeCircuit_1_4_i_get_all_ind_events(SpongeCircuit_1_4_i);

    % get GFPmut3 levels before any sponge induction and (the evntual level) after the second event
    [SpongeCircuit_1_4_i_GFPmut3_change, SpongeCircuit_1_4_i_GFPmut3_before, SpongeCircuit_1_4_i_GFPmut3_after] = before_and_after( ...
        SpongeCircuit_1_4_i(:,1), SpongeCircuit_1_4_i(:,5), ...   % time and trajectory
        SpongeCircuit_1_4_i_all_ind_events{3}(1), ... % before time - up to first sponge induction
        SpongeCircuit_1_4_i_all_ind_events{2}(2), ...  % % after time - up to second sRNA induction
        1800 ...       % time window to consider - 30 min = 1800s
        );

    % calculate the upregulation
    SpongeCircuit_1_4_i_GFPmut3_upreg_abs = SpongeCircuit_1_4_i_GFPmut3_change;
    SpongeCircuit_1_4_i_GFPmut3_upreg_rel = SpongeCircuit_1_4_i_GFPmut3_change./SpongeCircuit_1_4_i_GFPmut3_before;
end

%% FUNCTIONS - get GFPmut3 upregulations by the sponge - for models

% Sponge Circuit 1.3
function [SpongeCircuit_1_3_i_GFPmut3_upreg_abs, SpongeCircuit_1_3_i_GFPmut3_upreg_rel] = SpongeCircuit_1_3_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_3_i, ...  % experimental data array
    t_model_SpongeCircuit_1_3_i, y_model_SpongeCircuit_1_3_i, ...   % simulated times and trajectories
    SF_SpongeCircuit_1_3_i ... % SF factor for the particular experiment (from k)
    )
    % get induction events
    [SpongeCircuit_1_3_i_which_ind, SpongeCircuit_1_3_i_all_ind_events] = SpongeCircuit_1_3_i_get_all_ind_events(SpongeCircuit_1_3_i);

    % get GFPmut3 levels before any sponge induction and (the evntual level) after the second event
    [SpongeCircuit_1_4_i_GFPmut3_change, SpongeCircuit_1_4_i_GFPmut3_before, SpongeCircuit_1_3_i_GFPmut3_after] = before_and_after( ...
        t_model_SpongeCircuit_1_3_i, y_model_SpongeCircuit_1_3_i(:,3).*SF_SpongeCircuit_1_3_i, ...   % time and trajectory
        SpongeCircuit_1_3_i_all_ind_events{3}(1), ... % before time - up to first sponge induc
        SpongeCircuit_1_3_i_all_ind_events{2}(2), ...  % before and after times
        3600 ...       % time window to consider - 30 min = 3600s
        );

    % calculate the upregulation
    SpongeCircuit_1_3_i_GFPmut3_upreg_abs = SpongeCircuit_1_4_i_GFPmut3_change;
    SpongeCircuit_1_3_i_GFPmut3_upreg_rel = SpongeCircuit_1_4_i_GFPmut3_change./SpongeCircuit_1_4_i_GFPmut3_before;
end

% Sponge Circuit 1.4
function [SpongeCircuit_1_4_i_GFPmut3_upreg_abs, SpongeCircuit_1_4_i_GFPmut3_upreg_rel] = SpongeCircuit_1_4_i_GFPmut3_upreg_model( ...
    SpongeCircuit_1_4_i, ...  % experimental data array
    t_model_SpongeCircuit_1_4_i, y_model_SpongeCircuit_1_4_i, ...   % simulated times and trajectories
    SF_SpongeCircuit_1_4_i ... % SF factor for the particular experiment (from k)
    )
    % get induction events
    [SpongeCircuit_1_4_i_which_ind, SpongeCircuit_1_4_i_all_ind_events] = SpongeCircuit_1_4_i_get_all_ind_events(SpongeCircuit_1_4_i);

    % get GFPmut3 levels before any sponge induction and (the evntual level) after the second event
    [SpongeCircuit_1_4_i_GFPmut3_change, SpongeCircuit_1_4_i_GFPmut3_before, SpongeCircuit_1_4_i_GFPmut3_after] = before_and_after( ...
        t_model_SpongeCircuit_1_4_i, y_model_SpongeCircuit_1_4_i(:,3).*SF_SpongeCircuit_1_4_i, ...   % time and trajectory
        SpongeCircuit_1_4_i_all_ind_events{3}(1), ... % before time - up to first sponge induc
        SpongeCircuit_1_4_i_all_ind_events{2}(2), ...  % before and after times
        3600 ...       % time window to consider - 30 min = 3600s
        );

    % calculate the upregulation
    SpongeCircuit_1_4_i_GFPmut3_upreg_abs = SpongeCircuit_1_4_i_GFPmut3_change;
    SpongeCircuit_1_4_i_GFPmut3_upreg_rel = SpongeCircuit_1_4_i_GFPmut3_change./SpongeCircuit_1_4_i_GFPmut3_before;
end

%% FUNCTIONS - interpolation for inducers

function inducer = inducer_function(t_exp, I, t_model)
    inducer = interp1(t_exp, I, t_model, "previous");
end

%% FUNCTIONS - growth rate calculation

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
    Q, ... % vector with burden thresholds for GFPmut3 and mScarlet, respectively
    mum ... % maximum possible growth rate
    )

    % calculate the burdens from each protein - zero if Q set to negative value
    GFPmut3_burden = (GFPmut3_conc./Q(1)) .* (Q(1)>=0) .* 1e-6; % extra 1e-6 factor because protein concs are nM, Q values are in mM
    mScarlet_burden = (mScarlet_conc./Q(2)) .* (Q(2)>=0) .*1e-6; % extra 1e-6 factor because protein concs are nM, Q values are in mM
    
    % calculate the mu value according to the growth feedback model
    % mu_value_calc = mum./(1+GFPmut3_burden+mScarlet_burden);
    mu_value_calc = mum.*(1-GFPmut3_burden-mScarlet_burden)./3600; % extra 3600 division because mum is in 1/h but growth rate is in 1/s

    % clip the mu value from below at 0.5/h=1.39e-4/s to avoid unrealistic simulations
    mu_value = max(mu_value_calc, 1.39e-4);
end

%% FUNCTIONS - ODE models

% ODE Model for Sponge_1.4 Circuit
function Sponge_1_4 = Sponge_1_4_ODEs(t, y, k_fit, F_fit, Sponge_1_4_Array,... % original arguments for the model ODE
    Q, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(Sponge_1_4_Array(:,1), Sponge_1_4_Array(:,6), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = y(8) + y(9); % total (nascent and mature)  mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, Q, mum);
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

% ODE Model for Sponge_1.3 Circuit
function Sponge_1_3 = Sponge_1_3_ODEs(t, y, k_fit, F_fit, Sponge_1_3_Array,... % original arguments for the model ODE
    Q, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(Sponge_1_3_Array(:,1), Sponge_1_3_Array(:,6), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = 0; % circuit 1.3 has no mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, Q, mum);
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

% ODE Model for sRNA Circuit
function sRNA = sRNA_ODEs(t, y, k_fit, F_fit, sRNACircuit_Array_1,... % original arguments for the model ODE
    Q, mum...   % growth feedback arguments: burden thresholds of GFPmut3 and mScarlet, maximum growth rate
    )

    % find growth rate
    if(mum<0)   % if negative max. growth rate passed, interpolate from experimental data instead
        mu   = mu_from_exp(sRNACircuit_Array(:,1), sRNACircuit_Array(:,5), t);
    else    % otherwise, estimate using growth feedback relations
        GFPmut3_conc = y(2) + y(3);    % total (nascent and mature) GFPmut3 level
        mScarlet_conc = 0; % sRNA-only circuit has no mScarlet
        mu = mu_from_growth_fbck(GFPmut3_conc, mScarlet_conc, Q, mum);
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