%% Setting parameters
Ni = 1e4; % number of BEC atoms during compression stroke
Nf = 8e3; % number of BEC atoms during expansion stroke
gi = 1; % initial interaction strength
gf = 0.8; % final interaction strength
eta_ad = 1 - (gf/gi)^(2/5); % maximum adiabatic efficiency in 3D

%% Reading groundstate energy from cycle endpoints for calculating work and heat
input = readcell(fullfile(fileparts(pwd),sprintf('/data/groundstate_N%d_g%2.1f_params.dat',Ni,gi))).';
params = cell2struct(input(2,2:end), input(1,2:end), 2);
E1 = params.Ef; 

input = readcell(fullfile(fileparts(pwd),sprintf('/data/groundstate_N%d_g%2.1f_params.dat',Nf,gf))).';
params = cell2struct(input(2,2:end), input(1,2:end), 2);
E3 = params.Ef; 

clear input params

%% Importing data
data_comp_sta = importdata(fullfile(fileparts(pwd),sprintf('/data/compression_stroke_N%d_gi%2.1f_gf%2.1f_batch_STA.dat',Ni,gi,gf)));
data_comp_tra = importdata(fullfile(fileparts(pwd),sprintf('/data/compression_stroke_N%d_gi%2.1f_gf%2.1f_batch_TRA.dat',Ni,gi,gf)));
data_exp_sta = importdata(fullfile(fileparts(pwd),sprintf('/data/expansion_stroke_N%d_gi%2.1f_gf%2.1f_batch_STA.dat',Nf,gf,gi)));
data_exp_tra = importdata(fullfile(fileparts(pwd),sprintf('/data/expansion_stroke_N%d_gi%2.1f_gf%2.1f_batch_TRA.dat',Nf,gf,gi)));

%% Computing work, heat and resulting efficiency and power
tau = 2*data_comp_sta(:,1); 
WC_STA = Ni*(data_comp_sta(:,2) - E1); 
WC_TRA = Ni*(data_comp_tra(:,2) - E1); 
WE_STA = Nf*(data_exp_sta(:,2) - E3); 
WE_TRA = Nf*(data_exp_tra(:,2) - E3); 
QP_STA = Ni*E1 - Nf*data_exp_sta(:,2);
QP_TRA = Ni*E1 - Nf*data_exp_tra(:,2);

eff_sta = -(WC_STA + WE_STA)./QP_STA;
eff_sta = eff_sta/eta_ad;
eff_tra = -(WC_TRA + WE_TRA)./QP_TRA;
eff_tra = eff_tra/eta_ad;
pow_sta = -(WC_STA + WE_STA)./tau;
pow_tra = -(WC_TRA + WE_TRA)./tau;
pow_ratio = pow_sta./pow_tra;

%% Checking for unphysical values due to modulational instability
eff_sta(eff_sta<0) = 0; eff_sta(eff_sta>1) = 0;
pow_sta(pow_sta<0) = 0;
pow_ratio(pow_ratio<0) = 0;

%% Write data file
dataout = [tau eff_sta eff_tra pow_sta pow_tra pow_ratio];
header = ["tau","eff_sta","eff_tra","pow_sta","pow_tra","pow_ratio"];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/feshbach_engine_performance.dat'),'Delimiter','tab');