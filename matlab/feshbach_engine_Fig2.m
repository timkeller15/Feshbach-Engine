%% Setting parameters
Ni = 1e4; % number of BEC atoms during compression stroke
Nf = 8e3; % number of BEC atoms during expansion stroke
gi = 1; % initial interaction strength
gf = 0.8; % final interaction strength

%% Reading groundstate energy from cycle endpoints for calculating irreversible work 
input = readcell(fullfile(fileparts(pwd),sprintf('/data/groundstate_N%d_g%2.1f_params.dat',Ni,gf))).';
params = cell2struct(input(2,2:end), input(1,2:end), 2);
E2 = params.Ef; 

input = readcell(fullfile(fileparts(pwd),sprintf('/data/groundstate_N%d_g%2.1f_params.dat',Nf,gi))).';
params = cell2struct(input(2,2:end), input(1,2:end), 2);
E4 = params.Ef; 

clear input params

%% Importing data
data_comp_sta = importdata(fullfile(fileparts(pwd),sprintf('/data/compression_stroke_N%d_gi%2.1f_gf%2.1f_batch_STA.dat',Ni,gi,gf)));
data_comp_tra = importdata(fullfile(fileparts(pwd),sprintf('/data/compression_stroke_N%d_gi%2.1f_gf%2.1f_batch_TRA.dat',Ni,gi,gf)));
data_exp_sta = importdata(fullfile(fileparts(pwd),sprintf('/data/expansion_stroke_N%d_gi%2.1f_gf%2.1f_batch_STA.dat',Nf,gf,gi)));
data_exp_tra = importdata(fullfile(fileparts(pwd),sprintf('/data/expansion_stroke_N%d_gi%2.1f_gf%2.1f_batch_TRA.dat',Nf,gf,gi)));

%% Computing quantities
Tf = data_comp_sta(:,1);

% compression stroke
wirrc_sta = Ni*(data_comp_sta(:,2) - E2);
wirrc_tra = Ni*(data_comp_tra(:,2) - E2);
fidc_sta = data_comp_sta(:,3);
fidc_tra = data_comp_tra(:,3);

% expansion stroke
wirre_sta = Nf*(data_exp_sta(:,2) - E4);
wirre_tra = Nf*(data_exp_tra(:,2) - E4);
fide_sta = data_exp_sta(:,3);
fide_tra = data_exp_tra(:,3);

%% Fixing coincidental dips in Wirr
% Due to the underlying harmonic trapping potential both the irreversible
% work and the fidelity can be subject to coincidental dips and peaks
% respectively.
% In order to avoid unphysical negative values for the irreversible work
% and since we want to plot it on a logarithmic scale, we adjust the
% irreversible work by its minimum value plus some cosmetic offset for
% plotting, without altering its qualitative behavior.

if min(wirrc_sta) < 0
    wirrc_sta = wirrc_sta + abs(min(wirrc_sta)) + 1e-6;
end

if min(wirrc_tra) < 0
    wirrc_tra = wirrc_tra + abs(min(wirrc_tra)) + 1e-3;
end

if min(wirre_sta) < 0
    wirre_sta = wirre_sta + abs(min(wirre_sta)) + 1e-6;
end

if min(wirre_tra) < 0
    wirre_tra = wirre_tra + abs(min(wirre_tra)) + 1e-3;
end

%% Write data file
dataout = [Tf wirrc_sta wirrc_tra fidc_sta fidc_tra wirre_sta wirre_tra fide_sta fide_tra];
header = ["Tf","wirrc_sta","wirrc_tra","fidc_sta","fidc_tra","wirre_sta","wirre_tra","fide_sta","fide_tra"];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/feshbach_engine_wirr_fid_3D_gi1_gf08.dat'),'Delimiter','tab');
