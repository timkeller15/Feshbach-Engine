%% Setting parameters
N = 1e4; % number of BEC atoms
gi = 1; % initial interaction strength
gf = 0.8; % final interaction strength 
posmax = 30; % range of position grid, needs to be larger than Thomas-Fermi Radius RTF = sqrt(2*muTF)
Ngrid = 2048; % number of position grid points
[x,dx] = fftdef(posmax,Ngrid); % defines position grid

%% Find ground state
% Thomas-Fermi approximation for initial wave function guess
muTF = ((9/32)*(N*gi)^2)^(1/3);
wfi = real(sqrt((muTF - 0.5*x.^2)/gi));

dt = 1e-5; % time-step for imaginary-time evolution
Tf = 10; % total duration for imaginary-time evolution
mode = 'const'; % constant interaction term

tic; [~,d,~] = bec_interaction_ramp_1D(N,gi,gf,dt,Tf,wfi,posmax,Ngrid,mode); toc
wfi = d.wf;
clear d 

%% Interaction ramp using shortcut to adiabaticity
% For the given parameters in 1D any shortcut ramp duration Tf <~ 0.45 will
% trigger the modulational instability 
dt = 1i*1e-5; % time-step for real-time evolution
Tf = 0.4; % total duration of real-time evolution
mode = 'sta'; % shortcut interaction ramp

tic; [~,d,ani] = bec_interaction_ramp_1D(N,gi,gf,dt,Tf,wfi,posmax,Ngrid,mode); toc

% save snapshots of BEC density around t/Tf ~ 0.25 when modulational
% instability starts to show
dataout = [x abs(ani(:,d.time/Tf==0.24)).^2 abs(ani(:,d.time/Tf==0.25)).^2 abs(ani(:,d.time/Tf==0.26)).^2];

%% Write data file
header = ["x","wf1","wf2","wf3"];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/feshbach_engine_densities_instability_1D_gi1_gf08_Tf04.dat'),'Delimiter','tab');