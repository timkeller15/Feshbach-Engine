%% Setting parameters
N = 10000; % number of BEC atoms
tres = 300; % resolution in Tf for analytic curves
gres = 300; % resolution in gf for analytic curves
ai = 1; 
dt = 1e-5; % integration time step for analytic determination of stability 

% Combinations of dimension and initial interaction strength
tuples = [1 1; 1 2; 3 1];
dataout = [];

%% Calculating stability analytically according to Eq. (27) in the paper
% Setting the threshold for instability as \Delta = 1e14 yielded the best
% fit to our numerical data, but even changing it by several orders of 
% magnitude doesn't affect the resulting curves much 
for i = 1:length(tuples)
    dim = tuples(i,1); % dimension
    gi = tuples(i,2); % initial interaction strength
    gfs = gi*linspace(0.01,1,gres); % final interaction strengths

    switch dim
        case 1
            mui = (9*(N*gi)^2/32)^(1/3);
            thresh = 1e14; 
            Tfs = linspace(0.01,2,tres);
        case 2
            mui = sqrt(N*gi/pi);
            thresh = 1e14;
            Tfs = linspace(0.01,2,tres);
        case 3
            mui = (15*N*gi/(sqrt(2)*16*pi))^(2/5); 
            thresh = 1e14;
            Tfs = linspace(0.01,0.8,tres);
    end


    border = zeros(1,length(gfs));
    for j = 1:length(gfs)
        for l = 1:length(Tfs)
            gf = gfs(j);
            af = (gf/gi)^(1/(dim+2)); 
            Tf = Tfs(l);
            t = 0:dt:Tf;
            s = t/Tf;

            a = ai + (af - ai)*(10*s.^3 - 15*s.^4 + 6*s.^5); % smoother step polynomial for scaling function
            addot = (af - ai)*(60*s - 180*s.^2 + 120*s.^3)/Tf^2; % and its second derivative 

            % stability criterion from Eq. (27)
            f = mui*(a.*addot + a.^2);
            f(f>0) = 0; f = abs(f); 

            if exp(trapz(t,f)) < thresh 
                border(j) = Tf;
                break
            end
        end
    end
    
    dataout = [dataout gfs.' border.'];
end

%% Write data file
header = ["g1","T1","g2","T2","g3","T3"];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/feshbach_engine_stability_analytic.dat'),'Delimiter','tab');

%% Determining stability numerically in 1D 
N = 1e4; % number of BEC atoms
posmax = 60; % range of position grid, needs to be larger than Thomas-Fermi Radius RTF = sqrt(2*muTF)
Ngrid = 4096; % number of position grid points
[x,dx] = fftdef(posmax,Ngrid); % defines position grid

dim = 1;
dt = 1i*1e-5; % time-step for real-time evolution

dataout = [];

for gi = [1 2]
    
    cyclelen = 1.6:-.01:0.01; % range of ramp durations Tf
    interactions = gi*(0.1:.02:0.99); % range of final interactions gf
    border = NaN(1,length(interactions));
    
    % Thomas-Fermi approximation for initial wave function guess
    muTF = ((9/32)*(N*gi)^2)^(1/3);
    wfi = real(sqrt((muTF - 0.5*x.^2)/gi));
    tic; [~,d,~] = bec_interaction_ramp_1D(N,gi,gi,1e-5,10,wfi,posmax,Ngrid,'const'); toc
    wfi = d.wf;
    clear p d 

    for i = 1:length(interactions)
        gf = interactions(i);
        for j = 1:length(cyclelen)
            Tf = cyclelen(j);

            [~,d,~] = bec_interaction_ramp_1D(N,gi,gf,dt,Tf,wfi,posmax,Ngrid,'sta');

            % Stability treshold is determined from the ramp duration where
            % the irreversible work exceeds 1e4, which proved to be a
            % criterion yielding good agreement with the actual onset of
            % the modulational instability over the range of interactions
            % considered here.
            % Furthermore, in 1D the Thomas-Fermi approximation for the
            % energy of the target state is sufficient for calculating the
            % irreversible work. 

            muTF = ((9/32)*(N*gf)^2)^(1/3);
            ETF = 3*N*muTF/5;
            if abs(d.energy(end) - ETF) > 1e4
               border(i) = cyclelen(j-1); 
               cyclelen = cyclelen(j:end);
               break
            end
        end
    end
    dataout = [dataout interactions.' border.'];
end

%% Write data file
header = ["gf1","Tf1","gf2","Tf2"];
dataout = [header; dataout];
dataout = fillmissing(dataout,'constant',"NaN");
writematrix(dataout,fullfile(fileparts(pwd),'/data/feshbach_engine_stability_numeric1D.dat'),'Delimiter','tab');

%% Determining stability numerically in 3D 
N = 1e4; % number of BEC atoms
gi = 1; % initial interaction strength
interactions = 0.1:.1:0.9; % range of final interaction strength

% In contrast to the 1D case we found no reliable criterion for stability
% for the 3D simulations. Therefore, minimum ramp durations were extracted
% manually from the inflection points appearing when plotting the
% irreversible work vs ramp duration on a semilog scale. 

figure; hold all; 
for gf = [0.1:.1:.9]
   data = importdata(fullfile(fileparts(pwd),sprintf('/data/stability_N%d_gi%2.1f_gf%2.1f_batch_STA.dat',N,gi,gf))); 
   input = readcell(fullfile(fileparts(pwd),sprintf('/data/groundstate_N%d_g%2.1f_params.dat',N,gf))).';
   params = cell2struct(input(2,2:end), input(1,2:end), 2); 

   plot(data(:,1),N*(data(:,2) - params.Ef),'o-','DisplayName',sprintf('gf=%2.1f',gf)); 
   set(gca, 'YScale', 'log'); xlabel('Tf'); ylabel('$W_\mathrm{irr}$'); legend; grid on; 
   clear input params data gf 
end

% manually extracted values for T_f^\mathrm{min}
border = [0.37649 0.29836 0.23645 0.17886 0.14175 0.11233 0.077426 0.048626 0.020092];

%% Write data file
dataout = [interactions.' border.'];
header = ["gf","Tf"];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/feshbach_engine_stability_numeric3D.dat'),'Delimiter','tab');
