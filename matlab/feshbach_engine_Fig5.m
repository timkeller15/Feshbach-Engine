%% Setting parameters
dim = 1; % dimension
gi = 1; % initial interaction strength
gf = 0.8; % final interaction strength 
ai = 1;
af = (gf/gi)^(1/(dim + 2));
s = linspace(0,1,100);
dataout = [s.'];

%% Constructing Interaction Ramps
for Tf = [0.5 0.4 1e8]
    a = ai + (af - ai)*(10*s.^3 - 15*s.^4 + 6*s.^5); % smoother step polynomial for scaling function
    addot = (af - ai)*(60*s - 180*s.^2 + 120*s.^3)/Tf^2; % and its second derivative
    g = gi*a.^(dim+1).*(addot + a); % shortcut ramp 
    dataout = [dataout g.'];
end

%% Write data file
header = ["t","Tf1","Tf2","Tf3"];
dataout = [header; dataout];
writematrix(dataout,fullfile(fileparts(pwd),'/data/feshbach_engine_interaction_ramps_1D_gi1_gf08.dat'),'Delimiter','tab');