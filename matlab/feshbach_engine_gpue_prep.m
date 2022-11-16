%% Setting parameters
posmax = 10; % range of position grid, needs to be larger than Thomas-Fermi Radius RTF = sqrt(2*muTF)
Ngrid = 256; % number of position grid points

[x,dx] = fftdef(posmax,Ngrid); % defines position grid 
[xm,ym,zm] = meshgrid(x,x,x);

%% Write initial states for GPUE imaginary time evolution
% Fig. 2, 4 and 7 compression strokes
N = 1e4;
for g = [0.1:.1:1]
    muTF = (15*N*g/(sqrt(2)*16*pi))^(2/5); 
    % Thomas-Fermi approximation as an initial guess for finding the ground
    % state. Normalized to one for usage with GPUE. 
    wfTF = real(sqrt((muTF - 0.5*(xm.^2 + ym.^2 + zm.^2))/(N*g)));
    wfTF = reshape(wfTF,Ngrid*Ngrid*Ngrid,1); 
    fname = sprintf('/data/groundstate_N%d_g%2.1f',N,g);
    writematrix(wfTF,fullfile(fileparts(pwd),fname + "_real.dat"));
    writematrix(zeros(Ngrid*Ngrid*Ngrid,1),fullfile(fileparts(pwd),fname + "_imag.dat"));
end

%% Write initial states for GPUE imaginary time evolution
% Fig. 2 and 4 expansion strokes
N = 8e3;
for g = [0.8 1]
    muTF = (15*N*g/(sqrt(2)*16*pi))^(2/5); 
    % Thomas-Fermi approximation as an initial guess for finding the ground
    % state. Normalized to one for usage with GPUE. 
    wfTF = real(sqrt((muTF - 0.5*(xm.^2 + ym.^2 + zm.^2))/(N*g)));
    wfTF = reshape(wfTF,Ngrid*Ngrid*Ngrid,1); 
    fname = sprintf('/data/groundstate_N%d_g%2.1f',N,g);
    writematrix(wfTF,fullfile(fileparts(pwd),fname + "_real.dat"));
    writematrix(zeros(Ngrid*Ngrid*Ngrid,1),fullfile(fileparts(pwd),fname + "_imag.dat"));
end