This folder contains a modified barebone version of the GPU accelerated Gross-Pitaevskii equation solver GPUE (https://gpue-group.github.io/) as well as the Matlab scripts and resulting data used to create the figures in the publication

"Feshbach engine in the Thomas-Fermi regime"
Tim Keller, Thomás Fogarty, Jing Li, and Thomas Busch
Physical Review Research 2, 033335 (2020).

GPUE_Ramp can be built with CMake (requires version >= 3.8) from the CMakeLists.txt file via

	cmake .
	make 

or alternatively directly via make from the included Makefile after specifying the cuda path and architecture of your GPU device in the first two lines.
See also the GPUE documentation at https://gpue-group.github.io/. 

WARNING: At the time of writing, using CUDA v11 onwards leads to kernel errors during execution. Original results were obtained with CMAKE v3.11.4 and CUDA v10.0.130 on a Tesla P100-SXM2-16GB GPU and re-checked with CUDA v10.2. 

The 3D simulation data used in Figures 2 and 4 and Figure 7 can be obtained after preparing the necessary files containing the ramp durations e.g. in Matlab via

	Tf_arr = [0.01:.01:1.5 1.6:.1:3 3.01:.01:4 4.1:.1:10].';
	writematrix(Tf_arr,fullfile(fileparts(pwd),'/data/feshbach_engine_ramp_durations_Fig2.txt'));

for the data used in Figure 2 and 4 and via

	Tf_arr = logspace(-2,0,100).'; 
	writematrix(Tf_arr,fullfile(fileparts(pwd),'/data/feshbach_engine_ramp_durations_Fig7.txt'));

for the data used in Figure 7. 

The Thomas-Fermi approximations for the input wave functions for finding the ground states need to be prepared e.g. in Matlab via

	cd(fullfile(pwd,'/matlab/'))
	run('feshbach_engine_gpue_prep.m')

Finally, the GPUE_RAMP simulations can be run via

	bash feshbach_engine_Fig2.sh

and	

	bash feshbach_engine_Fig7.sh
	
respectively. 

The available ramp modes are 

	 --rampmode=0 = time-rescaled adiabatic reference ramp (TRA)
	 --rampmode=1 = shortcut to adiabaticity using a smoother step polynomial (STA)
	 --rampmode=2 = constant interaction 
	 
The general usage structure is for example

	- imaginary-time evolution with constant interaction for finding the ground state
		./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --Tf=1.4 --imaginarytime --atoms=10000 --gi=1.0 --rampmode=2 --input=groundstate --writewfc
		
	- real-time evolution for a single TRA ramp
		./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --atoms=10000 --gi=1.0 --gf=0.8 --Tf=1.0 --rampmode=0 --input=groundstate --target=groundstate
		
	- real-time evolution for multiple STA ramp durations Tf specified in the file ramp_durations.txt
		./gpue_ramp --Ngrid=256 --posmax=10 --dt=1e-4 --atoms=10000 --gi=1.0 --gf=0.8 --rampmode=1 --input=groundstate --target=groundstate --batchmode
		
For finding the ground states an imaginary-time evolution of around Tf = 1.4 is usually already enough to achieve convergence to a precision below 1e-6, using the Thomas-Fermi wave function as the initial guess.

Simulations for the 1D BEC as well as data analysis for the 3D GPUE simulations were performed with MATLAB R2019a and can be executed e.g. via

	cd(fullfile(pwd,'/matlab/'))
	run('feshbach_engine_Fig1.m')

External dependencies are the additional files 'fftdef.m' and 'v2struct.m' which are included in the matlab folder. 

The figures are created via tikz in the .tex file for the publication itself. 
Tikz settings need to be included in the preamble via

	\input{figures/tikz_settings.tex}

The figures are then created at the desired locations via

	\begin{figure}
	\centering
	\input{figures/feshbach_engine_Fig1.tikz}
	\caption{Caption}
	\label{fig:Fig1}
	\end{figure}
	
For questions please contact tim.keller@oist.jp 	