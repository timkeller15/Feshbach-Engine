function [p,d,ani] = bec_interaction_ramp_1D(N,gi,gf,dt,Tf,wfi,posmax,Ngrid,mode)
    %% Simulation parameters
    % they are chosen to ensure a reasonably fast and accurate convergence
    % under most conditions but may be adjusted of course 

    cutoff = 1e-6; % convergence criterion for the energy
    t = 0:abs(dt):Tf;
    steps = length(t);
    samples = 2000; 
    samplestep = floor(steps/samples); % interval for calculating observables
  
    %% Calculating Interation Ramp
    dim = 1; % dimension
    ai = 1; 
    af = (gf/gi)^(1/(dim + 2));

    s = t/Tf;
    a = ai + (af - ai)*(10*s.^3 - 15*s.^4 + 6*s.^5); % smoother step polynomial for scaling function
    addot = (af - ai)*(60*s - 180*s.^2 + 120*s.^3)/Tf^2; % and its second derivative
   
    switch mode
        case 'sta'
            g = gi*a.^(dim+1).*(addot + a); % shortcut ramp
        case 'tra'
            g = gi*a.^(dim+2); % time-rescaled adiabatic reference addot -> 0
        case 'lin'
            g = gi + (gf-gi)*linspace(0,1,steps); % linear ramp
        case 'const'
            g = gi*ones(1,steps); % constant interaction
        otherwise
            disp('Unknown ramp mode')
            return
    end

    %% Preparing operators and observables
    wf = wfi;

    samples = samples + 2;
    energy = zeros(1,samples); 
    Natom = zeros(1,samples); 
    time = zeros(1,samples);
    gramp = zeros(1,samples);
    ani = zeros(Ngrid,samples);

    [x,dx,p,~] = fftdef(posmax,Ngrid);
    Vtrap = 0.5*x.^2;
    Ekin = exp(-0.5*dt*p.^2);

    %% Initial values
    Ki = real(sum(0.5*conj(wfi).*ifft(p.^2.*fft(wfi))))*dx;
    Vi = real(sum(Vtrap.*abs(wfi).^2))*dx; 
    Ii = real(sum(0.5*gi*abs(wfi).^4))*dx;
    Ei = Ki + Vi + Ii;
    
    energy(1) = Ki + Vi + Ii;
    Natom(1) = sum(abs(wf).^2)*dx;
    gramp(1) = g(1);
    ani(:,1) = wfi;

    %% Fourier Split-Step method for time evolution
    count = 1;

    for i=1:steps

        % Step 1
        V = Vtrap + g(i)*abs(wf).^2;
        wf = exp(-0.5*dt*V).*wf;

        % Step 2
        wf = ifft(Ekin.*fft(wf));

        % Step 3
        V = Vtrap + g(i)*abs(wf).^2;
        wf = exp(-0.5*dt*V).*wf;

        % Normalize in case of imaginary time evolution
        if isreal(dt)
            wf = wf/sqrt(sum(conj(wf).*wf)*dx/N); 
        end
            
        %% Calculate observables at specified intervals
        if mod(i,samplestep) == 0
            count = count + 1; 
            time(count) = abs(dt)*i;

            K = real(sum(0.5*conj(wf).*ifft(p.^2.*fft(wf))))*dx;
            V = real(sum(Vtrap.*abs(wf).^2))*dx;  
            I = real(sum(0.5*g(i)*abs(wf).^4))*dx;

            energy(count) = K + V + I;
            Natom(count) = sum(abs(wf).^2)*dx;
            gramp(count) = g(i);
            ani(:,count) = wf;
            
            % check for energy convergence in case of imaginary-time
            % evolution
            if count > 2 && isreal(dt)
                diff = abs(energy(1,count) - energy(1,count-1));
                % return observables if convergence has been reached 
                if diff < cutoff 
                    time = time(1:count);
                    energy = energy(1:count);
                    Natom = Natom(1:count);
                    gramp = gramp(1:count);
                    ani = ani(1:count);
                    steps = i;
                    p = v2struct(steps,samples,dt,posmax,Ngrid,Ei,N,gi,gf,Tf,mode);
                    d = v2struct(wf,wfi,energy,time,Natom,gramp);
                    return
                end
            end 
        end

    end
    
    % final observable values
    
    K = real(sum(0.5*conj(wf).*ifft(p.^2.*fft(wf))))*dx;
    V = real(sum(Vtrap.*abs(wf).^2))*dx;  
    I = real(sum(0.5*g(end)*abs(wf).^4))*dx;
    
    energy(end) = K + V + I;
    Natom(end) = sum(abs(wf).^2)*dx;
    gramp(end) = g(end);
    time(end) = abs(dt)*i;
    ani(:,end) = wf;

    % return observables
    p = v2struct(steps,samples,dt,posmax,Ngrid,Ei,N,gi,gf,Tf,mode);
    d = v2struct(wf,wfi,energy,time,Natom,gramp);
end