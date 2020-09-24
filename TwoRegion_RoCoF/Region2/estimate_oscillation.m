function [a,A,T,phi] = estimate_oscillation(Oscillation,C,iterations)
% Author: Luis Badesa
%
% This function takes as input a time-series "Oscillations" where the
% underlying mathematical function is "exp(a*t)*A*sin(omega*t+phi) + C".
% "C" is given as an input to the function.
%
% It gives as output an estimate of parameters "a", "A", "T" (T=2*pi/omega) 
% and "phi". 
%
% NOTE: I am assuming that the attenuation term "a" is negative, otherwise
% this code might not work as expected.


    %% 1) Calculate the period of the oscillations "T":
    % The idea is to take the period from a maximum and a minimum of the 
    % oscillations, instead of using the zero-crossing of the oscillations
    % because it wouldn't work if there is an offset "C".
    %
    % This code also considers that the oscillation might have a phase 
    % shift "phi", and given that the oscillations are attenuated, the 
    % first maximum or minimum might be at sample t=0 but that might not be
    % a true peak in the oscillation. Therefore, that first "fake-peak" 
    % must not be considered for calculating the period of the oscillation.
    % 
    % NOTE: this methodology for calculating the period is
    % exact except if the time-series "Oscillation" does not exactly hit 
    % the oscillations peaks. This might happen if Simulink is set to run 
    % the simulation with adaptive time-step, which makes the simulation
    % faster but on the other hand uses a coarser time-step in some parts
    % of the simulation.
    % Due to this adaptive time-step, the period could be either slightly 
    % overestimated or slightly underestimated in some cases. 
    % To avoid this, manually set the Simulink to FIXED TIME-STEP, and use
    % a small value for the time-step so that the samples of the
    % oscillations are very accurate.

    [~,idx_max] = max(Oscillation.Data);
    if idx_max==1 % If the first sample is the max, discard it because it might not be a true max in the sinusoid due to the phase shift
        [~,idx_min] = min(Oscillation.Data);
        [~,idx_max] = max(Oscillation.Data(idx_min:end)); % Take the next maximum after the first minimum
        idx_max = idx_max + idx_min; % Reference the index of that second maximum to sample 0
    end 
    [~,idx_min] = min(Oscillation.Data);
    if idx_min==1 % If the first sample is the min, discard it because it might not be a true min in the sinusoid due to the phase shift
        [~,idx_min] = min(Oscillation.Data(idx_max1:end));
        idx_min = idx_min + idx_max;
    end

    T = abs(2*(Oscillation.Time(idx_max)-Oscillation.Time(idx_min)));

    
    %% 2) Estimate the phase-shift "phi"
    % The value of "phi" can be estimated from any zero-crossing of the
    % oscillations after removing the constant "C":
    % "exp(-a*t)*A*sin(omega*t+phi)+C", removing the constant becomes 
    % "exp(-a*t)*A*sin(omega*t+phi)", in any zero-crossing the sine term
    % is zero therefore "sin(omega*t+phi)=0", and as "omega" is already
    % known, you can calculate "phi" as "phi=asin(0)-omega*t_zero_crossing"
    
    Osc_without_C = Oscillation.Data - C;
    
    % Zero-crossings:
    tol=5e-3; % Use a time-step of 1e-5 in the simulation so that this tolerance works well
    zero_crossings = Osc_without_C((Osc_without_C<tol)&(Osc_without_C>-tol));
    t_zero_crossings = Oscillation.Time((Osc_without_C<tol)&(Osc_without_C>-tol));
    index_zero_crossings = find((Osc_without_C<tol)&(Osc_without_C>-tol));
    trend = sign(Osc_without_C(index_zero_crossings(2:end))-Osc_without_C(index_zero_crossings(2:end)-1)); % Discard the first "index_zero_crossings" in case it's 1, because it would make the index of the second term of this expression to be 0.

    phi_estimates = wrapTo2Pi(asin(0)-2*pi/T*t_zero_crossings);
    % If the zero-crossings correspond to a downwards trend of the
    % sinusoid, that is, if the samples before the zero-crossing are
    % positive and the samples after are negative, discard them. 
    for i=1:length(trend)
%         if trend(i)>0
%             phi_estimates(i) = phi_estimates(i) + pi;
%         end
        if trend(i)<0
            phi_estimates(i) = 0;
        end
    end
    phi_estimates = wrapTo2Pi(phi_estimates);
    phi_estimates = phi_estimates(phi_estimates>1e-3);
    phi = mean(phi_estimates);
       
    % NOTE THAT THIS ESTIMATION OF "PHI" considers a positive sin function.
    % For the area where the outage happens, the sin is actually negative,
    % so this estimated "phi" will be close to 2*pi
    
    
    %% 3) Estimate "a" and "A":

    % For estimating "a" and "A", it's better to not consider the just
    % obtained estimate for "phi" (because if that estimation has not been
    % accurate, it will introduce error in the estimation of "a" and "A").
    % Then, let's shift the oscillation to the first zero-crossing:
    Osc_shifted = Osc_without_C(index_zero_crossings(1):end);
    Osc_shifted_time = Oscillation.Time(1:length(Osc_shifted)); % Only works if you have chosen a fixed time-step in Simulink

    % The method uses a logarithmic transformation of the oscillation
    % samples so that the exponential term "a" and the log of the amplitude
    % term follow a linear relationship, which can then be estimated using
    % a linear regression. 
    
    % Remove the samples where the sin(omega*t) is almost zero (at least 
    % below the user-defined "threshold" below), because dividing by small 
    % numbers close to 0 would introduce numerical errors:
    threshold = 0.5;
    Osc_shifted_aboveThreshold = Osc_shifted(abs(sin(2*pi/T*Osc_shifted_time))>threshold);
    Osc_shifted_aboveThreshold_time = Osc_shifted_time(abs(sin(2*pi/T*Osc_shifted_time))>threshold);
    
    % Define the features matrix for the linear regression:
    X=[Osc_shifted_aboveThreshold_time ones(length(Osc_shifted_aboveThreshold_time),1)];
    
    log_forReg = log(Osc_shifted_aboveThreshold./(sin(2*pi/T*Osc_shifted_aboveThreshold_time))); % Samples for "log(A/exp(a*time_shift)) + a*t"
    opts = optimset('Display','off'); % Turn off the output message of "lsqlin"
    theta = lsqlin(X,log_forReg',[],[],[],[],[],[],[],opts); % unconstrained linear regression

    a = theta(1); % Estimate for "a" 
    A = exp(theta(2)); % Estimate for "A" but shifted, this estimate is corrected below
                       
    time_shift = Oscillation.Time(index_zero_crossings(1));
    A = A/exp(a*time_shift);
    
    A = abs(A);
    a = real(a);
    
    % Uncomment these lines if you want to visually check how good the
    % estimate for the oscillations was:   
    plot(Oscillation)
    hold on
    plot(Oscillation.Time,exp(a*Oscillation.Time)*A.*(sin(2*pi/T*Oscillation.Time))+C,'--')
    plot(Oscillation.Time,exp(a*Oscillation.Time)*A+C)
    legend('Original Oscillation','Estimated oscillation','Estimated "a" and "A"')
    
    cd('Figures Oscillations')
    print(['Fig' num2str(iterations) 'oscillations'],'-dpng')
    cd('..\')
    close all
    
    
end

