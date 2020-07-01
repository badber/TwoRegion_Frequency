clearvars
close all
clc

%%
% Comments:
%
% The name of the parameters in the Simulink model has to be the same as 
% the name of the variables that are being modified in the loops below in
% the code.
%
% Use the block "To workspace" in Simulink for saving the output of a
% signal. The value will be saved in a variable with the same name as the
% one given to block "To workspace". For example, in this case some of the 
% saved variables from Simulink are "Delta_f1" and "Delta_f2"

%% System condition and parameters
D = [0.5e-2 0.5e-2]; % Units 1/Hz
Demand = 25e3; % Units: MW
Td = 10; % Units: s. Delivery time of PFR.

X = 50; % Units: ohms. Line reactance
V = 345; % Units: kV (voltage of the transmission line) 
Line_term = 2*pi*V^2/X;

%%
% Some other system parameters, if you change these you will have to change
% the loops later in the code.
RoCoF_max = 0.5; % Relaxed RoCoF
Delta_fmax = 0.8;
Delta_fDB = 0;

%% Analyse the inter-area oscillations for many different system conditions:
% NOTE: this code only considers system conditions that meet the
% dynamic-frequency requirements for the single-area equivalent system, 
% because inter-area oscillations will only make things worse so it is not
% relevant to consider conditions that already violate the single-area
% requirements.
%
% For more details of single-area conditions for frequency security, refer
% to paper "Stochastic Scheduling With Inertia-Dependent Fast Frequency..." 

j=1; % index
X_data =[]; % For storing the features (the system conditions) for the regression



iterations=1;
for P_loss=1.2e3:300:1.8e3 % Units: MW. NOTE: the model considers the power outage in area 1, see the Simulink model
    for Demand = 20e3:25e3:45e3
    for fraction_D=0.25:0.25:0.75
        P_D = Demand*[fraction_D (1-fraction_D)]; 
        D_prime = D*P_D';
    
%         % Consider H that gives RoCoF between 0.125Hz/s and 1Hz/s for the
%         % single-area equivalent:
%         H = P_loss/(2*RoCoF_max);
%         original_H = H;

        for R=P_loss-D_prime*0.5:500:P_loss-D_prime*0.5+4e3
            % Now get the value of k* in constraint (13) of paper "Stochastic 
            % Scheduling With Inertia-Dependent...", so that I get a value of R
            % close to the amount needed for barely complying nadir in single-area 
            PL_prime = P_loss-D_prime*Delta_fDB;
            syms k
            y = vpasolve(2*k/Td*log(2*k/(Td*D_prime*PL_prime+2*k)) ==...
                D_prime^2*(Delta_fmax-Delta_fDB)-D_prime*PL_prime,...
                k, [0 Inf]);
            % This last term makes sure that the solution for k is >0. Be careful if
            % there are more than 1 solution for the equation
            sol_k = double(y);

            H=sol_k/R; 
            
            original_H=H;
            original_R=R;
            
            if H<(P_loss/(2*RoCoF_max)) % If H is very low because R is quite high, make H at least comply with RoCoF, so that I don't consider samples that would violate the RoCoF constraint (even if they would comply with the nadir)
                H = (P_loss/(2*RoCoF_max));
                % Now recalculate R because I want just the right amount to
                % comply with the nadir of the COI:
                syms k
                y = vpasolve(2*k/Td*log(2*k/(Td*D_prime*PL_prime+2*k)) ==...
                    D_prime^2*(Delta_fmax-Delta_fDB)-D_prime*PL_prime,...
                    k, [0 Inf]);
                % This last term makes sure that the solution for k is >0. Be careful if
                % there are more than 1 solution for the equation
                sol_k = double(y);
                R=sol_k/H;                
            end

            for fraction_R=0.25:0.25:0.75

                PFR1 = fraction_R*R;
                PFR2 = (1-fraction_R)*R;

                % Now split the inertia between the areas:
                for fraction_H=0.25:0.125:0.75

                    H1 = fraction_H*H;
                    H2 = (1-fraction_H)*H;

                    % Run the simulation:
                    sim('TwoRegion_swing_NotSoFineTimeStep')  %runs the simulink model
                    %This Simulink file uses a very small fixed time-step 
                    %for the simulation, so that the simulation hits
                    %exactly the peaks of the oscillation. This makes the
                    %simulation slower, but more accurate.

                    %% Compute numerical Rocof from simulation samples
                    % I use the 7th and 1st samples, because using the 2nd
                    % and 1st samples can give numerical errors because of
                    % dividing by a very small number
                    % "Delta_f1.Time(2)-Delta_f1.Time(1)"
                    nadir2_check = abs(min(Delta_f2.Data));
                    i=1;
                    while nadir2_check >= Delta_fmax
                        % Now get the value of k* in constraint (13) of paper "Stochastic 
                        % Scheduling With Inertia-Dependent...", so that I get a value of R
                        % close to the amount needed for barely complying nadir in single-area 
                        syms k
                        y = vpasolve(2*k/Td*log(2*k/(Td*D_prime*PL_prime+2*k)) ==...
                            D_prime^2*(Delta_fmax-0.00625*i-Delta_fDB)-D_prime*PL_prime,...
                            k, [0 Inf]);
                        i=i+1;
                        % This last term makes sure that the solution for k is >0. Be careful if
                        % there are more than 1 solution for the equation
                        sol_k = double(y);

                        % Now increase H and R to meet this new "sol_k":
                        H = sol_k/R;

                        H1 = fraction_H*H;
                        H2 = (1-fraction_H)*H;

                        % Run the simulation:
                        sim('TwoRegion_swing_NotSoFineTimeStep')

                        nadir2_check = abs(min(Delta_f2.Data));

                    end
                    clear i

                    nadir2(j) = nadir2_check;
                    [~,t_nadir_index] = min(Delta_f2.Data);
                    t_nadir2(j) = Delta_f2.Time(t_nadir_index);

                    %% Extract the samples for the oscillations:
            %                     % First, get the single-area equivalent frequency 
            %                     % deviation "Delta_f", so that I can use it to extract 
            %                     % the oscillations.
            %                     %
            %                     % Use the analytical solution from eq. (6) in paper 
            %                     % "Stochastic Scheduling With Inertia-Dependent...":
            %                     %
            %                     % NOTE: this expression is only valid until Td, because
            %                     % after that PFR is constant instead of ramping up
            %                     single_area = @(t) -((P_loss/(D*P_D')...
            %                         + 2*(PFR1+PFR2)*(H1+H2)/(Td*(D*P_D')^2)) * (1-exp(-(D*P_D')*t/(2*(H1+H2))))...
            %                         - (PFR1+PFR2)*t/(Td*(D*P_D')));
            % 
            %                     single_area_samples = single_area(Delta_f1.Time);



                    cd('OneArea')
                    T_d=Td;
                    sim('maxRamp_FromWorkspace')
                    cd('..\')      

                    % If you want to double-check that the simulation is
                    % right, plot the result:
                    fig=figure(1);
                    plot(Delta_f1,'LineWidth',1.5)
                    hold on
                    plot(Delta_f2,'LineWidth',1.5)  
                    plot(Delta_f,'LineWidth',1.5)
                    plot([0 1],[0 -RoCoF_max],'--')
                    %plot(Delta_f1.Time,single_area_samples) 
                    hold off
                    axis([0 10 -0.8 0.1])
                    set(findall(fig,'-property','FontSize'),'FontSize',12)
                    xlabel('time (s)','FontSize',14)
                    ylabel('$$\Delta$$f (Hz)','Interpreter','latex','FontSize',16)
                    cd('Figures Freq')
                    print(['Fig' num2str(iterations)],'-dpng')
                    cd('..\')
                    close all

                    iterations=iterations+1;

                    t_until_nadir = Delta_f2.Time(1:t_nadir_index);
                    Delta_f2_until_nadir = Delta_f2.Data(1:t_nadir_index);
                    
                    Delta_f1_until_nadir_in_region2 = Delta_f1.Data(1:t_nadir_index);
                    
                    %integral1 = trapz(t_until_nadir,Delta_f2_until_nadir);
                    
                    % CHANGE SIGN of Delta_f for the integrals:
                    Delta_f2_until_nadir = -Delta_f2_until_nadir;
                    Delta_f1_until_nadir_in_region2 = -Delta_f1_until_nadir_in_region2;
                    
                    
                    % Single and double integrals for Delta_f2
                    first_integral_f2 = [];
                    for k=2:length(Delta_f2_until_nadir)
                        first_integral_f2 = [first_integral_f2 trapz(t_until_nadir(1:k),Delta_f2_until_nadir(1:k))];
                    end
                    double_integral_f2 = [];
                    for k=2:length(first_integral_f2)
                        double_integral_f2 = [double_integral_f2 trapz(t_until_nadir(1:k),first_integral_f2(1:k))];
                    end
                    
                    integral1(j) = first_integral_f2(end);
                    
                    % Single and double integrals for Delta_f1 (only the
                    % double integral is needed, but the single integral
                    % must be computed to compute the double)
                    first_integral_f1 = [];
                    for k=2:length(Delta_f1_until_nadir_in_region2)
                        first_integral_f1 = [first_integral_f1 trapz(t_until_nadir(1:k),Delta_f1_until_nadir_in_region2(1:k))];
                    end
                    double_integral_f1 = [];
                    for k=2:length(first_integral_f1)
                        double_integral_f1 = [double_integral_f1 trapz(t_until_nadir(1:k),first_integral_f1(1:k))];
                    end
                    
                    integral2(j) = (double_integral_f1(end) - double_integral_f2(end)); % Energy exported to the faulted region, Region 1

    %                 plot(t_until_nadir,Delta_f2_until_nadir)
    %                 hold on
    %                 plot([0 t_until_nadir(end)],[0 -0.8])

                    %% Save the system condition considered, i.e. the "features" for the regression:
                    PLoss(j) = P_loss;
                    H_singleArea(j) = H1+H2;
                    R_singleArea(j) = PFR1+PFR2;
                    Rocof_singleArea(j) = P_loss/(2*(H1+H2)); 

                    H_1(j) = H1;
                    H_2(j) = H2;
                    R1(j) = PFR1;
                    R2(j) = PFR2;
                    D_prime_1(j) = D(1)*P_D(1);
                    D_prime_2(j) = D(2)*P_D(2);

                    % Save in matrix form for the multivariate
                    % regression
                    X_data = [X_data; H1 H2 PFR1 PFR2 D(1)*P_D(1)...
                        D(2)*P_D(2) P_loss];

                    j=j+1;

                    clear single_area_samples
                    
                    H = original_H;
                    %PFR1 = fraction_R*R;
                    %PFR2 = (1-fraction_R)*R;

                end
            end
            R = original_R;
        end
    end
    end
end

X_data = X_data*1e-3;
X_data = [ones(size(X_data,1),1) X_data]; % Add the y-intercept feature

clear Delta_f1 Delta_f2 fractionH1 fractionR1 H H1 H2 i idx_max1 idx_max2...
    idx_min1 idx_min2 j k logsout P_loss PFR1 PFR1_function PFR2 PFR2_function...
    R single_area sol_k tout y D_prime_1 D_prime_2
save('Energy_analysis_output.mat') 


