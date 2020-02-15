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
D = [1e-2 1e-2]; % Units 1/Hz
Demand = 25e3; % Units: MW
Td = 10; % Units: s. Delivery time of PFR.

X = 50; % Units: ohms. Line reactance
V = 345; % Units: kV (voltage of the transmission line) 
Line_term = 2*pi*V^2/X;

%%
% Some other system parameters, if you change these you will have to change
% the loops later in the code.
RoCoF_max = 1; % Relaxed RoCoF
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
    for fraction_D=0.25:0.25:0.75
        P_D = Demand*[fraction_D (1-fraction_D)]; 
        D_prime = D*P_D';
    
        % Consider H that gives RoCoF between 0.125Hz/s and 1Hz/s for the
        % single-area equivalent:
        H = P_loss/(2*RoCoF_max);
        original_H = H;

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

        R=sol_k/H; 
        if R<P_loss-D_prime*0.5
            R=P_loss-D_prime*0.5; % Make R respect the qss req
        end
        original_R = R;
        
        for fraction_R=0.25:0.25:0.75
        
            PFR1 = fraction_R*R;
            PFR2 = (1-fraction_R)*R;

            % Now split the inertia between the areas:
            for fraction_H=0.25:0.125:0.75

                H1 = fraction_H*H;
                H2 = (1-fraction_H)*H;

                % Run the simulation:
                sim('TwoRegion_swing_NotSoFineTimeStep.slx')  %runs the simulink model
                %This Simulink file uses a very small fixed time-step 
                %for the simulation, so that the simulation hits
                %exactly the peaks of the oscillation. This makes the
                %simulation slower, but more accurate.

                %% Compute numerical Rocof from simulation samples
                % I use the 7th and 1st samples, because using the 2nd
                % and 1st samples can give numerical errors because of
                % dividing by a very small number
                % "Delta_f1.Time(2)-Delta_f1.Time(1)"
                
                for index=1:length(Delta_f2.Time)/4
                    Rocof2_samples(index) = abs(Delta_f2.Data(index+10)-Delta_f2.Data(index))/(Delta_f2.Time(index+10)-Delta_f2.Time(index));
                end
                Rocof2_check=max(Rocof2_samples);
                i=1;
                while Rocof2_check >= RoCoF_max
                    H = P_loss/(2*(RoCoF_max-0.00625*i)); %Note that H is in MWs^2
                    i=i+1;


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

                    R=sol_k/H; 
                    if R<P_loss-D_prime*0.5
                        R=P_loss-D_prime*0.5; % Make R respect the qss req
                    end
                    PFR1 = fraction_R*R;
                    PFR2 = (1-fraction_R)*R;


                    H1 = fraction_H*H;
                    H2 = (1-fraction_H)*H;

                    % Run the simulation:
                    sim('TwoRegion_swing_NotSoFineTimeStep.slx')

                    for index=1:length(Delta_f2.Time)/4
                        Rocof2_samples(index) = abs(Delta_f2.Data(index+10)-Delta_f2.Data(index))/(Delta_f2.Time(index+10)-Delta_f2.Time(index));
                    end
                    Rocof2_check=max(Rocof2_samples);
                end
                clear i



                Rocof2(j) = Rocof2_check;
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
                plot([0 1],[0 -Rocof2_check],'--')
                plot([0 1],[0.2 -Rocof2_check+0.2],'--')
                plot([0 1],[0.4 -Rocof2_check+0.4],'--')
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

%                 % First area:
%                 [C_COIapprox,Constant_ofExponential_single_area,...
%                     exponent_single_area,C_oscillations,single_area_samples] =...
%                     COIapprox_terms_fromLaplace_f1([H1 H2],D,P_D,[PFR1 PFR2],Td,P_loss,Line_term,Delta_f1);
%                 % Remove the single-area frequency deviation
%                 Delta_f1.Data = Delta_f1.Data - single_area_samples;
%                 C1(j) = C_oscillations;
        
                % Second area:
                [C_COIapprox,Constant_ofExponential_single_area,...
                    exponent_single_area,C_oscillations,single_area_samples] =...
                    COIapprox_terms_fromLaplace_f2([H1 H2],D,P_D,[PFR1 PFR2],Td,P_loss,Line_term,Delta_f2);
                % Remove the single-area frequency deviation
                Delta_f2.Data = Delta_f2.Data - single_area_samples;
                C2(j) = C_oscillations;
        
%                 % If you want to double-check the oscillations, plot them:
%                 fig=figure(2);                
%                 plot(Delta_f1,'LineWidth',1.5)
%                 hold on
%                 plot(Delta_f2,'LineWidth',1.5)
%                 hold off                   
%                 set(findall(fig,'-property','FontSize'),'FontSize',12)
%                 xlabel('time (s)','FontSize',14)
        
        
                %% Estimate oscillation parameters
%                 %C1(j) = D(2)*P_D(2)*(PFR1*D(2)*P_D(2)-PFR2*D(1)*P_D(1))/((D(1)*P_D(1)+D(2)*P_D(2))^2*Td*Line_term);
%                 Delta_f1.Data = -Delta_f1.Data;
%                 C1(j) = -C1(j); 
%     
%                 [a1(j),A1(j),T1(j),phi1(j)] = estimate_oscillation(Delta_f1,C1(j),iterations);
%                 % Take the negative of the oscillations, because
%                 % function "estimate_oscillation" is better suited for 
%                 % a sine function, rather than a "-sin()" 
%                 %
%                 % I also need to input the -(-C) of the constant because as an input to the
%                 % funcion I have given the -Oscillations, so that the sine for this
%                 % faulted area was positive
%                 iterations = iterations+1;
        
                %C2(j) = D(1)*P_D(1)*(PFR1*D(2)*P_D(2)-PFR2*D(1)*P_D(1))/((D(1)*P_D(1)+D(2)*P_D(2))^2*Td*Line_term);
                [a2(j),A2(j),T2(j),phi2(j)] = estimate_oscillation(Delta_f2,C2(j),iterations);
                iterations = iterations+1;
        
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
                R = original_R;
                PFR1 = fraction_R*R;
                PFR2 = (1-fraction_R)*R;

            end
        end
    end
end

X_data = X_data*1e-3;
X_data = [ones(size(X_data,1),1) X_data]; % Add the y-intercept feature

clear Delta_f1 Delta_f2 fractionH1 fractionR1 H H1 H2 i idx_max1 idx_max2...
    idx_min1 idx_min2 j k logsout P_loss PFR1 PFR1_function PFR2 PFR2_function...
    R single_area sol_k tout y D_prime_1 D_prime_2
save('oscillation_analysis_output.mat') 


