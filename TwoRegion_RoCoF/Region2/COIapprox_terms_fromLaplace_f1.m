function [C_COIapprox,Constant_ofExponential_single_area,...
    exponent_single_area,C_oscillations,single_area_samples] =...
    COIapprox_terms_fromLaplace_f1(H,D,P_D,R,Td,P_loss,VV_times2pi_overXline,Delta_f1)

    warning off
    
    % First solve the 2-area system using the symbolic toolbox and the
    % Laplace transform:
    
    syms H1 D1 P1 R1 T_d 
    syms H2 D2 P2 R2
    syms P_Loss Transfer
    syms f1(t) f2(t) s

    syms t v

    df1(t) = diff(f1(t), t);
    intf1(t) = int(sym('f1(v)'),v,0,t); % I use "v" because the variable of integration 
                                        % and the extreme of integration cannot 
                                        % be the same.
                                        % I took this idea for the integral from:
                                        % http://stackoverflow.com/questions/5284457/how-to-find-laplace-of-intxt-t-0-t-in-matlab
    df2(t) = diff(f2(t),t);
    intf2(t) = int(sym('f2(v)'),v,0,t);

    % Define PFR:
    PFR1(t) = R1/T_d*t;
    PFR2(t) = R2/T_d*t;

    % Define the integro-differential eqs.:
    eq1(t) = 2*H1*df1(t) + D1*P1*f1(t) - PFR1(t) + P_Loss + Transfer*(intf1(t)-intf2(t));
    eq2(t) = 2*H2*df2(t) + D2*P2*f2(t) - PFR2(t) + Transfer*(intf2(t)-intf1(t));

    % Compute the Laplace transforms:
    F1(t) = laplace(eq1,t,s);
    F2(t) = laplace(eq2,t,s);
    syms LF1 LF2
    NF1 = subs(F1(t),{laplace(f1(t),t,s),laplace(f2(t),t,s)},{LF1,LF2});
    NF2 = subs(F2(t),{laplace(f1(t),t,s),laplace(f2(t),t,s)},{LF1,LF2});

    % Substitute the initial condition:
    NF1 = subs(NF1,f1(0),0);
    NF2 = subs(NF2,f2(0),0);

    % Substitute the values of the parameters:
    NF1 = subs(NF1,{H1,D1,P1,R1,T_d,P_Loss,Transfer,f1(0)}, ...
          {H(1),D(1),P_D(1),R(1),Td,P_loss(1),VV_times2pi_overXline,0});
    NF2 = subs(NF2,{H2,D2,P2,R2,T_d,Transfer,f2(0)}, ...
          {H(2),D(2),P_D(2),R(2),Td,VV_times2pi_overXline,0});

    % Collect items:
    NF1 = collect(NF1,LF1);
    NF2 = collect(NF2,LF2);

    % Solve the set of eqs defined by the Laplace transforms of the
    % differential eqs:
    [LF1, LF2] = solve(NF1, NF2, LF1, LF2);
    LF1 = partfrac(LF1);
    LF2 = partfrac(LF2);

    % Compute the inverse Laplace transform:
    f1 = ilaplace(LF1, s, t);
    f1 = vpa(f1,7);
    f2 = ilaplace(LF2, s, t);
    f2 = vpa(f2,7);

    %     % If you want to visualize the solution of the inverse Laplace,
    %     % uncomment these lines: 
    %     LF1_visualize = vpa(LF1,3)
    %     LF2_visualize = vpa(LF2,3)
    %     f1 = vpa(f1,7)
    %     f2 = vpa(f2,7)

    
    % Now that I have the solution from the inverse Laplace, get the values
    % of the parameters that I need:

    solution_terms = children(f1);
    
    C_total = double(solution_terms(end));
    
    t=0;
    Constant_ofExponential_single_area = double(subs(solution_terms(2)));
    t=1;
    exponent_single_area = log(double(subs(solution_terms(2)))/Constant_ofExponential_single_area);
    
    C_oscillations = C_total + Constant_ofExponential_single_area;
    C_COIapprox = -Constant_ofExponential_single_area;
    
    % Takes sample of the single-area solution, so that I can substract them
    % from the multi-area samples in order to obtain samples of just the
    % oscillations:
    t=1;
    Constant_ofRamp_single_area = double(subs(solution_terms(1)));
    t = Delta_f1.Time;
    single_area_samples = Constant_ofRamp_single_area*t+...
        Constant_ofExponential_single_area*exp(exponent_single_area*t)+C_COIapprox;
    
    warning on
    
end

