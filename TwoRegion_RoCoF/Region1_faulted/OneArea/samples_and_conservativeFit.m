clear
close all
clc

%%
% Comments:

%The name of the parameters in the simulink model has to be the same as 
% the name of the variables that are being modified in the loop.

%Use the block "to workspace" in Simulink for saving the output of a
%signal. The value will be saved in a variable with the same name as the
%one given to block "to workspace in the Simulink model. For example, in
%this case it's "Delta_f1" and "Delta_f2"

%%
%H = 4900; % Units: MW*s^2
D = 0.0000001e-2; % Units 1/Hz
P_D = 50e3; % Units: MW
%R = 1.7e3; % Units: MW
T_d = 10; % Units: s
P_loss = 1.8e3; % Units: MW. The model considers the power outage in area 1

%total_time = 300; % Units: s. Total time of the simulation

nadir_req = 0.8;
tol_nadir = 0.01;

%%
H_vector=3000:100:5500;
R_vector=1800:100:3500;

Xnadir = [];
for i=1:length(H_vector)
    for j=1:length(R_vector) 
        tic
        H = H_vector(i);
        R = R_vector(j);
        sim('FromWorkspace')  %runs the simulink model

        % Nadir:
        nadir = abs(min(Delta_f.Data));

        % if nadir is close to the acceptable limit, keep the
        % samples. If not, discard them:
        if (nadir<=nadir_req)&&(nadir>=nadir_req-tol_nadir)
            Xnadir = [Xnadir; H R];
        end
        
        toc
        clear nadir
    end
end
close_system('OneArea_test_1')

figure(1)
plot(Xnadir(:,1),Xnadir(:,2),'.b')
hold on

%% Constrained fit:
x = Xnadir(:,1); %reshape the data into a column vector
y = Xnadir(:,2);

% 'C' is not the Vandermonde matrix for 'x' in this case, as I am not
% fitting a polynomial but an inverse function (check the picture in folder
% "conservative above-all-points" for more info

V(:,2) = ones(length(x),1,class(x));
V(:,1) = 1./x;

C = V;

% 'd' is the vector of target values, 'y'.
d = y;

% Set that the fitted function must be above all data points: 
A = -C;
b = -d;

% No equality constraints:
Aeq = [];
beq = [];

p = lsqlin( C, d, A, b, Aeq, beq );

% Plot fitted data
x_plot = x(1):1:x(end);
plot(x_plot,1./x_plot*(P_loss^2*T_d/(4*nadir_req)),'r','linewidth',2)
plot(x_plot,p(1)*1./x_plot+p(2),'g','linewidth',2) 
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',12)
ylabel('R (MW)','FontSize',14)
xlabel('H (MW)','FontSize',14)
legend('Samples','Analytical constraint','Numerical-Conservative constraint')
print(figure(1),'-dpng', 'Fig2')
print(figure(1),'-depsc', 'Fig2')

%%
% Check that the constrained fit is right:
yhat2 = p(1)*1./x+p(2);
max_constrained = max(yhat2-y) % Should be higher than 0
min_constrained = min(yhat2-y) % Should be 0


        