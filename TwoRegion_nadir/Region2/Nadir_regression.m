clearvars
close all
clc

%%
load('Energy_analysis_output')

interval_extremes = [3 3.35];
first_interval = (t_nadir2>interval_extremes(1)).*(t_nadir2<interval_extremes(2));


% Separate the data by nadir time-intervals:

X_data_temp = [X_data integral1' integral2' nadir2' t_nadir2'];
X_data_temp = X_data_temp.*(first_interval'*ones(1,size(X_data_temp,2)));
X_data_temp = X_data_temp(any(X_data_temp,2),:); % Remove rows with all zeros, https://uk.mathworks.com/matlabcentral/answers/40390-remove-rows-with-all-zeros
X_data_1st_interval = X_data_temp(:,1:end-4);
integral1_1st_interval = X_data_temp(:,end-3);
integral2_1st_interval = X_data_temp(:,end-2);
nadir2_1st_interval = X_data_temp(:,end-1);
t_nadir2_1st_interval = X_data_temp(:,end);


X_data=X_data_1st_interval;
integral1=integral1_1st_interval';
integral2=integral2_1st_interval';
nadir2=nadir2_1st_interval';
t_nadir2=t_nadir2_1st_interval';


y1 = integral1'; % Labels for the first type of regression (regression for integral1)
y2 = integral2'; % Labels for the second type of regression (regression for integral2)

%% Use a multivariate linear regression for each time-interval

opts = optimset('Display','off'); % Turn off the output message of "lsqlin"

% Non-conservative regressions:
%
% theta_integral1 = regress(y1,X_data);
% theta_integral2 = regress(y2,X_data);
% 
% integral1_regression = X_data*theta_integral1;
% integral2_regression = X_data*theta_integral2;


% CONSTRAINED FIT (BELOW all data points, so that the energy "injected" 
% from the integral is underestimated)
% Set that the regression must be below all data points: 
A = X_data;
b = y1;
% No equality constraints:
Aeq = [];
beq = [];

theta_conservative_integral1 = lsqlin(X_data,y1,A,b,Aeq,beq,[],[],[],opts)

integral1_regression_CONSERVATIVE = X_data*theta_conservative_integral1;


% CONSTRAINED FIT (BELOW all data points):
% Set that the regression must be below all data points: 
A = X_data;
b = y2;
% No equality constraints:
Aeq = [];
beq = [];

theta_conservative_integral2 = lsqlin(X_data,y2,A,b,Aeq,beq,[],[],[],opts)
integral2_regression_CONSERVATIVE = X_data*theta_conservative_integral2;



%%
% %save('Rocof_reg_output.mat','X_data','y1','y2','theta1','theta2')
% save('Nadir_reg_output.mat','X_data','y1','theta1')

