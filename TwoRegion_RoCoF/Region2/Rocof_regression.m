clearvars
close all
clc

%%
load('oscillation_analysis_output')

%% Use a multivariate linear regression 
X_over_Htotal = X_data./(2*H_singleArea'*1e-3); % Divided by H_singleArea so that the final Rocof constraint is linear

opts = optimset('Display','off'); % Turn off the output message of "lsqlin"

% y1 = (A1.*(2*pi./T1))'; % A*omega
% theta1 = regress(y1,X_over_Htotal)

y2 = (A2.*(2*pi./T2))';
theta2 = regress(y2,X_over_Htotal)


%% Constrained fit (conservative above-all points):
% This regression makes sure that frequency security is guranteed in ALL
% CASES.
%
% The idea is to run a CONSTRAINED least-squares linear regression. The
% constraint makes sure that the regression is above all data points
% (therefore it estimates a conservative RoCoF, possibly higher than the
% actual RoCoF for that system condition).
%
% Some info here: https://uk.mathworks.com/matlabcentral/answers/94272-how-do-i-constrain-a-fitted-curve-through-specific-points-like-the-origin-in-matlab
%
% See also the picture explanation_conservative_regression.jpg

% Set that the regression must be above all data points: 
A = -X_over_Htotal;
%b = -y1;

% No equality constraints:
Aeq = [];
beq = [];

% theta1 = lsqlin(X_over_Htotal,y1,A,b,Aeq,beq,[],[],[],opts)

% % Check that the 2 extreme points are positive
% max(X_over_Htotal)*theta1-min(y2)
% min(X_over_Htotal)*theta1-min(y1)

%% Do it for area 2 now
% 
% % Set that the fitted polynomial must be above all data points: 
% b = -y2;
% 
% theta2 = lsqlin(X_over_Htotal,y2,A,b,Aeq,beq,[],[],[],opts)

% NOTE: for area 2 (the non-faulted area) it's not necessary to use the
% conservative regression, since the


%%
%save('Rocof_reg_output.mat','X_over_Htotal','y1','y2','theta1','theta2')
%save('Rocof_reg_output.mat','X_over_Htotal','y1','theta1')
save('Rocof_reg_output.mat','X_over_Htotal','y2','theta2')

