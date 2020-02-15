Rocof_osci1_regression = X_over_Htotal*theta1;
%X_data = [X_data; H1 H2 PFR1 PFR2 P_loss];
%X_data = [ones(size(X_data,1),1) X_data]; % Add the y-intercept feature
H_total=sum(X_data(:,2:3),2);
Ploss = X_data(:,end);
Rocof_COI = Ploss./(2*H_total);
                    

Rocof_estimated1 = Rocof_osci1_regression+Rocof_COI;

diff=Rocof_estimated1'-Rocof1;

%save('diff.mat','diff')

% Although "diff" has some slightly negative values, I believe this is due
% to calculating the numerical Rocof using the first 7 samples, which
% introduces some numerical errors
%   - Another possible reason for these negative values is that the Rocof
%   related to the COI is slightly higher than the one that comes from
%   Fei's Rocof constraint 