clc;
clearvars;
close all;
% *******************************************************************

%% load and modify the data
load("COVID_STL.mat");
cases_STL2 = [cases_STL(1),diff(cases_STL)];
deaths_STL2 = [deaths_STL(1),diff(deaths_STL)];

measured_deaths = deaths_STL/POP_STL;
measured_infected = cases_STL/POP_STL-measured_deaths;


plot(dates, [measured_infected;measured_deaths]);

% *******************************************************************
% about first 100 days.
x0 = [1-7/POP_STL, 7/POP_STL, 0, 0];

A = [0.9997 0.5    0.0005 0;
     0.0003 0.4    0      0;
     0      0.02   0.9995 0;
     0      0.08   0      1];

B = zeros(4,1);
t = 1:14;

y = model(A,B,t,x0);
figure;
plot(measured_infected(1:14));
hold on;
plot(y(:,[2,4]));
hold on;
plot(measured_deaths(1:14));
legend("Model I","Model D","Measured infected","measured deaths")
% *******************************************************************
% before 6/30/21
x0 = [1-7/POP_STL, 7/POP_STL, 0, 0];
A = [0.999 0.7   0.0005  0;
     0.001 0.2   0       0;
     0     0.05  0.9995  0;
     0     0.05  0       1];
B = zeros(4,1);
t = 1:67;

y = model(A,B,t,x0);
figure;
plot(cumsum(y(:,2)));
hold on;
plot(y(:,4));
plot(measured_infected(1:67));
plot(measured_deaths(1:67));
legend("Model I","Model D","Measured infected","measured deaths")
% *******************************************************************
x0 = [1-7/POP_STL, 7/POP_STL, 0, 0];

A = [0.9997 0.5    0.0005 0;
     0.0003 0.4    0      0;
     0      0.02   0.9995 0;
     0      0.08   0      1];
B = zeros(4,1);
t = 1:67;

new_cases_model = 



% *******************************************************************

%% Define functions
% Model function without input

function Y = model(A,B,t,x0)
size1 = size(t,2);
sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
Y = lsim(sys_sir_base,zeros(size1,1),t,x0);
end