clc;
clearvars;
close all;
% *******************************************************************

%% load and modify the data
load("COVID_STL.mat");
cases_STL2 = [cases_STL(1),diff(cases_STL)];
deaths_STL2 = [deaths_STL(1),diff(deaths_STL)];

measured_deaths = deaths_STL/POP_STL;
measured_infected = cases_STL/POP_STL;

plot(dates, [measured_infected;measured_deaths]);
% *******************************************************************
% From '2020-03-18' to '2021-06-23'
x0 = [1-7/POP_STL, 7/POP_STL, 0, 0];

infectious_rate = 0.0015;
immune_rate = 0.0921;
death_rate = 0.018;
recover_rate = 0.201;
rein_rate = 0.0005;

A = [1-infectious_rate  1-(recover_rate+immune_rate+death_rate) rein_rate   0;
     infectious_rate    recover_rate                            0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

B = zeros(4,1);
t = 1:67;

y = model(A,B,t,x0);
new_cases_model = cumsum([cases_STL(1)/POP_STL,(y(:,1)*infectious_rate).']);
new_cases_model = new_cases_model(1:67);

err1 = immse(y(:,4).',measured_deaths(1:67));
err2 = immse(new_cases_model,measured_infected(1:67));

figure;
plot(new_cases_model);
hold on;
plot(y(:,2:4));
plot(measured_infected(1:67));
plot(measured_deaths(1:67));
legend("new cases model","Model I","Model R","Model D","measured infected","measured deaths",'Location','northwest')

% *******************************************************************
% From '2021-06-30' to '2021-10-26'
x0 = y(size(y,1),:);

infectious_rate = 0.00155;
immune_rate = 0.092;
death_rate = 0.008;
recover_rate = 0.20;
rein_rate = 0.0005;

A = [1-infectious_rate  1-(recover_rate+immune_rate+death_rate) rein_rate   0;
     infectious_rate    recover_rate                            0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

t = (68:84)-67;

y1 = model(A,B,t,x0);
new_cases_model2 = new_cases_model(length(new_cases_model))+cumsum((y1(:,1)*infectious_rate).');

err3 = immse(y1(:,4).',measured_deaths(68:84));
err4 = immse(new_cases_model2,measured_infected(68:84));

figure;
plot(new_cases_model2);
hold on;
plot(y1(:,2:4));
plot(measured_infected(68:84));
plot(measured_deaths(68:84));
legend("new cases model","Model I","Model R","Model D","measured infected","measured deaths",'Location','northwest')
% *******************************************************************
% From '2021-10-27' to '2022-03-22'
x0 = y1(size(y1,1),:);

infectious_rate = 0.0044;
immune_rate = 0.092;
death_rate = 0.007;
recover_rate = 0.204;
rein_rate = 0.0005;


t = (85:105)-84;
A = [1-infectious_rate  1-(recover_rate+immune_rate+death_rate) rein_rate   0;
     infectious_rate    recover_rate                            0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

y2 = model(A,B,t,x0);
new_cases_model3 = new_cases_model2(length(new_cases_model2))+cumsum((y2(:,1)*infectious_rate).');

err5 = immse(y2(:,4).',measured_deaths(85:105));
err6 = immse(new_cases_model3,measured_infected(85:105));

figure;
plot(new_cases_model3);
hold on;
plot(y2(:,2:4));
plot(measured_infected(85:105));
plot(measured_deaths(85:105));
legend("new cases model","Model I","Model R","Model D","measured infected","measured deaths",'Location','northwest')
% *******************************************************************
% From '2021-10-27' to '2022-03-22'
x0 = y2(size(y2,1),:);

infectious_rate = 0.0015;
immune_rate = 0.092;
death_rate = 0.004;
recover_rate = 0.2;
rein_rate = 0.0005;


t = (106:158)-105;
A = [1-infectious_rate  1-(recover_rate+immune_rate+death_rate) rein_rate   0;
     infectious_rate    recover_rate                            0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

y3 = model(A,B,t,x0);
new_cases_model4 = new_cases_model3(length(new_cases_model3))+cumsum((y3(:,1)*infectious_rate).');

err7 = immse(y3(:,4).',measured_deaths(106:158));
err8 = immse(new_cases_model4,measured_infected(106:158));

figure;
plot(new_cases_model4);
hold on;
plot(y3(:,4));
plot(measured_infected(106:158));
plot(measured_deaths(106:158));
legend("new cases model","Model D","measured infected","measured deaths",'Location','northwest')
% *******************************************************************


%% Define functions
% Model function without input

function Y = model(A,B,t,x0)
size1 = size(t,2);
sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
Y = lsim(sys_sir_base,zeros(size1,1),t,x0);
end
