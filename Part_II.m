clc;
clearvars;
close all;
% *******************************************************************

%% load and modify the data
load("COVID_STL.mat");
cases_STL2 = [cases_STL(1),diff(cases_STL)];
deaths_STL2 = [deaths_STL(1),diff(deaths_STL)];

measured_deaths = deaths_STL;
measured_infected = cases_STL;

plot(diff(measured_infected));


% *******************************************************************
%% Pre-Delta Period
% From '2020-03-18' to '2021-06-23'
x0 = [POP_STL-7, 7, 0, 0];

infectious_rate = 0.0015;
immune_rate = 0.39;
death_rate = 0.016;
recover_rate = 0.3;
rein_rate = 0.0005;

A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

B = zeros(4,1);
t = 1:67;

y = model(A,B,t,x0);
new_cases_model = cumsum([cases_STL(1),(y(:,1)*infectious_rate).']);
new_cases_model = new_cases_model(1:67);

err1 = immse(y(:,4).',measured_deaths(1:67));
err2 = immse(new_cases_model,measured_infected(1:67));

figure
tiledlayout(1,2)

nexttile
plot(new_cases_model);
hold on;
plot(measured_infected(1:67));
hold off;
nexttile
plot(y(:,4));
hold on;
plot(measured_deaths(1:67));


% *******************************************************************
%% Delta-Period
% From '2021-06-30' to '2021-10-26'

infectious_rate = 0.0017;
immune_rate = 0.35;
death_rate = 0.0075;
recover_rate = 0.35;
rein_rate = 0.005;

A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

t = (68:84)-67;
x0 = A*y(size(y,1),:).';
y1 = model(A,B,t,x0);
new_cases_model2 = new_cases_model(length(new_cases_model))+cumsum((y1(:,1)*infectious_rate).');

err3 = immse(y1(:,4).',measured_deaths(68:84));
err4 = immse(new_cases_model2,measured_infected(68:84));

figure
tiledlayout(1,2)

nexttile
plot(new_cases_model2);
hold on;
plot(measured_infected(68:84));
hold off;
nexttile
plot(y1(:,4));
hold on;
plot(measured_deaths(68:84));
% *******************************************************************
%% Omicron Period
% From '2021-10-27' to '2022-03-22'
infectious_rate = 0.0044;
immune_rate = 0.25;
death_rate = 0.005;
recover_rate = 0.4;
rein_rate = 0.05;


t = (85:105)-84;
A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];
x0 = A*y1(size(y1,1),:).';
y2 = model(A,B,t,x0);
new_cases_model3 = new_cases_model2(length(new_cases_model2))+cumsum((y2(:,1)*infectious_rate).');

err5 = immse(y2(:,4).',measured_deaths(85:105));
err6 = immse(new_cases_model3,measured_infected(85:105));

figure
tiledlayout(1,2)

nexttile
plot(new_cases_model3);
hold on;
plot(measured_infected(85:105));
hold off;
nexttile
plot(y2(:,4));
hold on;
plot(measured_deaths(85:105));
% *******************************************************************
%% Post-Omicron Period
% From '2021-10-27' to '2022-03-22'

infectious_rate = 0.0017;
immune_rate = 0.2;
death_rate = 0.0045;
recover_rate = 0.5;
rein_rate = 0.05;


t = (106:158)-105;
A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];

x0 = A*y2(size(y2,1),:).';
y3 = model(A,B,t,x0);
new_cases_model4 = new_cases_model3(length(new_cases_model3))+cumsum((y3(:,1)*infectious_rate).');

err7 = immse(y3(:,4).',measured_deaths(106:158));
err8 = immse(new_cases_model4,measured_infected(106:158));

figure
tiledlayout(1,2)

nexttile
plot(new_cases_model4);
hold on;
plot(measured_infected(106:158));
hold off;
nexttile
plot(y3(:,4));
hold on;
plot(measured_deaths(106:158));
% *******************************************************************
x = [y;y1;y2;y3];
cases_model = [new_cases_model,new_cases_model2,new_cases_model3,new_cases_model4];
measured = [measured_infected;measured_deaths].';
figure
plot(x(:,4));
hold on;
plot(measured_deaths);
legend("S","I","R","D");
figure;
plot(cases_model);
hold on;
plot(measured_infected);
% *******************************************************************
infectious_rate = 0.0044;
immune_rate = 0.25;
death_rate = 0.005;
recover_rate = 0.4;
rein_rate = 0.05;


t = (85:105)-84;
A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];
x0 = A*y1(size(y1,1),:).';
y2 = model(A,B,t,x0);
new_cases_model3 = new_cases_model2(length(new_cases_model2))+cumsum((y2(:,1)*infectious_rate).');

netOcase1 = new_cases_model3(21)-new_cases_model3(1);
netOdeath1 = y2(21,4)-y2(1,4);


infectious_rate = 0.003;
immune_rate = 0.25;
death_rate = 0.005;
recover_rate = 0.4;
rein_rate = 0.05;


t = (85:105)-84;
A = [1-infectious_rate  recover_rate                            rein_rate   0;
     infectious_rate    1-(recover_rate+immune_rate+death_rate) 0           0;
     0                  immune_rate                             1-rein_rate 0;
     0                  death_rate                              0           1];
x0 = A*y1(size(y1,1),:).';
y2 = model(A,B,t,x0);
new_cases_model3 = new_cases_model2(length(new_cases_model2))+cumsum((y2(:,1)*infectious_rate).');

netOcase2 = new_cases_model3(21)-new_cases_model3(1);
netOdeath2 = y2(21,4)-y2(1,4);


percent1 = (netOcase2-netOcase1)/netOcase1;
percent2 = (netOdeath2-netOdeath1)/netOdeath1;

% *******************************************************************
%% Define functions
% Model function without input

function Y = model(A,B,t,x0)
size1 = size(t,2);
sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
Y = lsim(sys_sir_base,zeros(size1,1),t,x0);
end
