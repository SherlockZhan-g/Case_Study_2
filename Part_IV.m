clc;
close all;
clearvars;

load("mockdata2023.mat");



w = 30;
cuminfections = cumsum(newInfections);
itrend = movmean(cuminfections,w);
dtrend = movmean(cumulativeDeaths,w);
plot(dtrend.')

figure;
plot([cumulativeDeaths;newInfections].');
% *******************************************************************
x0 = [1-newInfections(1)-cumulativeDeaths(1),newInfections(1),0,0,cumulativeDeaths(1),0];

v_rate = 0;
inf_rate1 = 0.005;
inf_rate2 = 0;
immune_rate = 0.3;
death1 = 0.06;
death2 = 0;
recover1 = 0.2;
recover2 = 0;
rein_rate = 0.05;

A = [1-inf_rate1-v_rate recover1                      0                 rein_rate   0 0
     inf_rate1          1-recover1-death1-immune_rate 0                 0           0 0
     0                  0                             1-recover2-death2 0           0 inf_rate2
     0                  immune_rate                   0                 1-rein_rate 0 0 
     0                  death1                        death2            0           1 0
     v_rate             0                             recover2          0           0 1-inf_rate2 ];

B = zeros(6,1);
t = 1:100;


y = modelVaccine(A,B,t,x0);


model_newinf = [x0(2)+x0(3),(y(:,1)*inf_rate1 + y(:,6)*inf_rate2).'];
model_newinf = model_newinf(1:100);

 
figure
tiledlayout(1,2)

nexttile
plot(model_newinf);
hold on;
plot(newInfections(1:100).');
hold off;
nexttile
plot(y(:,5));
hold on;
plot(cumulativeDeaths(1:100).');
% *******************************************************************
v_rate = 0;
inf_rate1 = 0.01;
inf_rate2 = 0;
immune_rate = 0.3;
death1 = 0.06;
death2 = 0;
recover1 = 0.2;
recover2 = 0;
rein_rate = 0.05;

A = [1-inf_rate1-v_rate recover1                      0                 rein_rate   0 0
     inf_rate1          1-recover1-death1-immune_rate 0                 0           0 0
     0                  0                             1-recover2-death2 0           0 inf_rate2
     0                  immune_rate                   0                 1-rein_rate 0 0 
     0                  death1                        death2            0           1 0
     v_rate             0                             recover2          0           0 1-inf_rate2 ];

B = zeros(6,1);
t = (101:135)-100;
x0 = A*y(100,:).';

y1 = modelVaccine(A,B,t,x0);


model_newinf2 = model_newinf(100)+(y1(:,1)*inf_rate1 + y1(:,6)*inf_rate2).';

 
figure
tiledlayout(1,2)

nexttile
plot(model_newinf2);
hold on;
plot((newInfections(101:135).'));
hold off;
nexttile
plot(y1(:,5));
hold on;
plot(cumulativeDeaths(101:135).');
% *******************************************************************
v_rate = 0.3;
inf_rate1 = 0.0055;
inf_rate2 = 0.001;
immune_rate = 0.5;
death1 = 0.085;
death2 = 0.001;
recover1 = 0.215;
recover2 = 0.615;
rein_rate = 0.3;

A = [1-inf_rate1-v_rate recover1                      0                 rein_rate   0 0
     inf_rate1          1-recover1-death1-immune_rate 0                 0           0 0
     0                  0                             1-recover2-death2 0           0 inf_rate2
     0                  immune_rate                   0                 1-rein_rate 0 0 
     0                  death1                        death2            0           1 0
     v_rate             0                             recover2          0           0 1-inf_rate2 ];


x0 = A*y1(35,:).';
B = zeros(6,1);
t = (136:400)-135;


y2 = modelVaccine(A,B,t,x0);


model_newinf3 = model_newinf2(15)+(y2(:,1)*inf_rate1 + y2(:,6)*inf_rate2).';

 
figure
tiledlayout(1,2)

nexttile
plot(model_newinf3);
hold on;
plot(newInfections(136:400).');
hold off;
nexttile
plot(y2(:,5));
hold on;
plot(cumulativeDeaths(136:400).');
% *******************************************************************
%% Define functions
function Y = model(A,B,t,x0)
size1 = size(t,2);
sys_sir_base = ss(A,B,eye(4),zeros(4,1),1);
Y = lsim(sys_sir_base,zeros(size1,1),t,x0);
end

function Y = modelVaccine(A,B,t,x0)
size1 = size(t,2);
sys_sir_base = ss(A,B,eye(6),zeros(6,1),1);
Y = lsim(sys_sir_base,zeros(size1,1),t,x0);
end