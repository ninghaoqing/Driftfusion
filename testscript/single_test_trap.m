% close all
% clear all


% plot(t,g);
% load("temp.mat");
% load("data1_layer_Jiaxin_no_sp_10.mat");
load("data1_layer_Jiaxin_sp_5.mat");
% tspan=time;
y=zeros(1,length(time));

string1="data1_layer_Jiaxin_no_sp_"
% for i=10:1:10
%    string=string1+num2str(i)+".mat"
y0=[0 0];

% global g; 
% load(string);
% load('low_T.mat');
g=interface;
i=1;
k1=9e-13;
%  k2=2e7;  %pristine
% k2=1e9; % SP
% k3=0;
dt=time(2)-time(1);
dy=zeros(size(time));

kBT=0.259 %eV
dEVB=-280e-3;
NVB=1e20;
% kr=1e-10;
k2_new=NVB*exp(dEVB/kBT)*k1;
n0=1e18;

for j=2:1:length(time)
    dy(i,j)=k1*g(j-1)*(n0-y(i,j-1))*dt-k2_new*y(i,j-1)*dt;
    y(i,j)=y(i,j-1)+dy(i,j);
end
% plot(time*1e9,interface./max(interface))
% plot(time,y(i,:));
% plot(time,y(i,:));
% time=time(1:800);
% y=y(1:800);
% plot(time*1e9,y/max(y));
plot(time*1e9-90,y/max(y));
temp=y/max(y);
% hold on
% plot(time,g/max(g));

hold on 
% load("jiaxin_240k_no_sp.mat");
load("jiaxin_300k_sp.mat");
% plot(time_2,trap);
plot(time_2((117:230)),trap((117:230)))
%  hold on
%  plot(time,g/max(g));
% end

% t=time;
% tspan=time;
% [t,y]=ode45(@(t,y) ODE(t,y,time,g),tspan,y0);
% 
% plot(t,y(:,1));
% figure;
% plot(t,y(:,2));
