% input_cs v = 'Input_files/3_layer_unit_test.csv';
close all
input_csv = 'Input_files/spiro_mapi_tio2_jiaxin_mob.csv';
% input_csv = 'Input_files/1_layer_Jiaxin_no_sp_4.csv';
par = pc(input_csv);
% par.push_time=1200;
% par.pushed=false;
% pc.ntrap=1e24;
% par.RelTol_vsr = 0.01; 
soleq = equilibrate(par);
% sol.par=par;
% sol.u=0;
%  JVsol=doLightPulse(soleq.el, 1, 1e-6, 100, 10, false ,false);
JVsol=doLightPulse(soleq.ion,1e5, 1000e-9, 2000, 1e-9,2e-9, false, false);
% JVsol=doLightPulse(soleq.el,0e-3, 10, 2000, 100, false, false);
% JVsol=doJV(soleq.el, 1e8,1 000,0.1,0,0,0.1,3);                 %CW

% dfplot.Jt(JVsol.ill.f,1e-2)
dfplot.Jt(JVsol,1)
% hold oncpe);
% figure;
% imagesc(JVsol.par.gx1)
% figure;
% imagesc(JVsol.par.gx2)
g1_fun=fun_gen(JVsol.par.g1_fun_type);
g2_fun = fun_gen(JVsol.par.g2_fun_type);
g1t=g1_fun(JVsol.par.g1_fun_arg,JVsol.t');
g2t=g2_fun(JVsol.par.g2_fun_arg,JVsol.t');
gxt1=g1t*JVsol.par.gx1;
gxt2=g2t*JVsol.par.gx2;
  gxt3 = zeros(size(gxt2));
%  ntrap=load('trap.mat');
%     ntrap=ntrap.y;
%     push_time=200;
% gxt3(push_time-10:push_time+10,1:end)=ntrap(push_time).*1e20;

% figure;
imagesc(gxt2);
figure;
imagesc(JVsol.x.*1e4,JVsol.t.*1e9,gxt2);
title('Excitation');
xlabel('x(um)');
ylabel('Time(ns)');
% hold on 
% plot(JVsol.t,temp/1e2);
figure;
 plot(JVsol.t.*1e9,JVsol.u(:,end,2));
 xlabel('time(ns)');
 ylabel('a.u');
 title('n(t)');
figure;
imagesc(JVsol.x.*1e4,JVsol.t.*1e9,JVsol.u(:,:,1));
title('Voltage(V)');
xlabel('x(um)');
ylabel('Time(ns)');
figure;
imagesc(JVsol.x.*1e4,JVsol.t.*1e9,JVsol.u(:,:,2));
title('n densities');
xlabel('x(um)');
ylabel('Time(ns)');
figure;
imagesc(JVsol.x.*1e4,JVsol.t.*1e9,JVsol.u(:,:,3));
title('p densities');
xlabel('x(um)');
ylabel('Time(ns)');
figure;
plot(JVsol.t.*1e9,JVsol.u(:,end,2));
% plot(JVsol.t*1e6,navg);
xlabel('Time(ns)');
ylabel('p densities (cm-3)');

figure;
 plot(JVsol.x.*1e4,JVsol.u(end,:,3));
 xlabel('X');
 ylabel('p');
% [Ecb, Evb, Efn, Efp,nt,pt]=dfana.calcTrap(JVsol);

figure;
points=1240;
navg=zeros(size(JVsol.u(:,end,3)));
for i=1:1:points
    navg=navg+JVsol.u(:,i,3);
end
navg=navg/points;
plot(JVsol.t*1e9,navg);
xlabel('Time(ns)');
ylabel('p densities (cm-3)');
p=JVsol.u(:,end,3);
time=JVsol.t;
[J, j, xmesh] = dfana.calcJ(JVsol, "sub");
dt=5e-10;
jint=sum(J.p(:,end)*dt);
[Ecb, Evb, Efn, Efp]=dfana.calcEnergies(JVsol);

time=time;
close all
plot(JVsol.t.*1e9,JVsol.u(:,end,3));
interface=JVsol.u(:,end,3);