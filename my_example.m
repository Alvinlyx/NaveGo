clc
clear 
close all 
load ref.mat
% no errors 
%% initinalize variables 
Mi = max(size(ref.vel));
ti = ref.t;
vel = ref.vel;
R2D = 180/pi;
lat_e  = zeros(Mi,1);
lon_e = zeros(Mi,1);
h_e = zeros(Mi,1);
lat_e(1) =ref.lat(1);
lon_e(1) = ref.lon(1);
h_e(1) = ref.h(1);
%%
for i=2:Mi 
    dt = ti(i)-ti(i-1);
    pos = pos_update([lat_e(i-1,1) lon_e(i-1,1) h_e(i-1,1)] , vel(i-1,:), dt);
    lat_e(i)  =pos(1);
    lon_e(i) = pos(2);
    h_e(i)   = pos(3);
end 
plot3(lat_e*R2D,lon_e*R2D,h_e);
hold on 
plot3(ref.lat*R2D,ref.lon*R2D,ref.h);