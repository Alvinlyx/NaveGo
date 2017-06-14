clc
clear 
close all 
load ref.mat
load imu2.mat
imu = imu2;
%% initinalize variables 
Mi  = max(size(ref.vel));
ti  = ref.t;
vel = ref.vel;
R2D = 180/pi;
% real 
lat_e = zeros(Mi,1);
lon_e = zeros(Mi,1);
h_e   = zeros(Mi,1);

lat_e(1) = ref.lat(1);
lon_e(1) = ref.lon(1);
h_e(1)   = ref.h(1);
% add noise
lat_noise = zeros(Mi,1);
lon_noise = zeros(Mi,1);
h_noise   = zeros(Mi,1);

lat_noise(1) = ref.lat(1);
lon_noise(1) = ref.lon(1);
h_noise(1)   = ref.h(1);

vel_noise     = zeros(Mi,3);
vel_noise(1,:)= ref.vel(1,:);

roll_noise  = zeros(Mi,1);
pitch_noise = zeros(Mi,1);
yaw_noise   = zeros(Mi,1);
roll_noise(1)  = imu.ini_align(1);
pitch_noise(1) = imu.ini_align(1);
yaw_noise(1)   = imu.ini_align(1);
% Initialize biases variables
    gb_drift = imu.gb_drift';
    ab_drift = imu.ab_drift';
    gb_fix   = imu.gb_fix';
    ab_fix   = imu.ab_fix'; 

%%
% real 
for i=2:Mi 
    dt = ti(i)-ti(i-1);
    pos = pos_update([lat_e(i-1,1) lon_e(i-1,1) h_e(i-1,1)] , vel(i-1,:), dt);
    lat_e(i) = pos(1);
    lon_e(i) = pos(2);
    h_e(i)   = pos(3);
end
% noise
att_mode = 'dcm';
precision = 'double';
qua   = euler2qua([roll_noise(1) pitch_noise(1) yaw_noise(1)]);
DCMnb = euler2dcm([roll_noise(1); pitch_noise(1); yaw_noise(1);]);
DCMbn = DCMnb';
for i=2:Mi
    dt = ti(i)-ti(i-1);

    wb = imu2.wb(i,:)'+ gb_fix + gb_drift ;
    fb = imu2.fb(i,:)'+ ab_fix + ab_drift;
    % attitude update
    omega_ie_N = earthrate(lat_noise(i-1), precision);
    omega_en_N = transportrate(lat_noise(i-1), vel_noise(i-1,1), ...
                               vel_noise(i-1,2), h_noise(i-1));
    [qua,DCMbn_n, euler]= att_update(wb,DCMbn,qua,omega_ie_N, ...
        omega_en_N, dt, att_mode);
    roll_noise(i)  = euler(1);
    pitch_noise(i) = euler(2);
    yaw_noise(i)   = euler(3);
    % Gravity update
    g = gravity(lat_noise(i-1), h_noise(i-1));
    % Velocity update
    fn = (DCMbn_n * fb);
    vel_n = vel_update(fn, vel_noise(i-1,:), omega_ie_N, omega_en_N, g', dt); %
    vel_noise(i,:) = vel_n;
    pos = pos_update([lat_noise(i-1,1) lon_noise(i-1,1) h_noise(i-1,1)] ...
                     , vel_noise(i,:), dt);
    
    lat_noise(i) = pos(1);
    lon_noise(i) = pos(2);
    h_noise(i)   = pos(3);
end
figure;
plot3(lat_e*R2D,lon_e*R2D,h_e);
hold on 
%plot3(ref.lat*R2D,ref.lon*R2D,ref.h);
plot3(lat_noise(1:2000)*R2D,lon_noise(1:2000)*R2D,h_noise(1:2000));
figure 
plot(lat_e*R2D,lon_e*R2D);
figure
plot(lat_noise*R2D,lon_noise*R2D);
