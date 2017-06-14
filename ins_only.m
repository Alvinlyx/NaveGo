function [ins_e] = ins_only(imu,gps, ref,att_mode, precision)
% this function is based on NaveGo
% input
% imu , IMU data structure
%     t: lx1 time vector ,(seconds)
%     fb: lx3 accelerations vector in body frame XYZ(m/s^2)
%     wb: lx3 turn rates vector in body frame XYZ(radians/s)
%     arw:lx3 angle random walks (rad/s/root-Hz)
%     vrw:lx3 angle velocity walks(radians/s)
%     gstd:lx3 gyros standard deviations(radians/s)
%     astd:lx3 accrs standard deviations (m/s^2)
%     gb_fix:lx3 gyros static biases or turn-on biases(radians/s)
%     ab_fix:lx3 accrs static biases or turn-in biases(m/s^2)
%     gb_drift:lx3 gyros dynamic biases(radians/s)
%     ab_drift:lx3 accrs dynamic biases(m/s^2)
%     gb_corr: 1x3 gyros correlation times (seconds)
%     ab_corr: 1x3 accrs correlation times (seconds)
%     gpsd: 1x3 gyros dynamic biases PSD(rad/s/root-Hz)
%     apsd: 1x3 accrs dynamic biases PSD (m/s^2/root-Hz)
%     freq: 1x1 sampling frequency(Hz)
%     ini_align:1x3 initial attitude at t(1)
%     ini_align_err: 1x3 initial attitude errors at t(1)
% 
%     att_mode, attitude mode string
%     'quaternion' :: attitude updated in quaternion format. Default ...
%         value. Not support the parameter  in this version 
%     'dcm' : attitude updated in Direct Cosine Matrix format.
% 
%     precision, finite number precision string, Not Support in the
%     current version 
%     double : default value, 64 bits
%     single : 32bits
% 
% 
%  Output :
%     ins_e , INS estimates data structure
%     t: Ix1 time vector (seconds).
%     roll: Ix1 roll (radians).
%     pitch: Ix1 pitch (radians)
%     yaw: Ix1 yaw (radians).
%     vel: Ix3 NED velocities (m/s).
%     lat: Ix1 latitude (radians).
%     lon: Ix1 longitude (radians).
%     h: Ix1 altitude (m).

    if nargin < 4 
        precision = 'single';
    end
    if nargin < 3
        att_mode = 'dcm';
    end 
    ti = imu.t;
    Mi = max(size(imu.t));
    Mg = (max(size(gps.t)));
    % preallocate memory for estimates
    roll_e = zeros(Mi,1);
    pitch_e = zeros(Mi,1);
    yaw_e  = zeros(Mi,1);
    vel_e  = zeros(Mi,3);
    h_e    = zeros(Mi,1);
    
    % constant matrices
    I = eye(3);
    Z = zeros(3);

    % kalman matrices for later analysis
    % pass 
    
    % initialize biases variables
    gb_drift = imu.gb_drift';
    ab_drift = imu.ab_drift';
    gb_fix = imu.gb_fix';
    ab_fix = imu.ab_fix';
if strcmp(precision, 'single')  % single precision

    ti = single(imu.t);
    tg = single(gps.t);
    % Preallocate memory for estimates
    roll_e  = single(zeros (Mi, 1));
    pitch_e = single(zeros (Mi, 1));
    yaw_e   = single(zeros (Mi, 1));
    vel_e   = single(zeros (Mi, 3));
    h_e     = single(zeros (Mi, 1));
    
    % Constant matrices
    I = single(eye(3));
    Z = single(zeros(3));
    
    % Kalman matrices for later analysis
    In = single(zeros(Mg, 6));        % Kalman filter innovations
    Pi = single(zeros(Mg, 441));    % Elements from a priori covariance matrix, Pi
    Pp = single(zeros(Mg, 441));   % Elements from a posteriori covariance matrix, Pp
    A  = single(zeros(Mg, 441));    % Elements from transition-state matrix, A
    Xi = single(zeros(Mg, 21));      % Evolution of Kalman filter a priori states, xi
    Xp = single(zeros(Mg, 21));     % Evolution of Kalman filter a posteriori states, xp
    B  = single(zeros(Mg, 12));      % Biases compensantions after Kalman filter correction
    x  = single([ zeros(1,9), imu.gb_fix, imu.ab_fix, imu.gb_drift, imu.ab_drift ]');  % Kalman filter error vector state
        
    % Initialize biases variables
    gb_drift = single(imu.gb_drift');
    ab_drift = single(imu.ab_drift');
    gb_fix   = single(imu.gb_fix');
    ab_fix   = single(imu.ab_fix');

    
    % Initialize estimates at tti=1
    roll_e (1)   = single(imu.ini_align(1));
    pitch_e(1) = single(imu.ini_align(2));
    yaw_e(1)  = single(imu.ini_align(3));
    vel_e(1,:)  = single(gps.vel(1,:));    
    h_e(1)       = single(gps.h(1));
    
else % double precision
    
    ti = (imu.t);
    tg = (gps.t);
    
    % Preallocate memory for estimates
    roll_e  = zeros (Mi, 1);
    pitch_e = zeros (Mi, 1);
    yaw_e   = zeros (Mi, 1);
    vel_e   = zeros (Mi, 3);
    h_e     = zeros (Mi, 1);
    
    % Constant matrices
    I = eye(3);
    Z = zeros(3);
    
    % Kalman matrices for later analysis
    In = zeros(Mg, 6);           % Kalman filter innovations
    Pi = zeros(Mg, 441);       % Elements from a priori covariance matrices, Pi
    Pp = zeros(Mg, 441);      % Elements from a posteriori covariance matrices, Pp
    A  = zeros(Mg, 441);       % Elements from transition-state matrices, A
    Xi = zeros(Mg, 21);         % Evolution of Kalman filter a priori states, xi
    Xp = zeros(Mg, 21);        % Evolution of Kalman filter a posteriori states, xp
    B  = zeros(Mg, 12);         % Biases compensantions after Kalman filter correction
    x  = [ zeros(1,9), imu.gb_fix, imu.ab_fix, imu.gb_drift, imu.ab_drift ]';  % Kalman filter error vector state
        
    % Initialize biases variables
    gb_drift = imu.gb_drift';
    ab_drift = imu.ab_drift';
    gb_fix   = imu.gb_fix';
    ab_fix   = imu.ab_fix';
    
    % Initialize estimates at tti = 1
    roll_e(1)  = imu.ini_align(1);
    pitch_e(1) = imu.ini_align(2);
    yaw_e(1)   = imu.ini_align(3);
    vel_e(1,:) = gps.vel(1,:);
    h_e(1)     = gps.h(1);
end

% Lat and lon cannot be set in single precision. They need full (double) precision.
lat_e    = zeros (Mi,1);
lon_e    = zeros (Mi,1);
lat_e(1) = double(ref.lat(1));
lon_e(1) = double(ref.lon(1));

DCMnb = euler2dcm([roll_e(1); pitch_e(1); yaw_e(1);]);
DCMbn = DCMnb';
qua   = euler2qua([roll_e(1) pitch_e(1) yaw_e(1)]);

% Initialize Kalman filter matrices
S.R  = diag([gps.stdv, gps.stdm].^2);
S.Q  = diag([imu.arw, imu.vrw, imu.gpsd, imu.apsd].^2);
S.Pp = diag([imu.ini_align_err, gps.stdv, gps.std, imu.gb_fix, imu.ab_fix, imu.gb_drift, imu.ab_drift].^2);

% UD filter matrices
% [Up, Dp] = myUD(S.P);
% dp = diag(Dp);

% Initialize matrices for INS/GPS performance analysis
%Pp(1,:) = reshape(S.Pp, 1, 441);
%B(1,:)  = [gb_fix', ab_fix', gb_drift', ab_drift'];
    % Lat and  lon cannot be set in single precision. They need
    % full precision
%     lat_e = zeros(Mi,1);
%     lon_e = zeros(Mi,1);
%     lat_e(1) = double(gps.lat(1));
%     lon_e(1) = double(gps.lon(1));
%     h_e(1) = double(gps.h(1));
%     vel_e(1,:) = gps.vel(1,:);
%     DCMnb = euler2dcm([roll_e(1); pitch_e(1);yaw_e(1)]);
%     DCMbn = DCMnb';
% 
%     qua = euler2qua([roll_e(1); pitch_e(1);yaw_e(1)]);

    % initialize matrices for INS performance analysis 
    % pass 
    % initialize estimates at tti = 1
    for  i=2:Mi       % Index for INS navigation 
        if (mod(i,10000) == 0), fprintf('. ');  end
        % Print a return on console every 200,000 INS executions
        if (mod(i,200000) == 0), fprintf('\n'); end 
        % INS period
        dti = ti(i) - ti(i-1);
        
        % Correct inertial sensors
%      wb_corrected = (imu.wb(i,:)' + gb_fix + gb_drift );
%      fb_corrected = (imu.fb(i,:)' + ab_fix + ab_drift );
%          wb_corrected = imu.wb(i,:)';%-gb_fix-gb_drift;
%          fb_corrected = imu.fb(i,:)';%-ab_fix-ab_drift;
%         % Attitude update
%         omega_ie_N = earthrate(lat_e(i-1), precision);
%         omega_en_N = transportrate(lat_e(i-1), vel_e(i-1,1), vel_e(i-1,2), h_e(i-1));
        
       % [qua_n, DCMbn_n, euler] = att_update(wb_corrected, DCMbn, qua,omega_ie_N, omega_en_N, dti, att_mode);
%         roll_e(i) = euler(1);
%         pitch_e(i)= euler(2);
%         yaw_e(i)  = euler(3);
        
       % DCMbn = DCMbn_n;
     %   qua = qua_n;
        
        % Gravity update
%         g = gravity(lat_e(i-1), h_e(i-1));
        
        % Velocity update
%         fn = (DCMbn_n * fb_corrected);
%         vel_n = vel_update(fn, vel_e(i-1,:), omega_ie_N, omega_en_N, g', dti); %
%         vel_e (i,:) = vel_n;
        vel_e(i,:) =  ref.vel(i,:);
        % Position update
        pos = pos_update([lat_e(i-1) lon_e(i-1) double(h_e(i-1))], double(vel_e(i,:)), double(dti) );
        lat_e(i) = pos(1);
        lon_e(i) = pos(2);
        h_e(i)   = pos(3);
       % Magnetic heading update
       % yawm_e(i) = hd_update (imu.mb(i,:), roll_e(i),  pitch_e(i), D)
    end
    ins_e.lat = lat_e;
    ins_e.lon = lat_e;
    ins_e.h = h_e;
    fprintf('\n');
