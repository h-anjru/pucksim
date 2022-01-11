function puckPlan(base,begin,dir1,dir2,n_lines,length,spacing,...
    h,calib_speed,speed,file)
% PUCKPLAN Plan serpentine UAS data collect mission
%
%          base: 1x2 base of operations E,N (UTM 17N, m)   
%         begin: 1x2 beginning point in E,N (UTM 17N, m)
%          dir1: direction of first flight line [deg]
%          dir2: direction of first turn (string 'l' or 'r')
%       n_lines: number of flight lines
%        length: length of flight line
%       spacing: flight line spacing
%             h: flying height above ground [m]
%   calib_speed: forward speed during calibration [m/s]
%         speed: forward speed during scanning [m/s]
%          file: filename for KML (no extension)
%
%   Output is an array of mission waypoints in WGS84 lon/lat in dd.ddddd
%   format (precise to 5 places after decimal, roughly 1-m precision). 
%   Also writes KML file.
%   
%   To omit lidar calibration pattern, enter calib_speed of zero.

% initialize some vars
    waypoints_utm = zeros(2*n_lines,2);
    waypoints_utm(1,:) = begin;
    home_speed = 10; % speed from/to base [m/s]

% define the "turn" bearing
    if dir2 == 'l'
        turn = dir1 - 90;
    elseif dir2 == 'r'
        turn = dir1 + 90;
    end
       
% define calibration pattern (double reverse bowtie pattern)
if calib_speed ~= 0
    calib = zeros(9,2); % initialize
    calib_length = 50; % length of parallels [m]
    calib(1,:) = begin + [2 2]; % slight offset to help appearance in KML
    calib(2,:) = calib(1,:) + ...
        calib_length*[sind(turn) cosd(turn)];
    calib(3,:) = calib(1,:) + ...
        calib_length*[sind(dir1) cosd(dir1)];
    calib(4,:) = calib(2,:) + ...
        calib_length*[sind(dir1) cosd(dir1)];
    calib(5,:) = calib(1,:);
    calib(6,:) = calib(4,:);
    calib(7,:) = calib(3,:);
    calib(8,:) = calib(2,:);
    calib(9,:) = calib(1,:);
else
    calib = [];
end

% populate waypoints (UTM)
    for ii = 2:2*n_lines
        if mod(ii,4) == 2
            waypoints_utm(ii,:) = waypoints_utm(ii-1,:) + ...
                length*[sind(dir1) cosd(dir1)];
        elseif mod(ii,4) == 1 || mod(ii,4) == 3
            waypoints_utm(ii,:) = waypoints_utm(ii-1,:) + ...
                spacing*[sind(turn) cosd(turn)];
        elseif mod(ii,4) == 0
            waypoints_utm(ii,:) = waypoints_utm(ii-1,:) + ...
                length*[sind(dir1+180) cosd(dir1+180)];
        end
    end
    
% combine all points
    all_pts = [base;calib;waypoints_utm;base];
    
% generate array of speeds
    all_spd = [home_speed;calib_speed*ones(size(calib,1),1);...
        speed*ones(size(waypoints_utm,1),1);...
        home_speed];
    
% calc approx. flight time
    takeoff_land = 30; % [s]
    
    % base to calibration
    D = pdist2(base,begin);
    to_area_time = D/home_speed;
    
    % calib pattern time
    calib_time = 0;
    for ii = 2:size(calib,1)
        D = pdist2(calib(ii,:),calib(ii-1,:));
        calib_time = calib_time + D/calib_speed;
    end
    
    % scan pattern time
    scan_time = 0;
    for ii = 2:size(waypoints_utm,1)
        D = pdist2(waypoints_utm(ii,:),waypoints_utm(ii-1,:));
        scan_time = scan_time + D/speed;
    end
    
    % end to base
    D = pdist2(waypoints_utm(end,:),begin);
    to_base_time = D/home_speed;
    
    % all time
    mission_time = 2*takeoff_land + to_area_time + calib_time + ...
        scan_time + to_base_time;
    fprintf('\n   // Mission time = %.1f min\n',round(mission_time/60,1));

% convert to lon/lat
    waypoints = utm2ll(all_pts);
  
% round lon/lat waypoints to nearest 1e-6 (roughly 0.1-m precision)
    waypoints = round(rad2deg(waypoints),6);
    
% write KML and TXT file
    writeKML(waypoints,file);
    writeTXT(round(mission_time/60,1),waypoints,all_spd,h,file);
    writeAWM(waypoints,all_spd,h,file);
    open([file '.txt']) % check
end
    