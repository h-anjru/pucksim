function puckSim(height, speed, varargin)
%PUCKSIM(height, speed) simulates the point cloud of a UAS-mounted Velodyne PUCK
%   VLP-16 Laser Scanner.
%   
%   This function models the VLP-16 scanner as a point moving parallel to a
%   target plane (specified in puckInt.m) at the specified flying height. The 
%   laser returns are modeled by lines with a known direction passing through 
%   the scanner point and intersecting the target plane.
%   
%   The scanner is assumed to be oriented such that its +z axis is coincident
%   with the direction of travel, which has been arbitrarily chosen to be the +Y
%   direction in object space. The scanner's y-axis is coincident with the -Z
%   direction on object space. This can be manipulated as described below.
%   
%   PUCKSIM(height, speed) outputs a text file that records the intersections of 
%   the modeled laser returns with the target plane. These inputs are required. 
%   A speed of 0 (zero) will model a single rotation of the scanner head.
%   
%   PUCKSIM(height, speed, 'rotationRate', R) changes the rotation rate of the
%   scanner head. By default, the rotation rate is set to 10 Hz. The VLP-16 
%   manual indicates the rotation rate can vary between 5-20 Hz (p.20).
%   
%   PUCKSIM(height, speed, 'profile', L) outputs a profile that shows the area 
%   of the target plane that was "fully illuninated," i.e. an area that was 
%   struck by lasers from all 16 channels. This profile is a "strip" of points 
%   L-m long in the direction of flight. This profile can be used to represent 
%   theoretical point density of a strip for a flight over an unobstructed 
%   target plane.
%   
%   PUCKSIM(height, speed, 'maxRange', M) limits the maximum range of recorded 
%   returns to M meters. The VLP-16 records returns up to 120 m, the default 
%   value for maxRange.
%
%   PUCKSIM(height, speed, 'tilt', T) changes the orientation of the scanner. 
%   By default, the scanner is oriented as described above; this can also be 
%   thought of as the scanner "facing down." (More formally, this is a rotation 
%   of -90° about the scanner's x-axis.)
%  
%   PUCKSIM(height, speed, 'yaw', Y) changes the orientation of the scanner
%   about the airframe's +Z-axis. This is useful for simulating crab.
%   
%   Example:
%      % Get a 5-m profile of points collected from a flying height of 50 m, 
%      % forward speed of 10 m/s, and a rotation rate of 20 Hz, with the scanner
%      % inclined 15° from nadir (-75° about its x-axis), and the aircraft
%      % crabbing at -5°. Limit range to 60 m.
%      
%      puckSim(50, 10, 'rotationRate', 20, 'profile', 5, 'maxRange', 60, ...
%         'tilt', 15, 'yaw', -5)
%      
%      % Output saved to results_h50_s10_r20_m60_t15_y-5_p5.txt.
%   
%   Reference: VLP-16 User's Manual and Programming Guide.
%
%   See also: puckOverlap, puckStats.

    % input parser--defaults, type enforcing
    p = inputParser;
    defaultProfile = false;
    defaultRotationRate = 10;
    defaultMaxRange = 120;
    defaultTilt = 0;
    defaultYaw = 0;

    addRequired(p, 'height', @isnumeric);
    addRequired(p, 'speed', @isnumeric);
    addParameter(p, 'profile', defaultProfile, @isnumeric)
    addParameter(p, 'rotationRate', defaultRotationRate, @isnumeric);
    addParameter(p, 'maxRange', defaultMaxRange, @isnumeric);
    addParameter(p, 'tilt', defaultTilt, @isnumeric);
    addParameter(p, 'yaw', defaultYaw, @isnumeric);

    parse(p, height, speed, varargin{:});

    % change tilt to be in terms of rotation about scanner's x-axis [rad]
    omega_s = deg2rad(p.Results.tilt) - pi/2;

    % change yaw to be in terms of rotation about airframe's z-axis [rad]
    kappa_b = -deg2rad(p.Results.yaw);

    % vertical angles of the 16 channels (invariant); Manual, p.21
    % not to be confused with omega_s
    omega = [-15 1 -13 3 -11 5 -9 7 -7 9 -5 11 -3 13 -1 15]; % [°]
    omega = deg2rad(omega);

    % timing information (Manual, p.16)
    % a suspected sig fig issue in manual; overlooking that for now
    firing_interval = 18.43E-6;
    timing_cycle = 55.296E-6;

    % initial vector of epochs for single firing group
    t_cycle = 0:firing_interval:15 * firing_interval;

    % range of azimuths at which a laser could theoretically intersect with the
    % target plane is (3pi/2, pi/2). The PUCK has a maximum range of 100 meters,
    % so this range of azimuiths is limited further.
    range_max = p.Results.maxRange;  % [m]
    a_max = acos(p.Results.height / range_max);  % a_max at omega = 0
    a_range = mod(1.02 * [-a_max a_max], 2 * pi);  % 2% buffer
    a0 = a_range(1);  % initial value for a0

    % initialize azimuth as a function of time, range [0-2pi)
    azi = mod(a0 + 2 * pi * p.Results.rotationRate * t_cycle, 2 * pi);

    % locations [X,Y,Z] of scanner traveling in +Y-dir. during t_cycle
    XYZ = [
        zeros(1,16); p.Results.speed * t_cycle; ...
        ones(1,16) * p.Results.height
    ];

    % solve for the minimum time of flight needed to "fully illuminate" a 
    % section of the target plane, and then calculate the bounds of the "fully
    % illuminated" profile area. Tripped if 'profile' value is entered as an
    % input argument.
    if p.Results.profile ~= false
        % find minimum distance scanner must travel
        a_max_ = acos(
            (p.Results.height / range_max + cos(omega_s) * sin(max(omega))) / ...
            (-sin(omega_s) * cos(max(omega)))
        );

        must_reach = puckInt(XYZ(:,1), a_max_, max(omega), omega_s, abs(kappa_b));

        % how far must scanner travel = XYZ2
        vec = must_reach - XYZ(:,1);
        vec2 = -vec; vec2(2) = vec(2);
        XYZ2 = must_reach + vec2;

        % calculate bounds of profile
        % x_bounds = [-must_reach(1) must_reach(1)]; 
        y_bounds = [must_reach(2) must_reach(2) + p.Results.profile];

        % find flight time needed to cover the profile
        travel = XYZ2(2) - XYZ(1);
        flight_time = 1.1 * travel / p.Results.speed; % 10% buffer

        % name of file
        name = [
            'results_h' num2str(p.Results.height) ...
            '_s' num2str(p.Results.speed) ...
            '_r' num2str(p.Results.rotationRate) ...
            '_m' num2str(p.Results.maxRange) ...
            '_t' num2str(p.Results.tilt) ...
            '_y' num2str(p.Results.yaw) ...
            '_p' num2str(p.Results.profile) '.txt'
        ];

    else
        travel = 10*p.Results.speed; % arbitrary value
        % exception: if speed = 0, only need one rotation
        if p.Results.speed == 0
            flight_time = 1 / p.Results.rotationRate;
        else
            flight_time = travel / p.Results.speed; % [s]
        end
        name = ['results_h' num2str(p.Results.height) ...
            '_s' num2str(p.Results.speed) ...
            '_r' num2str(p.Results.rotationRate) ...
            '_m' num2str(p.Results.maxRange) ...
            '_t' num2str(p.Results.tilt) ...
            '_y' num2str(p.Results.yaw) '.txt'];
    end

    % open text file for recording results
    fprintf([name '\n']); % print file name in command window
    fid = fopen(name,'wt'); % open file for writing

    % write header row of text file
    fprintf(fid,[' X Y Z azimuth vert_angle time ' ...
        'Tx Ty Tz range dirx diry dirz\n']);

    % format string for text file; variable called in while loop
    str = [
        '%8.3f %8.3f %8.3f %11.8f %11.8f %12.9f' ...
        '%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n'
    ];

    % counter to trip progress display
    jj = 0;

    while t_cycle(1) < flight_time
        % progress display
        if mod(jj, 7768) == 0
            fprintf('%6.1f%% complete...\n',100*t_cycle(1) / flight_time)
        end
        jj = jj + 1; % increment counter

        % only solve for intersections if azimuth is in range (saves time)
        if azi(1) > a_range(1) || azi(16) < a_range(2)
            % solve for line-plane intersections
            [results, direction] = puckInt(XYZ, azi, omega, omega_s, kappa_b);
            range = puckRange(XYZ, results);

            % file info; each column is a modeled return
            file_info = [results; azi; omega; t_cycle; XYZ; range; direction];

            % look for ranges > 100 and mark for deletion
            columns_to_remove = file_info(10, :) > range_max;  % result: bool
            file_info(:, columns_to_remove) = [];  % remove columns range > 100

            % look for points not in profile and mark for deletion
            if p.Results.profile ~= false
                columns_to_remove2 = ...
                file_info(2,:) < y_bounds(1) | ...
                file_info(2,:) > y_bounds(2);
                % file_info(1,:) < x_bounds(1) | ...
                % file_info(1,:) > x_bounds(2) | ...
                file_info(:, columns_to_remove2) = [];
            end

            % write results to text file (fprintf will write in rows)
            if isempty(file_info)
                % do nothing
            else
                fprintf(fid, str, file_info);
            end
        end

        % update epochs for next firing group
        t_cycle = t_cycle + timing_cycle;
        % update azimuiths for next firing group
        azi = mod(a0 + 2 * pi * p.Results.rotationRate * t_cycle, 2 * pi);
        % update location
        XYZ(2, :) = p.Results.speed * t_cycle;
    end

    fprintf(' 100.0%% complete.\n') %  final progress update
    fclose(fid);  % close text file

    % assign name in workspace to prepare to use puckStats.m
    assignin('base', 'simfile', name)
end


function ranges = puckRange(scanner_point, points)
% PUCKRANGE Solves for the slant ranges between the modeled scanner positions 
%   and the modeled laser returns for the puckPointDensity function.
%   
%             Call: ranges = puckRange(scanner_point, points)
%   
%    scanner_point: 3 x n coordinates of points on line
%           points: 3 x n coordinates of line-plane intersection
%   
%           rangess: 1 x n vector of slant distances

    diff_squared = (scanner_point - points) .^ 2;
    ranges = sqrt(sum(diff_squared, 1));

end


function [points,varargout] = puckInt(scanner_point,azimuth,channel,os,kb)
% PUCKINT Solves for the point of intersection of a line and plane in 3-space
%   for the puckPointDensity function.
%  
%             Call: point = puckInt(plane_normal, plane_point, ...
%                      line_direction, line_point)
%    
%    scanner_point: 3 x n coordinates of points on line
%          azimuth: 1 x n azimuths (alpha) of lines in SOCS [rad]
%          channel: 1 x n vertical angles (omega) of channel in SOCS [rad]
%               os: rotation of scanner about its x-axis (i.e. omega_s)
%               kb: rotation of airframe about its z-axis (i.e. kappa_b)
%
%           points: 3 x n coordinates of line-plane intersection
%
%   Reference:
%      en.wikipedia.org/wiki/Line%E2%80%93plane_intersection#Algebraic_form

    n = size(scanner_point, 2);

    % define target plane
    plane_normal = [0 0 1]'; plane_normal = repmat(plane_normal, 1, n);
    plane_point = [0 0 0]'; plane_point = repmat(plane_point, 1, n);

    % transform azimuth and channel angles to Cartesian vector in SOCS
    x = cos(channel) .* sin(azimuth);
    y = cos(channel) .* cos(azimuth);
    z = sin(channel);
    line_direction_SOCS = [x; y; z];

    % rotation between SOCS and object space
    Ro = makehgtform('xrotate', os);
    Rk = makehgtform('zrotate', kb);
    R = Rk * Ro;
    R = R(1:3, 1:3);
    line_direction = R * line_direction_SOCS;

    % initialize array of results
    points = zeros(3, n); % points
    varargout = points;  % line directions (below)

    % solve for line-plane intersections
    for ii = 1:n
        d = dot(
            (plane_point(:, ii) - scanner_point(:, ii)),
            plane_normal(:, ii)
        ) / (dot(line_direction(:, ii), plane_normal(:, ii)));

        points(:, ii) = d * line_direction(:, ii) + scanner_point(:, ii);
    end

    % return line_direction info (this looks weird...)
    nout = max(nargout, 1) - 1;
    varargout = cell(nout, 1);
    for kk = 1:nout
        varargout{kk} = line_direction;
    end
end
