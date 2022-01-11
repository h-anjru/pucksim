function puckStats(infile)
% PUCKSTATS calculates spatial statistics on a point cloud generated by 
%   PUCKSIM or PUCKOVERLAP.
%   
%   It is recommended to use this function on a point cloud "profile" (for
%   more info on profile, type "help puckSim"). This will eliminate edge
%   effects and reduce runtime significantly. The input must be the output 
%   from PUCKSIM or PUCKOVERLAP or identically formatted.
%   
%   The name of the output file from the last run of PUCKSIM in the current
%   session will be saved in the base workspace as 'simfile'. The name of
%   the output file from the last run of PUCKOVERLAP in the current session
%   will be saved in the base workspace as 'overlapfile'. You can use
%   either of those variable names as input for PUCKSIM.
%   
%   WARNING: For each bin, a square distance matrix must be computed in
%   order to run the spatial analysis. Executing PUCKSTATS with too few
%   bins could lead to memory issues. At least 80 bins are recommended for
%   analyzing an overlap file.
%   
%   Reference available <a
%   href="http://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/h-how-average-nearest-neighbor-distance-spatial-st.htm">here</a>.
%   
%   See also: puckSim, puckOverlap.

% read first two columns of text file, skipping header row
in1 = dlmread(infile,'',[1 0 0 1]); % XY of returns
in2 = dlmread(infile,'',[1 9 0 9]); % range
in3 = dlmread(infile,'',[1 12 0 12]); % z of vector

pts = [
    in1 ...
    in2 ...
    in3 ...
    ];

% shift (good for large files and overlap sims)
dX = max(pts(:,1)) - 0.5*range(pts(:,1)); dY = mean(pts(:,2));
pts(:,1:2) = pts(:,1:2) - [dX dY];

% one bin per meter across-track
num_bins = floor(range(pts(:, 1)));

clear in1 in2 in3

% set up bins
bin_width = range(pts(:,1)) / num_bins;
bins = min(pts(:,1)) + bin_width*(1:num_bins);
bins = [-9999 bins]; % add to beginning to facilitate logical test below

% cell array to store matrices of points for each bin
bin_cells = cell(1,num_bins);

% add points to bins
for ii = 1:num_bins
    bin_holder = [];
    for jj = 1:length(pts)
        if pts(jj,1) > bins(ii) && pts(jj,1) < bins(ii + 1)
            bin_holder = [bin_holder;pts(jj,:)];
        end
        bin_cells{ii} = bin_holder;
    end
end

% individual plots for bins
% for ii = 1:num_bins
%     figure(ii+1)
%     C = bin_cells{ii}; % store cell as matrix
%     plot(C(:,1),C(:,2),'.');
%     title('bin '+string(bin_centers(ii)))
%     grid on;xlabel('X');ylabel('Y');
% end

% initialize some matrices to store results
znn = zeros(num_bins,1); ppm = znn; avg = ppm; sdv = avg; rng = sdv;

% average nearest neighbor + statistics for each bin
for kk = 1:num_bins
    B = bin_cells{kk};             % store cell as matrix
    C = B(:,1:2);                  % only keep X and Y
    R = B(:,3);                    % ranges
    Z = B(:,4);                    % z-components of vectors
    D = squareform(pdist(C));      % square distance matrix
    D(eye(size(D))~=0)=inf;        % replace zeroes on diagonal with inf
    min_distance = min(D);
    clear D
    area = range(C(:,1))*range(C(:,2)); % is this right?
    n = length(C);                 % number of points in bin
    d_exp = 0.5/sqrt(n/area);      % expected average NN distance
    SE = 0.26136/sqrt(n^2/area);   % standard error
    d_obs = mean(min_distance);    % observed average NN distance
    avg(kk) = d_obs;
    nn_diff = d_obs - d_exp; 
    znn(kk) = nn_diff / SE;        % standardized NN distance
    ppm(kk) = n / area;            % points per unit area [m^2]
    sdv(kk) = std(min_distance);
    rng(kk) = mean(R);
    ang(kk) = mean(acosd(abs(Z))); % angle from nadir (unsigned)
end

fprintf('%s\n',infile(1:end-4))

puckStatPlot(infile(1:end-4), bins, bin_width, num_bins, znn, ppm, ...
    'z-score', bone, pts)
% puckStatPlot(infile(1:end-4),bins,bin_width,num_bins,rng,ppm,...
%     'range',jet,pts)
% puckStatPlot(infile(1:end-4), bins, bin_width, num_bins, ang, ppm, ...
%     'scan angle', jet, pts)

% % p-values for average nearest neighbor method
% p = 2*(1-normcdf(abs(znn),0,1));
% 
% % display results in command window & assign in worksapce as 'results'
% fprintf('\n       avg       znn       ppm     p-val     stdev\n')
% results = [avg znn ppm p sdv];
% assignin('base','results',results)

end

function puckStatPlot(string, bin_val, bin_w, bin_n, x_val, y_val, ...
    name, map, pts)
    % figure handle, axes handle, title, etc
    f1 = figure(1); clf(f1);
    a1 = axes('parent', f1);
    title(string, 'Interpreter', 'none')
    grid on; hold all;
    f1.Position = [300 300 1200 350];
    
    % xlabel, ylabel (use directly)
    xlabel('X_{bin} [m]', 'FontSize', 18); 
    ylabel('pts/m^2', 'FontSize', 18); 

    % convert range to range of colormap
    colors = colormap(map);
    s_ = x_val - min(x_val) + 1;
    s_ = s_ * (length(colors) / max(s_));
    
    % find bin centers
    bin_centers = round(bin_val(2:end) - 0.5 * bin_w, 1);

    % plot each bar in histogram separately
    for mm = 1:length(y_val)
        H = bar(bin_centers(mm), y_val(mm), 'parent', a1, ...
            'facecolor', colors(ceil(s_(mm)), :));
        H.BarWidth = range(bin_centers) / bin_n;
    end

    % annotation
    % anno = ['width = ' num2str(width) ' m'];
    % annotation(f1,'textbox',[0.8 0.6 0.1 0.1],'String',anno,'FontSize',14)

    % colorbar properties
    c1 = colorbar('peer', a1); 
    c1.Label.String = name;
    c1.Ticks = 0:0.1:1;
    jj = round(linspace(min(x_val), max(x_val), 11), 1);
    c1.TickLabels = {
        num2str(jj(1)), num2str(jj(2)), ...
        num2str(jj(3)), num2str(jj(4)), ...
        num2str(jj(5)), num2str(jj(6)), ...
        num2str(jj(7)), num2str(jj(8)), ...
        num2str(jj(9)), num2str(jj(10)), ...
        num2str(jj(11))
    };
    % replace hyphen with minus sign
    c1.TickLabels = strrep(c1.TickLabels, '-', '−');
    c1.Location = 'northoutside'; c1.FontSize = 12;

    % tick label properties
    lTick = round(min(pts(:,1)),0);
    rTick = round(max(pts(:,1)),0);
    
%     a1.XTick = linspace(lTick, rTick, 9);
%     a1.XLim = [min(bin_centers)-H.BarWidth max(bin_centers)+H.BarWidth];
    
    a1.XTick = linspace(-70, 70, 15);
    a1.XLim = [-70 70];
    
    % replace hyphen with en dash (minus sign)
    a1.XTickLabels = strrep(a1.XTickLabels, '-', '−');
    
    a1.YLim = [0 300];
    
    a1.XAxis.FontSize = 14;
    a1.YAxis.FontSize = 14;
    
    % find height and speed for point density eq plot
    temp_cell = extractBetween(string,'_h','_');
    ht = temp_cell{1}; ht = str2num(ht);
    temp_cell = extractBetween(string,'_s','_');
    sp = temp_cell{1}; sp = str2num(sp);
    temp_cell = extractBetween(string,'_y','_');
    ya = temp_cell{1}; ya = str2num(ya);
    
    % plot point density
    lf = 3e5; x = linspace(a1.XLim(1), a1.XLim(2), 200);
    px = (lf .* ht) ./ (2 .* sp .* pi .* (ht^2 + x.^2));
    px_ = (lf * ht * cosd(ya)) ./ (2 .* sp .* pi .* (x.^2 + ht^2 * cosd(ya)^2));
    plot(x, px, 'k--', 'LineWidth', 4)
    plot(x, px_, 'k', 'LineWidth', 4)

    % assign variables (debugging)
    saveas(f1,[string '_' name '.png'])
    assignin('base','f1',f1)
    assignin('base','a1',a1)
    assignin('base','c1',c1)
end

%   o-------------------------------------------------o   %
%   |    H. Andrew Lassiter (halassiter@ufl.edu)      |   %
%   |        Created: 01 March 2017                   |   %
%   |       Modified: 17 October 2017                 |   %
%   o-------------------------------------------------o   %