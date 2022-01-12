% run puckStats
puckStats(simfile)

% input file for gap plot script below
infile = simfile;

% find height, speed, yawpoint density eq plot
temp_cell = extractBetween(infile, '_h', '_');
ht = temp_cell{1}; ht = str2num(ht);
temp_cell = extractBetween(infile, '_s', '_');
sp = temp_cell{1}; sp = str2num(sp);
temp_cell = extractBetween(infile, '_r', '_');
ro = temp_cell{1}; ro = str2num(ro);
temp_cell = extractBetween(infile, '_y', '_');
ya = temp_cell{1}; ya = str2num(ya);

x = gap_equation(ht, sp, ro);

pts = dlmread(infile, '', [1 0 0 1]);  % XY of returns

% plot(x, zeros(size(x)), 'b*', 'MarkerSize', 10, 'LineWidth',2)
% plot(-x, zeros(size(x)), 'b*', 'MarkerSize', 10, 'LineWidth',2)

if ya ~= 0
    x_ = x * cosd(ya);
    plot(x_, zeros(size(x)), 'ro', 'MarkerSize', 15, 'LineWidth', 2)
    plot(-x_, zeros(size(x)), 'ro', 'MarkerSize', 15, 'LineWidth', 2)
end

% save figure
saveas(f1, [simfile '_gaps' '.png'])

f2 = figure(2); clf(2);
a2 = axes('parent', f2);
hold on;
f2.Position = f1.Position;
h = plot(pts(:, 1), pts(:, 2), 'b.', 'MarkerSize', 3);
h.MarkerEdgeColor = [0 0.4470 0.7410];  % a kinder shade of blue

% lTick = round(min(pts(:,1)),0); rTick = round(max(pts(:,1)),0);
lTick = -100; rTick = 100;
a2.XTick = linspace(lTick, rTick, 11);
% a2.XLim = [min(pts(:,1)) max(pts(:,1))];
a2.XLim = [-100 100];
a2.YLim = [min(pts(:,2)) max(pts(:,2))];
a2.XAxis.FontSize = 14; a2.YAxis.FontSize = 14;

% replace hyphen with en dash (minus sign)
a2.XTickLabels = strrep(a2.XTickLabels, '-', '−');

title(infile,'Interpreter', 'none')
xlabel('across track [m]', 'FontSize', 14); 
ylabel('along track [m]', 'FontSize', 14); 

% name output
string = infile(1:end - 4);
saveas(f2, [string '_points.png'])

function x_gaps = gap_equation(height, speed, rotation)
        dw = 2;  % angular separation of channels [deg]
    
        ii = 1:10; ii = ii'; % vector for finding "iith" gap
        
        % gap equation (no yaw)
        xi = height*tan(acos((rotation * height * tand(dw)) ./ (ii * speed)));
    
        x_gaps = real(xi);
        x_gaps = x_gaps(x_gaps < 110);
    end