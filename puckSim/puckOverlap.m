function puckOverlap(infile, width)
%PUCKOVERLAP(infile, width) produces a point cloud that simulates three
%   overlapping passes by a UAS-mounted Velodyne PUCK(TM) VLP-16 Laser Scanner. 
%   The input must be the output from PUCKSIM, which is the result of simulating
%   a single pass by the UAS-mounted PUCK. The two passes will be parallel, 
%   opposite in direction, and separated by width W.
%   
%   The name of the output file from the last run of PUCKSIM in the current
%   session will be saved in the base workspace as 'simfile'.
%
%   Example:
%      % Simulate overlapping passes whose flight lines are parallel, opposite
%      % in direction, and separated by 50 m. The name of the output file from 
%      % PUCKSIM is in the workspace as 'simfile'.
%      
%      puckOverlap(simfile, 50)
%      
%      % Output saved to overlap_*.txt, where * = simfile (sans file ext.).
%
%   See also: puckSim, puckStats.

    % read infile, skipping header row
    in = dlmread(infile, '', 1, 0);
    in = in';  % pass 2: these will be rotated 180� and translated by width
    in_ = in;  % pass 3: these will be translated by 2*width

    % pass 2
    % 2D conformal transformation on points [X,Y] and [Tx,Ty]
    rot = [-1 0; 0 -1];

    % translation in Y to center points at origin
    T = (max(in(2, :)) + min(in(2, :))) / 2;
    in(1:2, :) = in(1:2, :) + [0; -T];
    in(7:8, :) = in(7:8, :) + [0; -T];

    % rotate and translate in X
    in(1:2, :) = rot * in(1:2, :) + [width; 0];
    in(7:8, :) = rot * in(7:8, :) + [width; 0];

    % undo translation in Y
    in(1:2, :) = in(1:2, :) + [0; T];
    in(7:8, :) = in(7:8, :) + [0; T];

    % add arbitrary constant to time
    in(6, :) = in(6, :) + 30;

    % pass 3: Tx and time
    in_(1:2, :) = in_(1:2, :) + [2 * width; 0];
    in_(7:8, :) = in_(7:8, :) + [2 * width; 0];
    in_(6, :) = in_(6, :) + 60;

    results = [in; in_];

    % output new text file
    outfile = ['overlap_' num2str(width) '_' infile];

    % start with existing file, append results
    copyfile(infile, outfile);
    fid = fopen(outfile, 'a');
    str = [
        '%8.3f %8.3f %8.3f %11.8f %11.8f %12.9f' ...
        '%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n'
    ];
    fprintf(fid, str, results);
    fclose(fid);

    % assign name in workspace to save some time
    assignin('base', 'overlapfile', outfile)

end
