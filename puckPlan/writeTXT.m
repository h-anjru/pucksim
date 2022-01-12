function writeTXT(t,wp,v,hgt,infile)
% write mission to human-readable text file
    
    fid = fopen([infile '.txt'],'w');
    
    fprintf(fid,'%s\n-------\n',infile);
    fprintf(fid,'Flying height = %.0f m\n',hgt);
    fprintf(fid,'Mission time = %.0f min\n\n',t);
    fprintf(fid,'WAYPOINTS [WGS84]:\n');
    fprintf(fid,'    LON         LAT       SPEED\n');
    for ii = 1:size(wp,1)
        fprintf(fid,'%10.6f %10.6f %5.0f m/s\n',...
            wp(ii,1),wp(ii,2),v(ii));
    end
    
    fclose(fid);
end