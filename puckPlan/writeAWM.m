function writeAWM(wp, v, hgt, infile)
% write mission to DJI GS AWM file

    fid = fopen([infile '.awm'], 'w');
    
    fprintf(fid, '<?xml version="1.0" encoding="ISO-8859-1" standalone="yes"?>\n');
    fprintf(fid, '<Mission MissionTimeLmt="65535" IsPatrol="Continuous" ');
    fprintf(fid, 'StartWayPointIndex="0" VerticalSpeedLimit="2">\n');

    for ii = 1:size(wp, 1)
        fprintf(fid, '  <WayPoint id="%i">\n', ii - 1);
        fprintf(fid, '    <Latitude>%.6f</Latitude>\n', wp(ii, 2));
        fprintf(fid, '    <Longitude>%.6f</Longitude>\n', wp(ii, 1));
        fprintf(fid, '    <Altitude>%.0f</Altitude>\n', hgt);
        fprintf(fid, '    <Speed>%.0f</Speed>\n', v(ii));
        fprintf(fid, '    <TimeLimit>1000</TimeLimit>\n');
        fprintf(fid, '    <YawDegree>360</YawDegree>\n');
        fprintf(fid, '    <HoldTime>3</HoldTime>\n');
        fprintf(fid, '    <StartDelay>0</StartDelay>\n');
        fprintf(fid, '    <Period>0</Period>\n');
        fprintf(fid, '    <RepeatTime>0</RepeatTime>\n');
        fprintf(fid, '    <RepeatDistance>0</RepeatDistance>\n');
        fprintf(fid, '    <TurnMode>Adaptive_Bank_Turn</TurnMode>\n');
        fprintf(fid, '  </WayPoint>\n');
    end

    fprintf(fid,'</Mission>\n');
    
    fclose(fid);
end