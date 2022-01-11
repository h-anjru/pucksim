function writeKML(wp,infile)
% write waypoints to KML

    fid = fopen([infile '.kml'],'w');

    fprintf(fid,'<?xml version="1.0" encoding="ISO-8859-1"?>\n');
    fprintf(fid,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
    fprintf(fid,'<Document id="root_doc">\n');
    fprintf(fid,'<Folder><name>%s</name>\n',infile);
    fprintf(fid,'  <Placemark>\n');
    fprintf(fid,'  <Style><LineStyle><color>aaa0dd00</color>');
    fprintf(fid,'<width>3.0</width></LineStyle>');
    fprintf(fid,'<PolyStyle><fill>0</fill></PolyStyle></Style>\n');
    fprintf(fid,'    <LineString><coordinates>\n');
    for ii = 1:length(wp)
        fprintf(fid,'      %.6f,%.6f \n',wp(ii,1),wp(ii,2));
    end
    fprintf(fid,'    </coordinates></LineString>\n');
    fprintf(fid,'  </Placemark>\n');
    fprintf(fid,'</Folder>\n');
    fprintf(fid,'</Document></kml>\n');

    fclose(fid); % close text file
end