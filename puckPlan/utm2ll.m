function lonlat = utm2ll(Xin)
% convert waypoints from UTM to lat/lon

    % WGS84 ellipsoid parameters
    a = 6378137;              % semimajor axis
    b = 6356752.314245;       % semiminor axis
    f = 1 / 298.257223563;    % flattening
    esq = f * (2 - f);        % first eccentricity squared
    eee = sqrt(esq);          % first eccentricity
    epsq = esq / (1 - f)^2;   % second eccentricity squared

    % UTM zone 17N parameters
    ku = 0.9996;         % scale factor at central meridian
    fu = 0;              % latitude of origin
    lu = degtorad(-81);  % central meridian for UTM Zone 17N
    eu = 500000;         % false easting 17N
    nu = 0;              % false northing 17N

    mu = a * ((1 - eee^2 / 4 - 3 * eee^4 / 64 - 5 * eee^6 / 256) * fu ...
        - (3 * eee^2 / 8 + 3 * eee^4 / 32 + 45 * eee^6 / 1024) * ... 
        sin(2 * fu) + (15 * eee^4 / 256 + 45 * eee^6 / 1024) * ... 
        sin(4 * fu) - (35 * eee^6 / 3072) * sin(6 * fu));

    e1 = (1 - sqrt(1 - esq)) / (1 + sqrt(1 - esq));
    
    % UTM to lat/lon
    M = mu + (Xin(:,2) - nu) ./ ku;
    MU = M ./ (a * (1 - eee^2 / 4 - 3 * eee^4 / 64 - 5 * eee^6 / 256));
    F1 = MU + (3 * e1 /2 - 27 * e1^3 / 32) * sin(2 .* MU) + ... 
        (21 * e1^2 / 16 - 55 * e1^4 / 32) * sin(4 .* MU) + ...
        (151 * e1^3 / 96) * sin(6 .* MU) + (1097 * e1^4 / 512) * ...
        sin(8 .* MU);
    R1 = a * (1 - esq) ./ (1 - esq .* sin(F1).^2).^1.5;
    C1 = epsq .* cos(F1).^2;
    T1 = tan(F1).^2;
    N1 = a ./ sqrt(1 - esq .* sin(F1).^2);
    D = (Xin(:,1) - eu) ./ (N1 .* ku);

    lat = F1 - (N1 .* tan(F1) ./ R1) .* (D.^2 ./ 2 - (5 + 3 .* T1 + ...
        10 .* C1 - 4 .* C1.^2 - 9 * epsq) .* (D.^4 ./ 24) + (61 + 90 .* ...
        T1 + 298 .* C1 + 45 .* T1.^2 - 252 * epsq - 3 * C1.^2) .* ...
        (D.^6 ./ 720));  % latitude output
    lon = lu + (1 ./ cos(F1)) .* (D - (1 + 2 .* T1 + C1) .* ...
        (D.^3 ./ 6) + (5 - 2 .* C1 + 28 .* T1 - 3 .* C1.^2 + 8 * epsq + ...
        24 .* T1.^2) .* (D.^5 ./ 120)); % longitude output
  
    lonlat = [lon lat];
end