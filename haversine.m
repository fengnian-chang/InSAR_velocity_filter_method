function distance_km = haversine(lat1, lon1, lat2, lon2)
%HAVERSINE  Great-circle distance (km) between two points given in degrees.
%
%   distance_km = haversine(lat1, lon1, lat2, lon2)
%
%   Uses the haversine formula on a sphere with radius 6371 km. Used to
%   convert velocity differences between neighboring pixels to gradients.

    % Degrees to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % Mean Earth radius (km)
    R = 6371;

    % Haversine formula
    dlon = lon2 - lon1;
    dlat = lat2 - lat1;
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distance_km = R * c;
end
