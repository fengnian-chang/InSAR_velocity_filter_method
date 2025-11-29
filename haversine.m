function distance_km = haversine(lat1, lon1, lat2, lon2)
    % 
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % 
    R = 6371;

    % Haversine
    dlon = lon2 - lon1;
    dlat = lat2 - lat1;
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    distance_km = R * c;
end
