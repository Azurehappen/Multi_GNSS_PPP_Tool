function spec_ecef = specInEnuToEcef(pos_lla, spec)

% pos_lla: position in Lat, Lng, H (rad, rad, m)
% spec: accuracy specification in ENU (3 x 3)

lat = pos_lla(1);
lng = pos_lla(2);

R = [-sin(lng), -sin(lat)*cos(lng), cos(lat)*cos(lng);
     cos(lng), -sin(lat)*sin(lng), cos(lat)*sin(lng);
     0, cos(lat), sin(lat)];
spec_ecef = R * spec * R';