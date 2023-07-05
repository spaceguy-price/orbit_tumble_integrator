function [r_lla] = ecef2lla(r_ecef)
%ecef2lla Converts ECEF coordinates to LLA coordinates
%   Converts Earth-centered, Earth-fixed coordinates to the Latitude
%   Longitude Altitude.
%   
%   See https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_geodetic_coordinates
%   for an overview.
%   See https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=303772
%   for the mathmetics
%   (Conversion of Earth-Centered Earth-Fixed Coordinates to Geodetic 
%   Coordinates) IEEE 1994 publication
%
%   Inputs:
%       -r_ecef [x,y,z] [Nx3]
%   Outputs:
%       -r_lla [long, lat, alt] [Nx3]

% Andrew Price
% 14 December 2021

x = r_ecef(:,1);
y = r_ecef(:,2);
z = r_ecef(:,3);

%% Ferrari's solution
a = 6378137.0;      % [m] Earth Equatorial Radius
b = 6356752.3142;   % [m] Earth Polar Radius
e2  = (a^2 - b^2)/(a^2);
er2 = (a^2 - b^2)/(b^2);
p = sqrt(x.^2 + y.^2);
F = 54*b^2.*(z.^2);
G = p.^2 + (1-e2)*(z.^2) - e2*(a^2 - b^2);
c = ((e2^2)*F.*(p.^2)) ./ (G.^3);
s = (1 + c + (c.^2 + 2*c).^(1/2)).^(1/3);
k = s + 1 + 1./s;
P = F ./ (3*(k.^2)*(G.^2));
Q = (1 + 2*(e2^2).*P).^(1/2);
r_o = (-P*e2.*p)./(1 + Q) + ((0.5*(a^2)*(1 + 1./Q)) - ((P.*(1-e2).*(z.^2))./(Q.*(1+Q))) - (0.5*P.*(p.^2))).^(1/2);
U = ((p - e2.*r_o).^2 + z.^2).^(1/2);
V = ((p - e2.*r_o).^2 + (1-e2)*(z.^2)).^(1/2);
z_o = ((b^2)*z)./(a*V);
h = U.*(1 - (b^2)./(a*V));
lat = atand((z + er2.*z_o)./p);
long = atan2d(y,x);

r_lla = [lat, long, h];
end

