function dva = pix2deg(pix,scrd,pixpercm)
%PIX2DEG takes a number of pixels with the screen distance (in cm)
% and the pixels per cm and returns the degrees of visual angle
%
% pix2deg(pix,scrd,pixpercm)

% Matthew A. Smith
% Revised: 20130805

d = pix./pixpercm;
angle = atan(d./scrd);

dva = (180/pi) * angle;