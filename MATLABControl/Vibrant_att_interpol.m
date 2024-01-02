% Function description: 
% Calculates the PA attenuator settings 
%
% Created by Krista Jansen, May 12, 2010.

function PA_wavelengtharray = Vibrant_att_interpol(xi)
x = [600 700 710 750 800 850 900 1000 1100 1200 1300 1400 1500 1600 1700];
y = [];
Y=[];
yi = interp1(x,Y,xi,'linear'); % linear interpolation
PA_wavelengtharray=[xi,yi];
