function [alpha, cl, cd] = viterna(alpha, cl, cd)
% Apply Viterna extrapolation to the aerofoil given the minimal polar data
% (alpha, cl, and cd)
% Usage: [alpha, cl, cd] = VITERNA(alpha, cl, cd);
% As everywhere else, alpha is in RADIANS

% define values at stall (as where cl is max)
[cl_st, i_st] = max(cl);
a_st = alpha(i_st);
cd_st = cd(i_st);

% define parameters for extrapolation
ar = 10;
cd_max = 1.11 + 0.018 * ar;
a1 = cd_max/2;
b1 = cd_max;
a2 = (cl_st - cd_max * sin(a_st) * cos(a_st)) * sin(a_st)/cos(a_st)^2;
b2 = cd_st - cd_max * sin(a_st)^2/cos(a_st);

% values from viterna
da = mean(diff(alpha));
a_vit = (-5/6*pi:da:5/6*pi).';
cl_vit = a1 * sin(2*a_vit) + a2 * cos(a_vit).^2./sin(a_vit);
cd_vit = b1 * sin(a_vit).^2 + b2 * cos(a_vit);

% values from the original polar
cl = interp1(alpha, cl, a_vit, 'linear', nan);
cd = interp1(alpha, cd, a_vit, 'linear', nan);

% outside of the bounds, use viterna's
cl(isnan(cl)) = cl_vit(isnan(cl));
cd(isnan(cd)) = cd_vit(isnan(cd));
alpha = a_vit;
