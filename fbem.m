function [outr, outp] = fbem(blade, polar, tsr, B)
%  [outr, outp] = FBEM(blade, polar, tsr, B)
%
%  This function computes the solution of the BEM problem for the blade given in
%  input, following the approach of Ning (2014) (doi: 10.1002/we.1636).
%  An implementation of Glauert's tip-loss factor is provided
%
%  Input:
%   * blade: a Nsec × 3 matrix: [radius(:), chord(:), twist(:)]
%     radius and chord are in the same units, twist is in RADIANS
%   * polar: a Nalpha × 3 × Nsec matrix: each page is [alpha(:), cl(:), cd(:)]
%     alpha is in RADIANS; the j-th page refers to the j-th section of the blade
%     the polar MUST BE defined over all values of alpha, ie from -pi to pi
%   * tsr: the tip-speed ratio at which the analysis is carried out
%   * B: the number of blades
%
%  Output:
%   * outr: a struct containing quantities varying along the blade span
%     these are: radius (r), induction factors (a, ap), local aoa (alpha),
%     forces (cfnorm, cftang, cl, cd)
%   * outp: a struct containing performance characteristics
%     these are: thrust (ct), torque (cq), power (cp), root-bending moment (cy)

% define short-hands
r = blade(:, 1);
c = blade(:, 2);
beta = blade(:, 3);

% assuming: r is sorted from root to tip
rtip = r(end);
ctip = c(end);

% solution of the BEM problem
eps = 1e-5;
a = nan(size(r));
ap = nan(size(r));
phi = nan(size(r));
cl = nan(size(r));
cd = nan(size(r));
for isec = 1:length(r)
	solid = B * c(isec) / (2*pi*r(isec));
	lsr = tsr * r(isec) / rtip;

	% tip-loss
	F = @(phi) max(2/pi * acos(exp(-B/2 * (tsr/lsr - 1) / sin(phi))), eps);

	clf = @(alpha) interp1(polar(:, 1, isec), polar(:, 2, isec), alpha);
	cdf = @(alpha) interp1(polar(:, 1, isec), polar(:, 3, isec), alpha);
	cnf = @(phi) clf(phi-beta(isec)) * cos(phi) + cdf(phi-beta(isec)) * sin(phi);
	ctf = @(phi) clf(phi-beta(isec)) * sin(phi) - cdf(phi-beta(isec)) * cos(phi);
	
	k = @(phi) solid * cnf(phi) ./ (4 * F(phi) * sin(phi).^2);
	kp = @(phi) solid * ctf(phi) ./ (4 * F(phi) * sin(phi) * cos(phi));
	
	gamma1 = @(phi) 2 * F(phi) * k(phi) - (10/9 - F(phi));
	gamma2 = @(phi) 2 * F(phi) * k(phi) - F(phi) * (4/3 - F(phi));
	gamma3 = @(phi) 2 * F(phi) * k(phi) - (25/9 - 2 * F(phi));

	af = @(phi) ...
		k(phi)/(1 + k(phi)) * (k(phi) <= 2/3) + ...
		(gamma1(phi) - sqrt(gamma2(phi))) / gamma3(phi) * (k(phi) > 2/3);
	apf = @(phi) kp(phi)/(1 - kp(phi));

	% objective function, regular
	f = @(phi) ...
		sin(phi)/(1 - af(phi)) - cos(phi)*(1 - kp(phi))/lsr;
	
	% objective function, propeller brake
	fpb = @(phi) ...
		sin(phi)/(1 - k(phi)) - cos(phi)*(1 - kp(phi))/lsr;

	% solve 
	if f(pi/2) > 0
		phistar = fzero(f, [eps, pi/2]);
	elseif fpb(-pi/4) > 0 && fpb(eps) > 0
		phistar = fzero(fpb, [-pi/4, -eps]);
	else
		phistar = fzero(f, [pi/2, pi]);
	end
	
	phi(isec) = phistar;
	a(isec) = af(phistar);
	ap(isec) = apf(phistar);
	cl(isec) = clf(phistar-beta(isec));
	cd(isec) = cdf(phistar-beta(isec));
end

lsr = tsr * r/rtip;

% reduced quantities
alpha = phi - beta;

% forces in the turbine frame of reference
cnorm = cl .* cos(phi) + cd .* sin(phi);
cfnorm = cnorm .* c/ctip .* ((1-a).^2 + lsr.^2 .* (1+ap).^2);

ctang = cl .* sin(phi) - cd .* cos(phi);
cftang = ctang .* c/ctip .* ((1-a).^2 + lsr.^2 .* (1+ap).^2);

% forces in the airfoil frame of reference (x: chordwise, y: chord-normal)
cx = cd .* cos(alpha) - cl .* sin(alpha);
cfx = cx .* c/ctip .* ((1-a).^2 + lsr.^2 .* (1+ap).^2);

cy = cd .* sin(alpha) + cl .* sin(alpha);
cfy = cy .* c/ctip .* ((1-a).^2 + lsr.^2 .* (1+ap).^2);

CT = B*ctip / (pi*rtip^2) * trapz(r, cfnorm);
CY = B*ctip / (pi*rtip^3) * trapz(r, r.*cfnorm);

CQ = B*ctip / (pi*rtip^3) * trapz(r, r.*cftang);
CP = tsr*CQ;

% package for output
% radial distributions
outr = struct('r', r, 'a', a, 'ap', ap, 'alpha', alpha, ...
	'cfnorm', cfnorm, 'cftang', cftang, 'cfx', cfx, 'cfy', cfy, ...
	'cl', cl, 'cd', cd);
% performance
outp = struct('ct', CT, 'cy', CY, 'cq', CQ, 'cp', CP);
