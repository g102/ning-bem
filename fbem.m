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
%   * polar: a struct containing two fields:
%     - polar.CL: an interpolant (gridded or scattered) CL(alpha, radius)
%     - polar.CD: an interpolant (gridded or scattered) CD(alpha, radius)
%     If polar.CL and .CD are functions of only one variable, this is assumed to
%     be the angle of attack (so, for constant airfoil, polar.CL(alpha))
%   * tsr: the tip-speed ratio at which the analysis is carried out
%   * B: the number of blades
%
%  Output:
%   * outr: a struct containing quantities varying along the blade span
%     these are: radius (r), induction factors (a, ap), local aoa (alpha),
%     forces (cfnorm, cftang, cl, cd)
%   * outp: a struct containing performance characteristics
%     these are: thrust (ct), torque (cq), power (cp), root-bending moment (cy)

eps = 1e-6;

% short-hands
r = blade(:, 1);
c = blade(:, 2);
beta = blade(:, 3);

rhub = r(1);
rtip = r(end);
ctip = c(end);

% geometrical functions
solid = griddedInterpolant(r, B*c ./ (2*pi*r));
beta = griddedInterpolant(r, beta);

% preallocation
a = nan(size(r));
ap = nan(size(r));
phi = nan(size(r));
loss = nan(size(r));
for isec = 1:length(r)
	if f(pi/2, r(isec)) > 0
		phistar = fzero(@(x) f(x, r(isec)), [eps, pi/2]);
	elseif fpb(-pi/4, r(isec)) > 0 && fpb(eps, r(isec)) > 0
		phistar = fzero(@(x) f(x, r(isec)), [-pi/4, eps]);
	else
		phistar = fzero(@(x) f(x, r(isec)), [pi/2, pi]);
	end

	phi(isec) = phistar;
	a(isec) = af(phistar, r(isec));
	ap(isec) = apf(phistar, r(isec));
	loss(isec) = F(phistar, r(isec));
end
lsr = tsr * r/rtip;

% reduced quantities
alpha = phi - beta(r);

% blade section coefficients
if length(polar.CL.GridVectors) == 2
	cl = polar.CL(alpha, r);
	cd = polar.CD(alpha, r);
else
	cl = polar.CL(alpha);
	cd = polar.CD(alpha);
end

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
outr = struct('r', r, 'phi', phi, 'a', a, 'ap', ap, 'alpha', alpha, ...
	'cfnorm', cnorm, 'cftang', ctang, 'cfx', cfx, 'cfy', cfy, ...
	'cl', cl, 'cd', cd, 'loss', loss);
% performance
outp = struct('ct', CT, 'cy', CY, 'cq', CQ, 'cp', CP);

% nested functions
	function y = F(phi, r)
		ftip = B/2 .* (rtip - r) ./ (r .* sin(phi));
		Ftip = 2/pi * acos(exp(-ftip));
		fhub = B/2 .* (r - rhub) ./ (r .* sin(phi));
		Fhub = 2/pi * acos(exp(-fhub));
		y = Ftip * Fhub;
	end

	function y = cnf(phi, r)
		aloc = phi - beta(r);
		if length(polar.CL.GridVectors) == 2
			y = polar.CL(aloc, r) .* cos(phi) + polar.CD(aloc, r) .* sin(phi);
		else
			y = polar.CL(aloc) .* cos(phi) + polar.CD(aloc) .* sin(phi);
		end
	end

	function y = ctf(phi, r)
		aloc = phi - beta(r);
		if length(polar.CL.GridVectors) == 2
			y = polar.CL(aloc, r) .* sin(phi) - polar.CD(aloc, r) .* cos(phi);
		else
			y = polar.CL(aloc) .* sin(phi) - polar.CD(aloc) .* cos(phi);
		end
	end

	function y = k(phi, r)
		y = solid(r) .* cnf(phi, r) ./ (4 * F(phi, r) .* (sin(phi)).^2);
		if isinf(y)
			y = realmax;
		end
	end

	function y = kp(phi, r)
 		y = solid(r) .* ctf(phi, r) ./ (4 * F(phi, r) .* sin(phi) .* cos(phi));
		if isinf(y)
			y = -realmax;
		end
	end

	function y = gamma1(phi, r)
		y = 2 * F(phi, r) .* k(phi, r) - (10/9 - F(phi, r));
	end

	function y = gamma2(phi, r)
		y = 2 * F(phi, r) .* k(phi, r) - F(phi, r) .* (4/3 - F(phi, r));
	end

	function y = gamma3(phi, r)
		y = 2 * F(phi, r) .* k(phi, r) - (25/9 - 2 * F(phi, r));
	end

	function y = af(phi, r)
		if k(phi, r) < 2/3
			y = k(phi, r)./(1 + k(phi, r));
		else
			y = (gamma1(phi, r) - sqrt(gamma2(phi, r))) ./ gamma3(phi, r);
		end
	end

	function y = apf(phi, r)
		y = kp(phi, r)./(1 - kp(phi, r));
	end

	function y = f(phi, r)
		y = sin(phi)./(1 - af(phi, r)) - cos(phi).*(1 - kp(phi, r))./(tsr * r/rtip);
	end

	function y = fpb(phi, r)
		y = sin(phi)./(1 - k(phi, r)) - cos(phi).*(1 - kp(phi, r))./(tsr * r/rtip);
	end

end
