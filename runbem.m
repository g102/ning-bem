turb = table2array(readtable("./data/turbine.csv"));
foil = table2array(readtable("./data/polar.csv"));

% reinterp blade to n sections
turb = interp1(1:size(turb, 1), turb, linspace(1, size(turb, 1), 100));

[foil360(:,1), foil360(:,2), foil360(:,3)] = viterna(foil(:,1), foil(:,2), foil(:,3));

% stretch the single polar to all sections
alpha = foil360(:, 1);
r = turb(:, 1).';
CL = repmat(foil360(:, 2), [1, length(r)]);
CD = repmat(foil360(:, 3), [1, length(r)]);
polar.CL = griddedInterpolant({alpha, r}, CL);
polar.CD = griddedInterpolant({alpha, r}, CD);

[outr, outp] = fbem(turb, polar, 7, 3);