# ning-bem
An implementation of Ning's BEM method in Matlab.

This function computes the solution of the BEM problem for the blade given in input, following the approach of Ning (2014) (`doi: 10.1002/we.1636`).
An implementation of Glauert's tip-loss factor is provided.

## Input
 * `blade`: a Nsec × 3 matrix: [radius(:), chord(:), twist(:)]. Radius and chord are in the same units, twist is in RADIANS
 * `polar`: a Nalpha × 3 × Nsec matrix: each page is [alpha(:), cl(:), cd(:)]. alpha is in RADIANS; the j-th page refers to the j-th section of the blade. The polar MUST BE defined over all values of alpha, ie from -pi to pi
 * `tsr`: the tip-speed ratio at which the analysis is carried out
 * `B`: the number of blades
 
## Output
 * `outr`: a struct containing quantities varying along the blade span: these are: radius (r), induction factors (a, ap), local aoa (alpha), forces (cfnorm, cftang, cl, cd)
 * `outp`: a struct containing performance characteristics: these are: thrust (ct), torque (cq), power (cp), root-bending moment (cy)

## Changelog:
v1.01 (2022/07/29):
 * implemented root-loss correction via Prandtl
 * added a function to perform Viterna extrapolation on a minimal foil polar
