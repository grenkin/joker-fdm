%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Physical data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# L -- optical thickness
# omega -- single scattering albedo
# n -- refractive index
# Nc -- radiation-to-conduction parameter
# A -- anisotropy coefficient

% data for layer 1
layer_data(1) = struct(
  "L",       1,
  "omega", 0.2,
  "n",     1.8,
  "Nc", 0.0001,
  "A",       0);
% data for layer 2
layer_data(2) = struct(
  "L",       1,
  "omega", 0.9,
  "n",       1,
  "Nc",    0.1,
  "A",       0);
% data for layer 3
layer_data(3) = struct(
  "L",       1,
  "omega", 0.2,
  "n",     1.8,
  "Nc", 0.0001,
  "A",       0);

# R -- reflectance
# thetab -- boundary temperature

% data for the left boundary
boundary_data(1) = struct(
  "R",      0.3,
  "thetab",   1);
% data for the right boundary
boundary_data(2) = struct(
  "R",      0.3,
  "thetab", 0.5);
boundary_data(1).n = layer_data(1).n;
boundary_data(2).n = layer_data(3).n;
