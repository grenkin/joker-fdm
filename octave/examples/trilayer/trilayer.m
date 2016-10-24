# A complex heat transfer problem for 3 layers

clear all
more off
format long

input_data

printf("Compute with taking into account refraction and reflection\n")
tic
[grid_info, theta, phi] = trilayer_calc(layer_data, boundary_data, true);
toc
printf("Compute without taking into account these effects\n")
tic
[grid_info, theta_wo, phi_wo] = trilayer_calc(layer_data, boundary_data, false);
toc

M = 3;
xgrid = get_grid_layers(grid_info);
L_total = sum(grid_info.L);
xgrid_all = concat_cells(xgrid, 1, M);

figure(1)
plot(xgrid_all, concat_cells(theta, 1, M), "-", ...
  xgrid_all, concat_cells(theta_wo, 1, M), "--");
xlabel('Optical depth')
ylabel('Normalized temperature')
filename = "fig/theta";
saveas(1, [filename ".eps"]);
print_func_layers([filename ".txt"], xgrid, theta);
print_func_layers([filename "_wo.txt"], xgrid, theta_wo);

figure(2)
plot(xgrid{1}, phi{1}, "-", xgrid{2}, phi{2}, "-", xgrid{3}, phi{3}, "-", ...
  xgrid{1}, phi_wo{1}, "--", xgrid{2}, phi_wo{2}, "--", xgrid{3}, phi_wo{3}, "--");
xlabel('Optical depth')
ylabel('Normalized radiative intensity')
filename = "fig/phi";
saveas(2, [filename ".eps"]);
print_func_layers([filename ".txt"], xgrid, phi);
print_func_layers([filename "_wo.txt"], xgrid, phi_wo);
