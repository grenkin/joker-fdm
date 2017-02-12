# Comparison of models for different wavelengths
# in the absence and in the presence of nanoparticles
# with taking into account refraction and reflection at the interfaces
# and with the term theta^4

clear all
more off

calc

freq = FREQ_IR;

[grid_info3, temp, intens, Q] = skin_calc(freq, PLANCK_ON, REFR_ON, NANO_OFF);
[grid_info4, temp_nano, intens_nano, Q_nano] = skin_calc(freq, PLANCK_ON, REFR_ON, NANO_ON);

xgrid3 = get_grid_layers(grid_info3);
xgrid4 = get_grid_layers(grid_info4);
L_total = sum(grid_info3.L);
xgrid3_all = concat_cells(xgrid3, 1, 3);
xgrid4_all = concat_cells(xgrid4, 1, 4);


figure(1)
plot(1e6 * xgrid3_all, concat_cells(temp, 1, 3), "b", ...
  1e6 * xgrid4_all, concat_cells(temp_nano, 1, 4), "m");
xlabel('z, micrometers')
ylabel('Temperature, C')
filename = ["4_planck/" int2str(freq) "_temp"];
saveas(1, [filename ".jpg"]);
print_func_layers([filename ".txt"], xgrid3, temp);
print_func_layers([filename "_nano.txt"], xgrid4, temp_nano);

figure(2)
plot(1e6 * concat_cells(xgrid3, 1, 2), 1e-6 * concat_cells(Q, 1, 2), "b", ...
  1e6 * concat_cells(xgrid4, 1, 3), 1e-6 * concat_cells(Q_nano, 1, 3), "m");
xlabel("z, micrometers");
ylabel("Absorbed energy density Q, W / cm^3");
filename = ["4_planck/" int2str(freq) "_Q"];
saveas(2, [filename ".jpg"]);
print_func_layers([filename ".txt"], xgrid3, Q);
print_func_layers([filename "_nano.txt"], xgrid4, Q_nano);

figure(3)
plot(1e6 * concat_cells(xgrid3, 1, 2), concat_cells(intens, 1, 2), "b", ...
  1e6 * concat_cells(xgrid4, 1, 3), concat_cells(intens_nano, 1, 3), "m");
xlabel("z, micrometers");
ylabel('Radiative intensity I_0, W / (m^2 * srad)', "interpreter", "tex");
filename = ["4_planck/" int2str(freq) "_intens"];
saveas(3, [filename ".jpg"]);
print_func_layers([filename ".txt"], xgrid3, intens);
print_func_layers([filename "_nano.txt"], xgrid4, intens_nano);
