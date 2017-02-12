# Comparison of models with taking into account refraction
# and reflection at the interfaces and without this
# without the term theta^4 and in the presence of nanoparticles

clear all
more off

calc

M = 4;
freq = FREQ_IR

[grid_info, temp_refr, intens_refr, Q_refr] = skin_calc(freq, PLANCK_OFF, REFR_ON, NANO_ON);
[grid_info, temp, intens, Q] = skin_calc(freq, PLANCK_OFF, REFR_OFF, NANO_ON);


xgrid = get_grid_layers(grid_info);
L_total = sum(grid_info.L);
xgrid_all = concat_cells(xgrid, 1, M);

figure(1)
plot(1e6 * xgrid_all, concat_cells(temp_refr, 1, M), "b", ...
  1e6 * xgrid_all, concat_cells(temp, 1, M), "g");
xlabel('z, micrometers')
ylabel('Temperature, C')
xlim([0, 1e6 * L_total])
# ylim([36, 38.5])
filename = ["2_nano/" int2str(freq) "_temp"];
saveas(1, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, temp_refr);
print_func_layers([filename ".txt"], xgrid, temp);

figure(2)
plot(1e6 * concat_cells(xgrid, 1, 3), concat_cells(intens_refr, 1, 3), "b", ...
  1e6 * concat_cells(xgrid, 1, 3), concat_cells(intens, 1, 3), "g");
xlabel("z, micrometers");
ylabel('Radiative intensity I_0, W / (m^2 * srad)', "interpreter", "tex");
filename = ["2_nano/" int2str(freq) "_intens"];
saveas(2, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, intens_refr);
print_func_layers([filename ".txt"], xgrid, intens);

figure(3)
plot(1e6 * concat_cells(xgrid, 1, 3), 1e-6 * concat_cells(Q_refr, 1, 3), "b", ...
  1e6 * concat_cells(xgrid, 1, 3), 1e-6 * concat_cells(Q, 1, 3), "g");
xlabel("z, micrometers");
ylabel("Absorbed energy density Q, W / cm^3");
filename = ["2_nano/" int2str(freq) "_Q"];
saveas(3, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, Q_refr);
print_func_layers([filename ".txt"], xgrid, Q);