# Comparison of models with the term theta^4 and without it
# with taking into account refraction and reflection at the interfaces
# in the presence of nanoparticles

clear all
more off

calc

M = 4;
freq = FREQ_IR

[grid_info, temp_pl, intens_pl, Q_pl] = skin_calc(freq, PLANCK_ON, REFR_ON, NANO_ON);
[grid_info, temp, intens, Q] = skin_calc(freq, PLANCK_OFF, REFR_ON, NANO_ON);



xgrid = get_grid_layers(grid_info);
L_total = sum(grid_info.L);
xgrid_all = concat_cells(xgrid, 1, M);

figure(1)
plot(1e6 * xgrid_all, concat_cells(temp_pl, 1, M), "r", ...
  1e6 * xgrid_all, concat_cells(temp, 1, M), "b");
xlabel('z, micrometers')
ylabel('Temperature, C')
xlim([0, 1e6 * L_total])
ylim([36.5, 38.5])
filename = ["1_nano/" int2str(freq) "_temp"];
saveas(1, [filename ".jpg"]);
print_func_layers([filename "_pl.txt"], xgrid, temp_pl);
print_func_layers([filename ".txt"], xgrid, temp);
                                                     
figure(2)
plot(1e6 * concat_cells(xgrid, 1, 3), concat_cells(intens_pl, 1, 3), "r", ...
  1e6 * concat_cells(xgrid, 1, 3), concat_cells(intens, 1, 3), "b");
xlabel("z, micrometers");
ylabel('Radiative intensity I_0, W / (m^2 * srad)', "interpreter", "tex");
filename = ["1_nano/" int2str(freq) "_intens"];
saveas(2, [filename ".jpg"]);
print_func_layers([filename "_pl.txt"], xgrid, intens_pl);
print_func_layers([filename ".txt"], xgrid, intens);

figure(3)
plot(1e6 * concat_cells(xgrid, 1, 3), 1e-6 * concat_cells(Q_pl, 1, 3), "r", ...
  1e6 * concat_cells(xgrid, 1, 3), 1e-6 * concat_cells(Q, 1, 3), "b");
xlabel("z, micrometers");
ylabel("Absorbed energy density Q, W / cm^3");
filename = ["1_nano/" int2str(freq) "_Q"];
saveas(3, [filename ".jpg"]);
print_func_layers([filename "_pl.txt"], xgrid, Q_pl);
print_func_layers([filename ".txt"], xgrid, Q);
