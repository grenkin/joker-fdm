# Comparison of models with taking into account refraction
# and reflection at the interfaces and without this
# with the term theta^4 and in the absence of nanoparticles

clear all
more off

calc

M = 3;
freq = FREQ_UV

[grid_info, temp_refr, intens_refr, Q_refr] = skin_calc(freq, PLANCK_ON, REFR_ON, NANO_OFF);
[grid_info, temp, intens, Q] = skin_calc(freq, PLANCK_ON, REFR_OFF, NANO_OFF);

n = [1.53, 1.4, 1.4];
I_norm_coeff = 1;
phi_refr_sup = phi_sup = cell(1, M);  
for j = 1 : M
  phi_refr{j} = intens_refr{j} * I_norm_coeff / n(j) ^ 2;
  phi{j} = intens{j} * I_norm_coeff / n(j) ^ 2;
endfor

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
filename = ["2_planck/" int2str(freq) "_temp"];
saveas(1, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, temp_refr);
print_func_layers([filename ".txt"], xgrid, temp);

figure(2)
plot(1e6 * concat_cells(xgrid, 1, 2), concat_cells(intens_refr, 1, 2), "b", ...
  1e6 * concat_cells(xgrid, 1, 2), concat_cells(intens, 1, 2), "g");
xlabel("z, micrometers");
ylabel('Radiative intensity I_0, W / (m^2 * srad)', "interpreter", "tex");
filename = ["2_planck/" int2str(freq) "_intens"];
saveas(2, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, intens_refr);
print_func_layers([filename ".txt"], xgrid, intens);

figure(3)
plot(1e6 * concat_cells(xgrid, 1, 2), 1e-6 * concat_cells(Q_refr, 1, 2), "b", ...
  1e6 * concat_cells(xgrid, 1, 2), 1e-6 * concat_cells(Q, 1, 2), "g");
xlabel("z, micrometers");
ylabel("Absorbed energy density Q, W / cm^3");
filename = ["2_planck/" int2str(freq) "_Q"];
saveas(3, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, Q_refr);
print_func_layers([filename ".txt"], xgrid, Q);

figure(4)
plot(1e6 * concat_cells(xgrid, 1, 2), concat_cells(phi_refr, 1, 2), "b", ...
  1e6 * concat_cells(xgrid, 1, 2), concat_cells(phi, 1, 2), "g");
xlabel("z, micrometers");
ylabel('Normalized intensity');
filename = ["2_planck/" int2str(freq) "_phi"];
saveas(4, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, phi_refr);
print_func_layers([filename ".txt"], xgrid, phi);