# Weighted sums for different wavelengths
# with taking into account refraction and reflection at the interfaces
# and without this in the absence of nanoparticles and without the term theta^4

clear all
more off

calc

M = 3;
freqs = [FREQ_UV, FREQ_VIS, FREQ_IR];
freqs_num = length(freqs);
global freqs_w = [0.05, 0.5, 0.45];

temp_refr = intens_refr = Q_refr = cell(1, freqs_num);
for fi = 1 : freqs_num
  freq = freqs(fi);
  [grid_info, temp_refr{fi}, intens_refr{fi}, Q_refr{fi}] = skin_calc(freq, PLANCK_OFF, REFR_ON, NANO_OFF);
  [grid_info, temp{fi}, intens{fi}, Q{fi}] = skin_calc(freq, PLANCK_OFF, REFR_OFF, NANO_OFF);
endfor

function ret = get_sup (v)
  global freqs_w
  M = length(v{1});
  ret = cell(1, M);
  for j = 1 : M
    ret{j} = v{1}{j} * freqs_w(1) + v{2}{j} * freqs_w(2) + v{3}{j} * freqs_w(3);  
  endfor
endfunction

temp_refr_sup = get_sup(temp_refr);
intens_refr_sup = get_sup(intens_refr);
Q_refr_sup = get_sup(Q_refr);
temp_sup = get_sup(temp);
intens_sup = get_sup(intens);
Q_sup = get_sup(Q);

n = [1.53, 1.4, 1.4];
# StBol = 5.67e-8;
# Tmax = 1
# phi(x) = ( pi / (StBol * n^2) ) * 1 / 2 * \int_{-1}^1 I(x,\nu) d\nu
# I_norm_coeff = pi / StBol;  # normalization coefficient

I_norm_coeff = 1;

phi_refr_sup = phi_sup = cell(1, M);
for j = 1 : M
  phi_refr_sup{j} = intens_refr_sup{j} * I_norm_coeff / n(j) ^ 2;
  phi_sup{j} = intens_sup{j} * I_norm_coeff / n(j) ^ 2;
endfor

xgrid = get_grid_layers(grid_info);
L_total = sum(grid_info.L);
xgrid_all = concat_cells(xgrid, 1, M);

figure(1)
plot(1e6 * xgrid_all, concat_cells(temp_refr_sup, 1, M), "b", ...
  1e6 * xgrid_all, concat_cells(temp_sup, 1, M), "g");
xlabel('z, micrometers')
ylabel('Temperature, C')
xlim([0, 1e6 * L_total])
# ylim([36, 38.5])
filename = "3/temp";
saveas(1, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, temp_refr_sup);
print_func_layers([filename ".txt"], xgrid, temp_sup);

figure(2)
plot(1e6 * concat_cells(xgrid, 1, 2), concat_cells(intens_refr_sup, 1, 2), "b", ...
  1e6 * concat_cells(xgrid, 1, 2), concat_cells(intens_sup, 1, 2), "g");
xlabel("z, micrometers");
ylabel('Radiative intensity I_0, W / (m^2 * srad)', "interpreter", "tex");
filename = "3/intens";
saveas(2, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, intens_refr_sup);
print_func_layers([filename ".txt"], xgrid, intens_sup);

figure(3)
plot(1e6 * concat_cells(xgrid, 1, 2), 1e-6 * concat_cells(Q_refr_sup, 1, 2), "b", ...
  1e6 * concat_cells(xgrid, 1, 2), 1e-6 * concat_cells(Q_sup, 1, 2), "g");
xlabel("z, micrometers");
ylabel("Absorbed energy density Q, W / cm^3");
filename = "3/Q";
saveas(3, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, Q_refr_sup);
print_func_layers([filename ".txt"], xgrid, Q_sup);

figure(4)
plot(1e6 * concat_cells(xgrid, 1, 2), concat_cells(phi_refr_sup, 1, 2), "b", ...
  1e6 * concat_cells(xgrid, 1, 2), concat_cells(phi_sup, 1, 2), "g");
xlabel("z, micrometers");
ylabel('Normalized intensity');
filename = "3/phi";
saveas(4, [filename ".jpg"]);
print_func_layers([filename "_refr.txt"], xgrid, phi_refr_sup);
print_func_layers([filename ".txt"], xgrid, phi_sup);