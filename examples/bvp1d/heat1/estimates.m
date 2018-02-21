L = 0.5;
k = 0.0521;
rho = 0.524;
c_p = 1068;
n = 1.0;
h = 10;
eps = 0.8;
Tmax = 773;
sigma = 5.67e-8;
a = k / (rho * c_p);
b = 4 * sigma * n^2 * Tmax^3 / (rho * c_p);
kappa = 10;
kappa_a = 1;
alpha = 1 / (3 * kappa);
beta = h / (rho * c_p);
gamma = eps / (2 * (2 - eps));
c = 3e8;

M = 1;
A = 0.3;
B = 0.3;

d = 0.0113;

k1_lower = min(a / (L^2 / 4 + 1), beta / L);
k1_upper = min(pi^2 * a / (L^2 + pi^2), 2 * beta / L);
k2_lower = min(4 * alpha / L^2, gamma / L);
k2_upper = min(pi^2 * alpha / L^2, 2 * gamma / L);
norm_A1_upper = a + 4 * beta * max(L, 1 / L);
norm_A1_lower = max(pi^2 * a / (L^2 + pi^2), 2 * beta / L);
norm_zeta0_sq = 3 * L * A^2 / 2 + 2 * pi^2 * A^2 / L;  # ||\zeta_0||_V^2
A1_zeta0_zeta0 = 2 * pi^2 * A^2 * a / L;
C1_upper = 1 / (2 * k1_lower) * norm_A1_upper^2 * norm_zeta0_sq ...
  + 4 * b * kappa_a * M^5 * L;
C1_lower = 1 / (2 * k1_upper) * norm_A1_lower^2 * norm_zeta0_sq ...
  + 4 * b * kappa_a * M^5 * L;
C2_upper = 4 * M^6 * kappa_a * C1_upper / k2_lower;
C2_lower = 4 * M^6 * kappa_a * C1_lower / k2_upper;

t_star_upper = 1 / (2 * k2_lower * c) * log(k2_upper * c * d / C2_lower)
norm_diff_sq_upper = 2 * C2_upper / (k2_lower * c) ...
  + 2 * C2_upper / (k2_lower * c) * log(k2_upper * c * d / C2_lower)
r_upper = sqrt(norm_diff_sq_upper / L)

C3_upper = 1 / k1_lower * A1_zeta0_zeta0 + 2 * b * kappa_a * M^5 * L / k1_lower;
C4_upper = 4 * M^6 * kappa_a / k2_lower * ...
  (1 / (2 * k1_lower) * norm_A1_upper^2 * C3_upper + 4 * b * kappa_a * M^5 * L);

s = 0.328;
y_star = solve_eq(s);
norm_diff_sq_upper2 = C4_upper * (2 + s + 1 / s) * y_star / (2 * k2_lower * c)
r_upper2 = sqrt(norm_diff_sq_upper2 / L)