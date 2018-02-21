format long

s = linspace(0.1, 1, 50);
y_star = arrayfun(@solve_eq, s);
psi = (2 + s + 1 ./ s) .* y_star;

[s' y_star' psi']
