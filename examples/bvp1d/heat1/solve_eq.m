function ret = solve_eq (s)
  N = 10;
  y = 10;
  f = @(s, y) (1 + s) * (1 + log(y)) - y;
  df = @(s, y) (1 + s) / y - 1;
  for i = 1 : N
    y = y - f(s, y) / df(s, y);
  endfor
  ret = y;
endfunction