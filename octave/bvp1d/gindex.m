# Index of a grid node
# j -- subdomain index
# n -- grid node index within the subdomain
# 1 <= i <= N, 1 <= j <= M, 0 <= n <= K(j)
# N -- number of equations, M -- number of subdomains
# K(j) -- number of grid nodes in the j-th subdomain
function ret = gindex (grid_info, j, n)
  ret = grid_info.cum_index(j) + n + 1;
endfunction
