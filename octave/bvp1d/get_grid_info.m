# Get additional information on the grid
# L -- lengths of subdomains
# K -- numbers of grid intervals
# L, K -- row vectors
# length(L) == length(K) == M (number of subdomains)
# L(j) > 0, K(j) is integer >= 2
function grid_info = get_grid_info (L, K)
  grid_info.K = K;
  grid_info.h = L ./ K;  # step size
  grid_info.nodes = sum(K + 1);  # total number of grid nodes
  # Calculate the first node index in each subdomain
  v = shift(cumsum(K + 1), 1);
  v(1) = 0;
  grid_info.first_index = v + 1;
endfunction
