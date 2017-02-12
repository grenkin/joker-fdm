# Concatenate cells from a to b of cell array c
function ret = concat_cells (c, a, b)
  ret = [];
  for i = a : b
    ret = [ret c{i}];
  endfor  
endfunction