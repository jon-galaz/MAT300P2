


function cubicspline
  
  # Clear console
  clc;
  # Get source data
  source("data.m");
  
  #compute the Ts mesh and the output mesh
  n = length(PX);
  
  Ts = linspace(0, n, n+1);
  outMesh = linspace(0, n, outputNodes);
  
  matrix = zeros(n + 2, n + 3);
  
  for i = 1 : n
    for j = 1 : n + 3
      
    endfor
  endfor
  
  endfunction