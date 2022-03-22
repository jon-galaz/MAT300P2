


function cubicspline
  
  # Clear console
  clc;
  # Get source data
  source("data.m");
  
  #compute the Ts mesh and the output mesh
  n = length(PX);
  
  Ts = linspace(0, n, n+1);
  outMesh = linspace(0, n, outputNodes);
  Ps = ones(1, n + 2);
  
  #matrix for x and y 
  if(dimension == 2)
    matrix = zeros(n + 2, n + 4);
    
    #rows
    for i = 1 : n
      #columns
      for j = 1 : n + 2
        if(j <=4 )
          matrix(i, j) = Ps(1,j) * Ts(1, i)^(j-1);
        else
          matrix(i, j) = Ps(1,j) * max(0,(Ts(1, i) - Ts(1, j - 3))^3);
        endif
      endfor
    endfor
    
    #"hardcode" the x and y values
    matrix(:, n+3) = resize(PX,1, n+2);
    matrix(:, n+4) = resize(PY,1, n+2);
    
    #"hardcode" derivatives
    matrix(n+1,3) = 2;
    #Last derivative:
    dV = zeros(1, n+4);
    dV(1, 3) = 1/3;
    for i = 1 : (n-1)
      dV(1, i+3) = Ts(n) - Ts(i);
    endfor
    dV = dV * 6;
    
    matrix(n+2, :) = dV;
    
    matrix
    
    #RREF matrix
    matrix = rref(matrix);
    
    # Output matrix to check result
    matrix
    
  #matrix for x, y and z
  else
    matrix = zeros(n + 2, n + 5);
  endif
  
  temp = matrix(:, (n+3))
  temp = matrix(:, (n+4))
  plot(matrix(:, (n+3)), matrix(:, (n+4)))
  
  
  endfunction