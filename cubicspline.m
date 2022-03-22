


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
    matrix
    
  #matrix for x, y and z
  else
    matrix = zeros(n + 2, n + 5);
  endif
  
  
  
  
  endfunction