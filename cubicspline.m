


function cubicspline
  
  # Clear console
  clc;
  # Get source data
  source("data.m");
  
  #compute the Ts mesh and the output mesh
  n = length(PX);

  
  Ts = linspace(0, n, n);
  outMesh = linspace(0, n, outputNodes);
  Ps = ones(1, n + dimension);
  evaluatedN = length(outMesh);
  
  if(dimension == 2)
    matrix = zeros(n + 2, n + 4);
  else
    matrix = zeros(n + 2, n + 5);
  endif
    
  #rows
  for i = 1 : n
    #columns
    for j = 1 : n + dimension
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
  if(dimension == 3)
    matrix(:, n+5) = resize(PZ,1, n+2);
  endif
  
  #"hardcode" derivatives
  matrix(n+1,3) = 2;
  #Last derivative:
  dV = zeros(1, n+2+dimension);
  dV(1, 3) = 1/3;
  for i = 1 : (n-1)
    dV(1, i+3) = Ts(n) - Ts(i);
  endfor
  dV = dV * 6;
  
  matrix(n+2, :) = dV;
  
  #RREF matrix
  matrix = rref(matrix);
  
  #Get the coefficients for the polynomial in each range
  if(dimension == 2)
  xCoef = matrix(:, (n+3))';
  yCoef = matrix(:, (n+4))';
  else
    xCoef = matrix(:, (n+3))';
    yCoef = matrix(:, (n+4))';
    zCoef = matrix(:, (n+5))';
  endif
  
  finalMatrix = zeros(3, evaluatedN-1);
  
  #evaluate the polynomial for each output node
  for i = 1 : evaluatedN-1
    finalMatrix(1, i) = xCoef(1,1);
    finalMatrix(2, i) = yCoef(1,1);
    if(dimension == 3)
      finalMatrix(3, i) = zCoef(1,1);
    endif
    for j = 2 : n + 2
      if(j <= 4)
        finalMatrix(1, i) += xCoef(1,j) * outMesh(1, i)^(j-1);
        finalMatrix(2, i) += yCoef(1,j) * outMesh(1, i)^(j-1);
        if(dimension == 3)
          finalMatrix(3, i) = zCoef(1,j)* outMesh(1, i)^(j-1);
        endif
      else
        finalMatrix(1, i) += xCoef(1,j) * max(0,(outMesh(1, i) - Ts(1, j - 3))^3);
        finalMatrix(2, i) += yCoef(1,j) * max(0,(outMesh(1, i) - Ts(1, j - 3))^3);
        if(dimension == 3)
          finalMatrix(3, i) = zCoef(1,j)* max(0,(outMesh(1, i) - Ts(1, j - 3))^3);
        endif
      endif
    endfor
  endfor

  if(dimension == 2)
    plot(finalMatrix(1, :), finalMatrix(2,:));
    hold on;
    plot(PX, PY, "or");
  else
    plot3(finalMatrix(1, :), finalMatrix(2, :), finalMatrix(3,:));
    hold on;
    plot3(PX, PY, PZ, "or");
  endif
  
  
  
  
  endfunction