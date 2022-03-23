# DigiPen Institute of Technology Europe Bilbao
#
# Pablo Riesco, Jon Galaz
# p.riesco@digipen.edu, jon.galaz@digipen.edu
# 3/25/2021
#
# MAT300 Cubic splines project
# cubicspline: curves in 2D and 3D using cubic splines
#===============================================================================
function cubicspline
  # Clear console
  clc;
  # Get source data
  source("data.m");
  
  # n = number of points
  n = length(PX);
  # Mesh of nodes Ts going from 0 to n-1
  Ts = linspace(0, n-1, n);
  # Output mesh from 0 to n-1 with outputNodes values
  outMesh = linspace(0, n-1, outputNodes);
  
  # Matrix used for 
  matrix = zeros(n + 2, n + 2 + dimension);
  
  # Output matrix, where the x, y and z are stored
  finalMatrix = zeros(dimension, outputNodes);
  
  #==============Code Block (1)===================
  # First 4 columns are the Ts^0, Ts^1 ...
  for i = 1 : 4
    matrix(:, i) = resize(Ts.^(i-1), 1, n+2);
  endfor
  #================================================
  
  #==============Code Block (2)===================
  # For the rest of the columns, compute manually
  for i = 3 : n
    for j = 5 : n + 2
      # Only apply addition if the result is positive
      if(Ts(1, i-1) - Ts(1, j-4) > 0)
        matrix(i, j) = (Ts(1, i-1) - Ts(1, j-4))^3;
      else
        continue;
      endif
    endfor
  endfor
  #================================================
  
  #==============Code Block (3)===================
  # Put the x, y and z values in the last columns
  matrix(:, n+3) = resize(PX,1, n+2);
  matrix(:, n+4) = resize(PY,1, n+2);
  if(dimension == 3)
    matrix(:, n+5) = resize(PZ,1, n+2);
  endif
  #================================================
  
  #==============Code Block (4)===================
  # Put in the derivatives on the last rows
  # First row: second derivative using t0 = 0
  matrix(n+1,3) = 2;
  
  # Second row: second derivate using tn = n-1
  matrix(n+2,3) = 2;
  for i = 1 : (n-1)
    matrix(n+2, i+3) = 6 * (Ts(n) - Ts(i));
  endfor
  #================================================
  
  # RREF matrix
  matrix = rref(matrix);
  
  # Get the coefficients for the polynomial
  xCoef = matrix(:, (n+3))';
  yCoef = matrix(:, (n+4))';
  if(dimension == 3)
    zCoef = matrix(:, (n+5))';
  endif
  
  # Evaluate the polynomial for each output node
  #==============Code Block (5)===================
  # Firstly add the first value
  finalMatrix(1, :) = xCoef(1,1);
  finalMatrix(2, :) = yCoef(1,1);
  if(dimension == 3)
    finalMatrix(3, :) = zCoef(1,1);
  endif
  #================================================
  
  #==============Code Block (6)===================
  # Add the second, third and fourth values
  for i = 2 : 4
    finalMatrix(1, :) += xCoef(1,i) * outMesh(1, :).^(i-1);
    finalMatrix(2, :) += yCoef(1,i) * outMesh(1, :).^(i-1);
    if(dimension == 3)
      finalMatrix(3, :) += zCoef(1,i)* outMesh(1, :).^(i-1);
    endif
  endfor
  #================================================
  
  #==============Code Block (7)===================
  # Add the remaining values manually
  for i = 1 : outputNodes
    for j = 5 : n + 2
      if((outMesh(1, i) - Ts(1, j - 3)) > 0)
        finalMatrix(1, i) += xCoef(1,j) * (outMesh(1, i) - Ts(1, j - 3))^3;
        finalMatrix(2, i) += yCoef(1,j) * (outMesh(1, i) - Ts(1, j - 3))^3;
        if(dimension == 3)
          finalMatrix(3, i) += zCoef(1,j)* (outMesh(1, i) - Ts(1, j - 3))^3;
        endif
      else
        continue;
      endif
    endfor
  endfor
  #================================================
  
  # Plot
  if(dimension == 2)
    plot(finalMatrix(1, :), finalMatrix(2,:), "-b");
    hold on;
    plot(PX, PY, "or");
  else
    plot3(finalMatrix(1, :), finalMatrix(2, :), finalMatrix(3,:));
    hold on;
    plot3(PX, PY, PZ, "or");
  endif
  
  endfunction