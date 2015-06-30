function [m] = baseAlfa4Functions(h)
  a = zeros(6, 6);
  
  h2 = h / 3;
  h3 = h / 6;
  h4 = 8 * h / 15;
  
  a(1, 1) = h2;
  a(1, 3) = h2;
  a(2, 2) = h2;
  a(2, 3) = h2; 
  a(3, 1) = h2; 
  a(3, 2) = h2; 
  a(4, 4) = h2; 
  a(4, 6) = h2; 
  a(5, 5) = h2; 
  a(5, 6) = h2; 
  a(6, 4) = h2; 
  a(6, 5) = h2;
  a(1, 2) = h3;
  a(2, 1) = h3;
  a(4, 5) = h3;
  a(5, 4) = h3;
  a(3, 3) = h4;
  a(6, 6) = h4;

  m=a;
end