function w=rectangleplate(x,y,t, a, b, D, rho, h)
  n=5;
  m=5;
  amplitude=0.01;
  freq=pi*pi*((m/a)*(m/a)+(n/b)*(n/b))*sqrt(D/(rho*h));
  w=amplitude*sin(m*pi*x/a)*sin(n*pi*y/b)*cos(freq*t);
endfunction
