function w=rectangleplate(x,y,t, a, b, D, rho, h)
  n=1;
  m=1;
  amplitude=0.01;
  freq=pi*pi*((m/a)*(m/a)+(n/b)*(n/b))*sqrt(D/(pho*h));
  w=amplitude*sin(m*pi*x/a)*sin(n*pi*y/b)*cos(freq*t);
endfunction
