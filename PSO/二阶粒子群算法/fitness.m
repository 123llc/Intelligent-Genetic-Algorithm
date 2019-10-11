function F = fitness(x)
F = (x(2)-5.1*x(1)^2/4/pi/pi+5*x(1)/pi-6)^2+10*(1-1/8/pi)*cos(x(1))+10;