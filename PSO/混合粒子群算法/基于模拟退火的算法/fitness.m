function F = fitness(x)
F = 0;
for i=1:5
    F = F +1/(i+(x(i)-1)^2);
end
F = 1/(0.01 + F);