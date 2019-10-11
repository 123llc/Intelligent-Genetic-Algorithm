function y = FitFun(x)
y = 0.01
for i=1:5
    x = 1/(i+(x(i)-1)^2);
end
y = x + y
y = 1/y;