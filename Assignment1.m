%%Exercise 5
x0 = .5;
h  = .05; T = 2;
times = linspace(0, T, T/h + 1);
x  = zeros(size(times, 2), 1);
x(1) = x0; 
fxt = @(x,t) x - t ^ 2 + 1;
%RK2
tic
for n=1:size(times,2) - 1
    
    xp = x(n) + h * fxt(x(n), times(n));
    x(n+1) = x(n) + h / 2 * (fxt(x(n), times(n)) + fxt(xp, times(n + 1)));
end
cpuRK2 = toc;
exact = times .^ 2 + 2 .* times + 1 - .5 .* exp(times);
hold on
title('RK2')
error = abs(exact' - x);
%plot(error); 

%RK4
xRK4  = zeros(size(times, 2), 1);
xRK4(1) = x0; 
tic
for n=1:size(times, 2) - 1
    
    xp1 = xRK4(n) + .5 * h * fxt(xRK4(n), times(n));
    xp2 = xRK4(n) + .5 * h * fxt(xp1, times(n) + .5 * h);
    xp3 = xRK4(n) +  h * fxt(xp2, times(n) + .5 * h);
    xRK4(n+1) = xRK4(n) + h * (1/6 * fxt(xRK4(n), times(n)) + 1/3 * fxt(xp1, times(n) + .5 * h) + ...
             1/3 * fxt(xp2, times(n) + .5 * h) + 1/6 * fxt(xp3, times(n) +  h));
end
cpuRK4 = toc;
%figure
%title('RK4')
hold on
errorRK4 = abs(exact' - xRK4);
%plot(errorRK4);


