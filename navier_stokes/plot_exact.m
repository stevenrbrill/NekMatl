function plot_exact(num)

figure(num);
hold on
pts = linspace(-1,1,1000);
plot(0.5*(1-pts.^2),pts,'r-')

end