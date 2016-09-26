# Make Plot


using PyPlot

# For level 1 facilities:
x1 = linspace(0,10,1000);
x2 = linspace(10,25,1000);
y1 = 2.72*ones(size(x1)[1], 1);
y2 = 2.39*ones(size(x2)[1], 1);

# For level 2 facilities:
x3 = linspace(0,10,1000);
x4 = linspace(10, 25, 1000);
x5 = linspace(25, 50, 1000);
y3 = 2.53*ones(size(x3)[1], 1);
y4 = 1.88*ones(size(x4)[1], 1);
y5 = 1.08*ones(size(x5)[1], 1);

# For level 3A

x6 = linspace(0,25,1000);
x7 = linspace(25,50,1000);
x8 = linspace(50,100,1000);
y6 = 1.69*ones(size(x6)[1], 1)
y7 = 1.78*ones(size(x7)[1], 1)
y8 = 1.08*ones(size(x8)[1], 1)

# For level 3B, 3C

x9 = linspace(0,25, 1000);
x10 = linspace(25,50, 1000);
y9 = 1.51*ones(size(x9)[1], 1)
y10 = 1.30ones(size(x10)[1], 1)

# Levels 3B, 3C, 3D

x11 = linspace(50,100,1000);
x12 = linspace(100, 150, 1000);
y11 = 1.19*ones(size(x11)[1], 1);
y12 = 1.0*ones(size(x12)[1], 1);

plot(x1, y1, linestyle = "--", color="blue", label = "Level 1")
plot(x2, y2, linestyle = "--", color = "blue")

plot(x3, y3, linestyle = ":", color = "red", label = "Level 2")
plot(x4, y4, linestyle = ":", color = "red")
plot(x5, y5, linestyle = ":", color = "red")

plot(x6, y6, linestyle = "-.", color = "green", label = "Level 3A")
plot(x7, y7, linestyle = "-.", color = "green")
plot(x8, y8, linestyle = "-.", color = "green")

plot(x9, y9, linestyle = "--", color = "black", label = "Level 3B or 3C")
plot(x10, y10, linestyle = "--", color = "black")
plot(x11, y11, linestyle = "-", color = "black", label = "Level 3B, 3C, 3D")
plot(x12, y12, linestyle = "-", color = "black")
legend()
title("Mortality Odds Ratio Relative to High Volume, High Level Facility 3B, 3C, 3D\nBased on Baker and Phibbs NEJM 2007")
xlabel("Annual Patient Volume")
ylabel("Odds Ratio for Mortality Relative to \nHighest Volume (>100) 3B, 3C, 3D")
axis([0, 150, 0.75, 3]);
savefig("/Users/austinbean/Google Drive/Current Projects/Neonatal Intensive Care Project/Progress Reports/RelativeRiskPlot.png")
