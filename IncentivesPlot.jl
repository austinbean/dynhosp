using PyPlot


function funcs(inp::LinSpace)
  outp = -0.5*(inp - 5).^2 + 10
  return outp
end

x1 = linspace(0,5, 1000)
y1 = funcs(x1)
x2 = linspace(5,10, 1000)
y2 = funcs(x2)

int1 = collect(0:5)
out1 = -0.5*(int1-5).^2 + 10

int2 = collect(8:12)
out2 = -0.5*(int2 - 8).^2 + 10

int3 = collect(6:7)
out3 = [10.2, 10.2]

#plot(x1, y1, linestyle = "--", color = "blue")
plot(int1, out1, marker="o", color = "blue", linestyle = "none")
plot(int2, out2, marker="o", color = "blue", linestyle = "none")
plot(int3, out3, marker="o", color = "blue", linestyle = "none")
#plot(x2, y2, linestyle = "--", color = "red")


ylabel("Consumer\nWelfare")
xlabel("Facilities in Market")
ax = axes()
#ax[:yaxis][:set_visible](false)
ax[:yaxis][:set_tick_params](which="major",length=0,width=0)
ax[:xaxis][:set_tick_params](which="major",length=0,width=0)
ax[:yaxis][:set_ticks]([])
ylabel("Consumer\nWelfare")
title("Facility Number vs. Outcome Quality\nTradeoffs")
ax[:text](1, 10, "Quality\nCompetition\nLow,\nVolume\nConcentrated", style="italic")
ax[:text](10, 10, "Quality\nCompetition\nHigh,\nVolume\nDispersed", style="italic")
ax[:text](5.5, 9, "''Optimality''", style="italic")
axis([0, 13, 0, 13]);
savefig("/Users/austinbean/Google Drive/Current Projects/Neonatal Intensive Care Project/Progress Reports/ConsumerWelfare.png")
