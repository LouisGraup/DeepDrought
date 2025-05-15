using LWFBrook90
example = LWFBrook90.run_example()


using Plots, Measures
plotamounts(example, :above_and_belowground, :showRWUcentroid)
plotisotopes(example, :d18O)
plotforcingandstates(example)
    
plot(example.ODESolution;
    idxs = [1, 2, 3, 4, 5, 6],
    label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

x = example.ODESolution_datetime
y = cumsum(example.ODESolution.prob.p.p_soil.p_THICK)
#z = get_theta(example)
z = get_soil_(:theta, example)
zz = Matrix(z[:,2:6])

heatmap(x, y, zz,
    yflip = true,
    xlabel = "Date",
    ylabel = "Depth [mm]",
    c = cgrad(:blues,  rev = false),
    colorbar_title = "θ [m3/m3]\n(of fine soil volume)")#, rightmargin = 10mm)