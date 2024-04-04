#=
xxxdocumentation
script for openFOAM-LES-runs to read blockMesh points from binary file
and return boundary files for runs
=#

using Dates, NCDatasets, DataFrames, Statistics, StatsBase, CSV, LsqFit, LaTeXStrings

include("/home/haugened/Documents/ibl_patch_snow/code/src/turb_data.jl")
import .turb

#dir of openFOAM-case
ofcasedir = "/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#location of binary file containing points (in constant/polyMesh)
pointsreadfile = joinpath(ofcasedir, "constant/polyMesh/points")

#location of output file containing the points for the boundary condition
pointsoutfile = joinpath(ofcasedir, "constant/boundaryData/inlet/points")

#location of output folder containing bc conditions
outfolder = joinpath(ofcasedir, "constant/boundaryData/inlet/")

#target for .nc-file containing forcing wind speed
forcingufile = joinpath(ofcasedir, "forcing.nc")

#location of sonic-input files (turbulence data)
tjkinfile =     joinpath(duerrpath, "scripts", "src", "data", "tjktmp.nc")
t1irginfile =   joinpath(duerrpath, "scripts", "src", "data", "t1irgtmp.nc")
t2irginfile =   joinpath(duerrpath, "scripts", "src", "data", "t2irgtmp.nc")
t2lcsatinfile = joinpath(duerrpath, "scripts", "src", "data", "t2lcsattmp.nc")
t2ucsatinfile = joinpath(duerrpath, "scripts", "src", "data", "t2ucsattmp.nc")

#location of slow data input file
tjkslowfile = joinpath(ofcasedir, "scripts", "src", "data", "tjk_data.csv")
t2slowfile = joinpath(ofcasedir, "scripts", "src", "data", "t2slow.nc")

#which input data should be used to force (currently implemented: "tjk", "t2")
datasource = "t2"

#length of double rotation blocks
drtime = Second(300)

#flux calculation (Reynolds averaging time) and avg. time
fxavgtime = Second(50)

#u wind speed for spin up time
spinupspeed = 10.0

#spinuptime = spinupfactor *
#[time it takes the flow with spinupspeed to cross the domain once]
spinupfactor = 0

#starting time of turbulence data (only change, if /scripts/src/tjktmp.nc changes!)
starttime = DateTime(2021,05,31,12,00,00)

#period to use. from t=0 to t=endsec [s], add with spinup time!!
endsec = 300

#calculated end time
endtime = starttime + Second(endsec)

#period to use (from t=0 to endidx (20Hz data))
endidx = endsec*20 +1

#timestep for forcing [s] (has to be 0.05 (=20Hz) to properly read measurements)
forcingtimestep = 0.05 #[s]

#vector containing the forcing times (has to be 20Hz!), spinup time gets added!
forcingtimes = collect(0:forcingtimestep:endsec) #collect(0:0.05:3)

#roughness length (grass) for log-wind profile and temperature
z0 = 5e-3

#integral length scale at height of TJK (5m) [m]
#estimated from integral time scale and mean wind (code snipplet at end)
integralL = 100

#u mult factor
umult = 1

#inlet surface temperature
t_srf_inlet = 293-5

#below only for datasource=="t2"

#height of sensors T2IRG, T2LCSAT, and T2UCSAT
heights = [0.9, 1.9, 2.85] #height of the sensors [m]

#constraints for log-fit (used, if datasource=="t2)
#friction velocity
lower_fric = 1.0e-3
upper_fric = 1.0
#aerodynamic roughness length (z_0)
lower_z0 = 1.0e-5
upper_z0 = 5.0 * z0
#stability (Obukhov length); not used at the moment
#lower_L = 0.0
#upper_L = 100.0

#should fit plot be shown?
plotfit = true

#should scatterplot of inlet points be shown
plotinlet = true

#which sensor to use for scaling after fitting log wind profile to means
#variable is also used for getting Reynolds stresses (assumed constant)
scalevar = "T2IRG" #possible: "T2IRG", "T2LCSAT", and "T2UCSAT"
if scalevar == "T2IRG"
    scalevarheight = heights[1]
elseif scalevar == "T2LCSAT"
    scalevarheight = heights[2]
elseif scalevar == "T2UCSAT"
    scalevarheight = heights[3]
else
    @error("Variable 'scalevar' has wrong value! Exiting...")
    exit()
end

#remove previous files
#rm(uoutfolder, recursive=true, force=true)
##########################################################
#input checks
randidx = rand(collect(2:length(forcingtimes)))
if round(forcingtimes[randidx] - forcingtimes[randidx-1], digits=6) != 0.05
    @error("Timestep needs to be 0.05 (20Hz)!")
end

##########################################################
#deal with points
pts = func.readbinarypoints(pointsreadfile)

frontpts = func.filterfrontpoints(pts)

if plotinlet
    import PyPlot
    PyPlot.pygui(true)
    fig = PyPlot.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.scatter(frontpts[:,2], frontpts[:,3], s=1)
    ax.set_xlabel("y [m]")
    ax.set_ylabel("z [m]")
end

write(pointsoutfile, func.vectoroutstring(frontpts))

#########################################################
#spin up

#determine parameters for spin up
domainlength = round(func.getdomainlength(pts), digits=2)
spinuptime = (domainlength / spinupspeed) * spinupfactor

spinuptimes = collect(0:forcingtimestep:spinuptime)
if length(spinuptimes)==1
    spinuptimes = Float64[]
end

totaltimes = collect(0:forcingtimestep:spinuptime+endsec+forcingtimestep)

if length(spinuptimes) == 0
    totaltimes = totaltimes[1:end-1]
end

#create output arrays
uout = zeros(Float64, size(frontpts, 1), 3, size(spinuptimes, 1) + size(forcingtimes, 1))
reystress = zeros(Float64, size(frontpts, 1), 6, size(spinuptimes, 1) + length(forcingtimes))
intLout = zeros(Float64, size(frontpts, 1), size(spinuptimes, 1) + size(forcingtimes, 1))
Tout = fill(300.0, size(frontpts, 1), size(spinuptimes, 1) + size(forcingtimes, 1))

#fill u for spinup
for itime in 1:length(spinuptimes)
    uout[:, 1, itime] .= func.logscalewindspeed(spinupspeed, 5.0, z0, frontpts[:, 3])
end

#########################################################
#deal with U, R, L, T

if datasource == "tjk"
    @info("TJK data used for wind speed forcing")

    #read in measured turbulence data
    turbdata = func.readtotalturbasnetcdf(tjkinfile)

    disallowmissing!(turbdata)

    #just take part that is needed
    turbdata = turbdata[1:endidx,:]

    #double rotation
    func.drdf!(turbdata; blockdur=drtime)

    print("mean u = ", mean(turbdata.u))

    turbdata = turb.despiking(turbdata)
    turbdata = turb.interpolatemissing(turbdata)

    for itime in 1:length(forcingtimes)
        #uout[:, 1, itime + length(spinuptimes)] .= 5.0 .+ abs.(sin.(pi/length(forcingtimes)*itime)) #func.logscalewindspeed(turbdata.u[itime], 5.0, z0, frontpts[:, 3])
        uout[:, 1, itime + length(spinuptimes)] = round(umult .* func.logscalewindspeed(turbdata.u[itime], 5.0, z0, frontpts[:, 3]), digits=2)
        uout[:, 2, itime + length(spinuptimes)] = round(func.logscalewindspeed(turbdata.v[itime], 5.0, z0, frontpts[:, 3]), digits=2)
        uout[:, 3, itime + length(spinuptimes)] = round(func.logscalewindspeed(turbdata.w[itime], 5.0, z0, frontpts[:, 3]), digits=2)
    end

    func.timedepvectoroutstring(uout, totaltimes, outfolder, "U")

    scaledata = turbdata

elseif datasource == "t2"
    #=
    first fit a logarithmic wind profile (neutral) through the mean over the whole period
    with data from T2IRG, T2LCSAT, and T2UCSAT
    then scale the log wind profile with data used according to variable "scalevar" and "scalevarheight"
    set v- and w- wind-component forcings to 0
    =#
    @info("Tower 2 data used for wind speed forcing")

    #read in measured turbulence data
    turbdata1 = func.readtotalturbasnetcdf(t2irginfile)
    turbdata2 = func.readtotalturbasnetcdf(t2lcsatinfile)
    turbdata3 = func.readtotalturbasnetcdf(t2ucsatinfile)

    disallowmissing!(turbdata1)
    disallowmissing!(turbdata2)
    disallowmissing!(turbdata3)

    #just take part that is needed
    turbdata1 = turbdata1[1:endidx,:]
    turbdata2 = turbdata2[1:endidx,:]
    turbdata3 = turbdata3[1:endidx,:]

    #double rotation
    func.drdf!(turbdata1; blockdur=drtime)
    func.drdf!(turbdata2; blockdur=drtime)
    func.drdf!(turbdata3; blockdur=drtime)

    turbdata1 = turb.despiking(turbdata1)
    turbdata2 = turb.despiking(turbdata2)
    turbdata3 = turb.despiking(turbdata3)
    turbdata1 = turb.interpolatemissing(turbdata1)
    turbdata2 = turb.interpolatemissing(turbdata2)
    turbdata3 = turb.interpolatemissing(turbdata3)

    #convert times from LT (UTC+2h) to UTC+1h
    turbdata1.time .-= Hour(1)
    turbdata2.time .-= Hour(1)
    turbdata3.time .-= Hour(1)

    #fit logarithmic wind profile to mean of data
    logwindprof(z, p) = p[1] / 0.4 * (log.(z / p[2]))# .+ 4.7 * (z / p[3])) #p[1]=u⋆, p[2]=z₀
    p0 = [0.3, 1e-2]#, 0.1] #initial parameters

    means_u = [mean(turbdata1.u), mean(turbdata2.u), mean(turbdata3.u)]

    logwindfit = curve_fit(logwindprof, heights, [mean(turbdata1.u), mean(turbdata2.u), mean(turbdata3.u)], p0)#; lower=[lower_fric,lower_z0,lower_L], upper=[upper_fric, upper_z0, upper_L])
    logwindparam = logwindfit.param
    println("----------------------------------------------------------")
    println("Parameters for log wind profile fit to means of T2-sonics:")
    println(string("u* = ", round(logwindparam[1], digits=2), " m/s"))
    println(string("z0 = ", round(logwindparam[2]*1e3, digits=2), " mm"))
    println("----------------------------------------------------------")
    zdata = sort!(unique(frontpts[:,3]))

    meanlogprof(z) = logwindprof(z, logwindparam)

    if plotfit == true
        import PyPlot
        PyPlot.pygui(true)
        fig = PyPlot.figure()
        fig.suptitle("Fit of logarithmic wind profile to mean of T2-sonics")
        ax = fig.add_subplot(111)
        ax.plot(meanlogprof(zdata), zdata, label="fit")
        ax.plot([mean(turbdata1.u), mean(turbdata2.u), mean(turbdata3.u)], heights, ".", label="sonic means")
        ax.set_xlabel(L"u~\mathrm{[m~s^{-1}]}")
        ax.set_ylabel(L"h~\mathrm{[m]}")
        ax.grid()
        ax.legend()
    end

    if scalevar == "T2IRG"
        scaledata = turbdata1
    elseif scalevar == "T2LCSAT"
        scaledata = turbdata2
    elseif scalevar == "T2UCSAT"
        scaledata = turbdata3
    end

    #evaluate log-wind-profile fit at height of scaling variable
    logwindatscaleheight = meanlogprof(scalevarheight)

    #calculate multiplication factor for each time step
    multfct_u = scaledata.u ./ logwindatscaleheight
    multfct_w = scaledata.w ./ logwindatscaleheight

    for itime in 1:length(forcingtimes)
        uout[:, 1, itime + length(spinuptimes)] = round.(umult .* multfct_u[itime] .* meanlogprof(frontpts[:,3]), digits=2)
        uout[:, 3, itime + length(spinuptimes)] = round.(multfct_w[itime] .* meanlogprof(frontpts[:,3]), digits=4)
    end

    #replace -Inf values
    replace!(uout, -Inf => 0.0)
    replace!(uout, Inf => 0.0)

    #set v- and w-component of forcing wind speed to 0 (double rotated!)
    uout[:,2,:] .= 0.0
    #uout[:,3,:] .= 0.0

    func.timedepvectoroutstring(uout, totaltimes, outfolder, "U")

else
    @error("wrong datasource! Exiting.")
    exit()
end

#calculate fluxes (compontents for Reynolds-stress tensor)
flx = func.turbflux(scaledata, fxavgtime)
flx = func.avgflux(flx, fxavgtime)

#create array containing the Reynolds-stresses
#(1st-dim: points, 2nd Rey-stresses, 3rd time)
#Rey-stresses: uu, uv, uw, vv, vw, ww

for itime in 1:length(totaltimes)
    reystress[:, 1, itime] .= flx.uu[itime] #mean(flx.uu)#[itime]
    reystress[:, 2, itime] .= flx.uv[itime] #mean(flx.uv)#[itime]
    reystress[:, 3, itime] .= flx.uw[itime] #mean(flx.uw)#[itime]
    reystress[:, 4, itime] .= flx.vv[itime] #mean(flx.vv)#[itime]
    reystress[:, 5, itime] .= flx.vw[itime] #mean(flx.vw)#[itime]
    reystress[:, 6, itime] .= flx.ww[itime] #mean(flx.ww)#[itime]
end

func.timedepvectoroutstring(reystress, totaltimes, outfolder, "R")

#scale integral length scale with log-scaling
for itime in 1:length(totaltimes)
    intLout[:,  itime] = func.logscalewindspeed(integralL, 5.0, z0, frontpts[:, 3])
end

func.timedepvectoroutstring(intLout, totaltimes, outfolder, "L")

#read in slow data
tjkmeteodata = func.csvtodataframe(tjkslowfile)
t2vent = turb.readturbasnetcdf(t2slowfile, starttime + Hour(1), starttime + Hour(2)) #t2 is in LT, not in UTC+1
t2vent.time .-= Hour(1)
height_t2vent = 1.00 #[m]
height_hygroVUE = 3.7 #[m]

tjkmeteodata_excerpt = tjkmeteodata[starttime .<= tjkmeteodata.time .<= starttime + Second(endsec), :]
t2vent_excerpt = t2vent[starttime .<= t2vent.time .<= starttime + Second(endsec), :]

meanTair_tjk = mean(tjkmeteodata_excerpt.tair_HygroVUE10)
meanTair_t2vent = mean(t2vent_excerpt.vent_air_temp)

scaledTair = t2vent_excerpt.vent_air_temp .* (meanTair_tjk/meanTair_t2vent)

#=
#correct sonic temperature with humidity (from slow sensor)
#possible: "humidity_rel_CS215", "humidity_rel_ClimaVUE50", "humidity_rel_HygroVUE10"
humidity_var = "humidity_rel_CS215" #checked: they all look pretty similar
#possible: "tair_ClimaVUE50_avg", "tair_CS215_avg", "tair_HygroVUE10"
temperature_var = "tair_CS215_avg"
humidity_data = tjkmeteodata[starttime .<= tjkmeteodata.time .<= endtime, ["time", humidity_var, temperature_var, "barom_pressure_absolute_ClimaVUE50"]]
if size(humidity_data, 1) == 0 #take closest value
    meantime = starttime + (endtime-starttime)/2
    minvector = abs.(tjkmeteodata.time .- meantime)
    minvalue = minimum(minvector)
    closestidx = findfirst(x->x==minvalue, minvector)
    humidity_data = tjkmeteodata[closestidx, ["time", humidity_var, temperature_var, "barom_pressure_absolute_ClimaVUE50"]] 
end

#correction according to Schotanus et al. (1983) in the revision described in van Dijk et al. (2004, eq. 353)
#accessed from: https://www.licor.com/env/support/EddyPro/topics/calculate-flux-level-123.html (equ 6-88)

#following stuff from: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
sat = 611 .* exp.((17.67 .*(humidity_data[:,3]))./(humidity_data[:, 3].-29.65.+273.15))
ws = 0.622 .* (sat ./ (humidity_data[:,4] .*100))
w = humidity_data[:,2]./100 .* ws
q = w ./ (w .+ 1)
q = mean(q)
Tcorrfactor = 1/(1+0.51*q)
turbdata.T = Tcorrfactor .* turbdata.T
=#

#create T field
datetimeforcing1 = starttime + Millisecond(round(Int, forcingtimes[1]*1e3))
for itime in 1:length(forcingtimes)
    if itime == 1
        scaledidx = findfirst(x->x==minimum(abs.(t2vent_excerpt.time .- datetimeforcing1)), abs.(t2vent_excerpt.time .- datetimeforcing1))
    end
    currtime = starttime + Millisecond(round(Int, forcingtimes[itime]*1e3))
    if scaledidx < length(t2vent_excerpt.time)
        if abs(t2vent_excerpt.time[scaledidx+1] - currtime) < abs(t2vent_excerpt.time[scaledidx] - currtime)
            scaledidx += 1
        end
    end
    Tout[:, itime+length(spinuptimes)] = round.(func.logscaletemperature(scaledTair[scaledidx]+273.15, height_hygroVUE, t_srf_inlet, z0, frontpts[:, 3]), digits=2)#height_t2vent, t_srf_inlet, z0, frontpts[:, 3]), digits=2)
end

func.timedepvectoroutstring(Tout, totaltimes, outfolder, "T")

func.saveforcingasnetcdf(frontpts, uout, Tout, totaltimes, forcingufile)

println("Done. Input created.")

#=
#Code for Determining the integral length scale from
#integral time scale and mean wind speed (Taylor hypothesis)

#double rotation
func.drdf!(turbdata; blockdur=Hour(1))

autocoru = autocor(turbdata.u, collect(0:20*1200))
autocorv = autocor(turbdata.v, collect(0:20*1200))
autocorw = autocor(turbdata.w, collect(0:20*1200))

integralautocoru = zeros(Float64, length(autocoru))
integralautocorv = zeros(Float64, length(autocorv))
integralautocorw = zeros(Float64, length(autocorw))

integralautocoru[1] = 1/20*autocoru[1]
integralautocorv[1] = 1/20*autocorv[1]
integralautocorw[1] = 1/20*autocorw[1]

import PyPlot
PyPlot.pygui(true)
a = PyPlot.figure()
ax = a.add_subplot(111)
ax.plot(autocoru[1:20:end], label="u")
ax.plot(autocorv[1:20:end], label="v")
ax.plot(autocorw[1:20:end], label="w")
ax.grid()
ax.legend()
ax.set_xlabel("time [s]")
ax.set_ylabel("autocorrelation")

#numerical integration to yield integral time scale
for i in 2:length(autocoru)
    integralautocoru[i]  = integralautocoru[i-1]+ autocoru[i]/20
end

for i in 2:length(autocorv)
    integralautocorv[i]  = integralautocorv[i-1]+ autocorv[i]/20
end

for i in 2:length(autocorw)
    integralautocorw[i]  = integralautocorw[i-1]+ autocorw[i]/20
end

b = PyPlot.figure()
bx = b.add_subplot(111)
bx.plot(integralautocoru[1:20:end], label="u")
bx.plot(integralautocorv[1:20:end], label="v")
bx.plot(integralautocorw[1:20:end], label="w")
bx.grid()
bx.legend()
bx.set_xlabel("time [s]")
bx.set_ylabel("integral time scale [s]")

=#