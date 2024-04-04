#=
xxxdocumentation
script to determine integral length scale necessary for turbulentDFSMInlet in openFOAM
=#

using Dates, NCDatasets, DataFrames, Statistics, StatsBase

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

#location of output folder containing bc conditions
outfolder = joinpath(ofcasedir, "constant/boundaryData/inlet/")

#location of TJK-input file
sonicinfile = joinpath(ofcasedir, "scripts/src/tjktmp.nc")

#length of double rotation blocks
drtime = Second(3600)

#flux calculation and avg. time
fxavgtime = Second(60)

#period to use. from t=0 to t=endsec [s], add with spinup time!!
endsec = 3600

#period to use (from t=0 to endidx (20Hz data))
endidx = endsec*20 +1

#timestep for forcing [s] (has to be 0.05 (=20Hz) to properly read measurements)
forcingtimestep = 0.05 #[s]

#vector containing the forcing times (has to be 20Hz!), spinup time gets added!
forcingtimes = collect(0:forcingtimestep:endsec) #collect(0:0.05:3)

#roughness length for log-wind profile and temperature
z0 = 1e-2

#integral length scale at height of TJK (5m) [m]
#estimated from integral time scale and mean wind (code snipplet at end)
integralL = 100

#read in measured turbulence data
turbdata = func.readtotalturbasnetcdf(sonicinfile)

disallowmissing!(turbdata)

#just take part that neededIbe
turbdata = turbdata[1:endidx,:]

#double rotation
func.drdf!(turbdata; blockdur=drtime)

turbdata = turb.despiking(turbdata)
turbdata = turb.interpolatemissing(turbdata)

#calculate fluxes (compontents for Reynolds-stress tensor)
flx = func.turbflux(turbdata, fxavgtime)
flx = func.avgflux(flx, fxavgtime)

#Code for Determining the integral length scale from
#integral time scale and mean wind speed (Taylor hypothesis)

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
ax.set_title("T_E either integral until 0-cross, or time when cross 1/3-line")
ax.plot(autocoru[1:20:end], label="u")
ax.plot(autocorv[1:20:end], label="v")
ax.plot(autocorw[1:20:end], label="w")
PyPlot.axhline(1/ℯ, label="1/e")
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