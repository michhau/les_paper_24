#=
script to generate the figures for the LES paper
=#

####################################################################################
####################################################################################
###  Fig. 1:            OVERVIEW MAP WITH FLUX FOOTPRINTS                        ###
####################################################################################
####################################################################################
#=
using Dates, PyCall, DataFrames, Statistics, LaTeXStrings, ProgressMeter, Distributed
import PyPlot, CSV
pydates = pyimport("matplotlib.dates")
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
LogNorm = pyimport("matplotlib.colors")
mpimg = pyimport("matplotlib.image")
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    datapath = "/home/michi/Documents/slf/data/"
    tjkpath = "/home/michi/Documents/slf/data/tjk/tjk_data.csv"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    datapath = "/home/haugened/Documents/data/"
    tjkpath = "/home/haugened/Documents/data/tjk/tjk_data.csv"
end
include(joinpath(importdir, "src", "turb_data.jl"))
include(joinpath(importdir, "src", "general.jl"))
include(joinpath(importdir, "src", "kljun_ffp.jl"))
import .turb
import .gen
import .kljun
PyPlot.pygui(true)

timestep = Millisecond(50)
ρ_air = 1.2 #kg m^{-3}
c_p = 1004 #J kg^{-1} K^{-1}
L_v = 2450e3 #J kg^{-1} (approx @0°C) 

#timestep between single measurements, 1/measurement frequency
timestep = Millisecond(50)
kaijo_period_file = joinpath(datapath, "kaijo", "21_kaijo_periods.txt")
kaijo_outfile_stam = joinpath(datapath, "kaijo", "kaijo")
outfile_stam = joinpath(datapath, "WindSonic_Duerrboden_2005xx/tjk_sonic_200504-200514")
tower_outfile_stam = joinpath(datapath, "tower", "preproc")
ventair_stam = joinpath(datapath, "tower", "vent_air")

#select data and measurement period to be evaluated
evalstart = DateTime(2021, 05, 31, 13, 00, 00)
evalend   = DateTime(2021, 05, 31, 13, 05, 00)
#evalend = evalstart + Day(10)

evaldf1 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t1irgdb.nc"), evalstart, evalend)
evaldf2 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2irgdb.nc"), evalstart, evalend)
evaldf3 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2lcsatdb.nc"), evalstart, evalend)
evaldf4 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "t2ucsatdb.nc"), evalstart, evalend)
evaldf5 = turb.readturbasnetcdf(string(kaijo_outfile_stam, ".nc"), evalstart, evalend)
evaldf6 = turb.readturbasnetcdf(joinpath(tower_outfile_stam, "tjkdf.nc"), evalstart, evalend)

#rotate kaijo data due to mounting. Be careful and check for every new position!
if @isdefined evaldf5
    @warn("Careful!! Rotate Kaijo data due to measurement position. Check pictures of setup!")
    new_w = evaldf5.u
    new_u = -evaldf5.w
    evaldf5.u = new_u
    evaldf5.w = new_w
end

#apply NaN-mask to T1 & T2 when repositioned (for DR)
turb.repositionnanmask!(evaldf1)
turb.repositionnanmask!(evaldf2)
turb.repositionnanmask!(evaldf3)
turb.repositionnanmask!(evaldf4)

#show statistics about missing data
turb.printmissstats(evaldf1)
turb.printmissstats(evaldf2)
turb.printmissstats(evaldf3)
turb.printmissstats(evaldf4)
turb.printmissstats(evaldf5)
turb.printmissstats(evaldf6)
#interpolate missing
evaldf1 = turb.interpolatemissing(evaldf1)
evaldf2 = turb.interpolatemissing(evaldf2)
evaldf3 = turb.interpolatemissing(evaldf3)
evaldf4 = turb.interpolatemissing(evaldf4)
try
    evaldf5 = turb.interpolatemissing(evaldf5)
catch LoadError
    @warn("Kaijo DataFrame empty")
end
evaldf6 = turb.interpolatemissing(evaldf6)
#double rotation
turb.drdf!(evaldf1)
turb.drdf!(evaldf2)
turb.drdf!(evaldf3)
turb.drdf!(evaldf4)
turb.drdf!(evaldf5)
turb.drdf!(evaldf6, periodwise=false)


tjkmeteodata = turb.csvtodataframe(tjkpath)
tjkmeteodata = tjkmeteodata[evalstart.<=tjkmeteodata.time.<=evalend, :]


turb.missing2nan!(evaldf1)
turb.missing2nan!(evaldf2)
turb.missing2nan!(evaldf3)
turb.missing2nan!(evaldf4)
turb.missing2nan!(evaldf5)
turb.missing2nan!(evaldf6)

#Reynolds averaging times
ra1 = Millisecond(2^11*timestep) #Second(330)
ra2 = Millisecond(2^11*timestep) #Second(330)
ra3 = Millisecond(2^11*timestep) #Second(320)
ra4 = Millisecond(2^11*timestep) #Second(210)
ra5 = Millisecond(2^9*timestep) #Second(300)
ra6 = Millisecond(2^10*timestep) #Second(55)

fx1_raw = turb.turbflux(evaldf1, ra1)
fx2_raw = turb.turbflux(evaldf2, ra2)
fx3_raw = turb.turbflux(evaldf3, ra3)
fx4_raw = turb.turbflux(evaldf4, ra4)
fx5_raw = turb.turbflux(evaldf5, ra5)
fx6_raw = turb.turbflux(evaldf6, ra6)

#averaging
fx1 = turb.avgflux(fx1_raw, Second(600))
fx2 = turb.avgflux(fx2_raw, Second(600))
fx3 = turb.avgflux(fx3_raw, Second(600))
fx4 = turb.avgflux(fx4_raw, Second(600))
fx5 = turb.avgflux(fx5_raw, Second(600))
fx6 = turb.avgflux(fx6_raw, Second(600))

#Obukhov-length
L1 = turb.obukhov(fx1_raw, Minute(30))
L2 = turb.obukhov(fx2_raw, Minute(10))
L3 = turb.obukhov(fx3_raw, Minute(10))
L4 = turb.obukhov(fx4_raw, Minute(10))
L5 = turb.obukhov(fx5_raw, Minute(30))
L6 = turb.obukhov(fx6_raw, Minute(10))

#variables
names = [:evaldf1, :evaldf2, :evaldf3, :evaldf4, :evaldf5, :evaldf6]
meas_heights = [1.2, 0.9, 1.9, 2.8, 0.3, 5]
pbl_height = 1000
Ls = [:L1, :L2, :L3, :L4, :L5, :L6]
nrelemgrid = 1000
fluxes = [:fx1, :fx2, :fx3, :fx4, :fx5, :fx6]
rs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] #levels for plotting
winddir = 330#mean(tjkmeteodata.wind_mean_vector_direction)#0 #rotating the footprint
rslayer = false #measurement within roughness sublayer (theory not working properly)
crop = true #crop output to maximum defined rs (max 0.9)
outnames = [:ffp1, :ffp2, :ffp3, :ffp4, :ffp5, :ffp6]

for ix in [1, 2, 5] #1:size(names, 1)
    println("Calculating footprint for ", String(names[ix]))
    turb.missing2nan!(@eval $(names[ix]))
    (x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d, rs, fr, xr, yr, flag_err) = kljun.ffp(
        meas_heights[ix], nothing, 5, pbl_height, mean(filter(!isnan, @eval $(Ls[ix]).L)),
        std(filter(!isnan, @eval $(names[ix]).v)), mean(filter(!isnan, @eval $(fluxes[ix]).u_star)), winddir, rs, rslayer, nrelemgrid, crop)
    if flag_err == true
        @warn("Error flag set to true! There is an error!")
    end
    @eval $(outnames[ix]) = $(Dict("x_ci_max" => x_ci_max, "x_ci" => x_ci, "f_ci" => f_ci,
        "x_2d" => x_2d, "y_2d" => y_2d, "f_2d" => f_2d, "rs" => rs, "fr" => fr,
        "xr" => xr, "yr" => yr, "flag_err" => flag_err))
end

#plotting the footprint on the ortho-mosaic

fileorthomosaic = joinpath("/home/haugened/Documents/les_paper_24/figures/site_overview_raw", "210530_bea_new.jpg")
orthomosaic = mpimg.imread(fileorthomosaic)
#PyPlot.imshow(orthomosaic)
#location of flux measurements 1-6 in original image
#[row-location, col-location]

##
fluxloc = [710 495; 804 542; 804 542; 804 542; 714 402; 249 472]

#extend of background [row, col]
#bgextend_m = [280, 280] #in m from measuring in GIS: 279.9
bgextend_pxl = [orthomosaic.shape[1], orthomosaic.shape[2]] #in pxl

#calculate m/pxl from it
meterperpxl_row = 0.1 #bgextend_m[1] / bgextend_pxl[1]
meterperpxl_col = meterperpxl_row #bgextend_m[2] / bgextend_pxl[2]

#origin of figure
figorigin = [714 402] #tower 2

#calculate fluxloc in new coordinates [m]
fluxloc_final = Array{Float64}(undef, size(fluxloc, 1), size(fluxloc, 2))
fluxloc_final[:, 1] = (figorigin[1] .- fluxloc[:, 1]) .* meterperpxl_row
fluxloc_final[:, 2] = (fluxloc[:, 2] .- figorigin[2])  .* meterperpxl_col

#calculate extend in new coordinates
lft = (-figorigin[2]) * meterperpxl_col
rght = (bgextend_pxl[2]-1-figorigin[2]) *meterperpxl_col
tp = (bgextend_pxl[1]-(figorigin[1]- bgextend_pxl[1])-1) * meterperpxl_row
btm = (figorigin[1] - bgextend_pxl[1]) * meterperpxl_row

bgextend_final = (-figorigin[2], bgextend_pxl[2]-1-figorigin[2], -(bgextend_pxl[1]-figorigin[1]), bgextend_pxl[1]-(bgextend_pxl[1]-figorigin[1])-1).*meterperpxl_col
#(lft, rght, btm, tp)


ctab10 = PyPlot.cm.tab10
ffp_fig = PyPlot.figure(figsize=(15,15))
ax1 = ffp_fig.add_subplot(111)
#ax1.set_title("Flux footprints 80% - Kljun et al. (2015)")
bg = ax1.imshow(orthomosaic, extent=bgextend_final)
#bg = ax1.pcolormesh(orthomosaic)
ax1.set_xlabel("meter")
ax1.set_ylabel("meter")
locfx1 = ax1.plot(fluxloc_final[1, 2], fluxloc_final[1, 1], ".", markersize=40, color=ctab10(0), alpha=0.7)#, label="T1IRG")
#orig = ax1.plot(0,0, ".")
locfx2 = ax1.plot(fluxloc_final[2, 2], fluxloc_final[2, 1], ".", markersize=40, color=ctab10(2), alpha=0.7)#, label="T2IRG")
#locfx3 = ax1.plot(fluxloc_final[3, 2], fluxloc_final[3, 1], ".", color=ctab10(2))#, label="T2LCSAT")
#locfx4 = ax1.plot(fluxloc_final[4, 2], fluxloc_final[4, 1], ".", color=ctab10(3))#, label="T2UCSAT")
#locfx5 = ax1.plot(fluxloc_final[5, 2], fluxloc_final[5, 1], ".", markersize=40, color=ctab10(3), alpha=0.7)#, label="Kaijo")
locfx6 = ax1.plot(fluxloc_final[6, 2], fluxloc_final[6, 1], ".", markersize=40, color=ctab10(3), alpha=0.7)#, label="TJK")
fp1 = ax1.plot(ffp1["xr"][:, end-1] .+ fluxloc_final[1, 2], ffp1["yr"][:, end-1] .+ fluxloc_final[1, 1], color=ctab10(0), label = "T1IRG", lw=6, alpha=0.7)
fp2 = ax1.plot(ffp2["xr"][:, end-1] .+ fluxloc_final[2, 2], ffp2["yr"][:, end-1] .+ fluxloc_final[2, 1], color=ctab10(2), label = "T2IRG", lw=6, alpha=0.7)
#fp3 = ax1.plot(ffp3["xr"][:, end-1] .+ fluxloc_final[3, 2], ffp3["yr"][:, end-1] .+ fluxloc_final[3, 1], color=ctab10(2), label = "T2LCSAT")
#fp4 = ax1.plot(ffp4["xr"][:, end-1] .+ fluxloc_final[4, 2], ffp4["yr"][:, end-1] .+ fluxloc_final[4, 1], color=ctab10(3), label = "T2UCSAT")
#fp5 = ax1.plot(ffp5["xr"][:, end-1] .+ fluxloc_final[5, 2], ffp5["yr"][:, end-1] .+ fluxloc_final[5, 1], color=ctab10(3), label = "Kaijo", lw=6, alpha=0.7)
#fp6 = ax1.plot(ffp6["xr"][:, end-1] .+ fluxloc_final[6, 2], ffp6["yr"][:, end-1] .+ fluxloc_final[6, 1], color=ctab10(5), label = "TJK")
#ax1.legend()
PyPlot.axis("off")
#PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/210530_new_ffp.jpg")
##
=#
####################################################################################
####################################################################################
###  Fig. 2:         MAP WITH LES DOMAIN AND SENSOR LOCATIONS                    ###
####################################################################################
####################################################################################
#see GIMP file

####################################################################################
####################################################################################
###  Fig. 3:   TIME SERIES AND SPECTRA COMPARISON: EC MEASUREMENT VS LES         ###
####################################################################################
####################################################################################
#xxx check time steps!
#=
using NCDatasets, DataFrames, StatsBase, LaTeXStrings, Dates, PyCall
import PyPlot
PyPlot.pygui(true)
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")

#outer directory of case
casedir1 = "/home/haugened/Documents/openfoam/real_run_state/cscs/third_run/"
casedir2 = "/home/haugened/Documents/openfoam/real_run_state/cscs/first_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#netCDF4-file containing the profile (length 3)
file1 = joinpath.(joinpath(casedir1, "postProcessing/surfaces/yNormal/"), "yNormal_x5cm.nc")
file2 = joinpath.(joinpath(casedir2, "postProcessing/surfaces/yNormal/"), "yNormal_new_total_time10th.nc")
#location of .nc-file containing the forcing wind speed
uforcingfile1 = joinpath(casedir1, "forcing.nc")
uforcingfile2 = joinpath(casedir2, "forcing.nc")

#axis along which profile is taken ("x", "y", or "z")
profilealong = "z"

#value of x where bare to snow transition is
snowtransition1 = 8.0 #m
snowtransition2 = 20.0 #m

#seconds to skip (spinuptime)
skipsec1 = 10
skipsec2 = 10

#reynolds averaging time
reyavgtime = Second(50)

#flag for filtering of profiles
postfilter = false
postfilterquantity = "u" #u, v, w, or T
#location of point to evaluate for post-filtering (fetch distance coordinates for x) [m]
postfilterloc = [-1.0, 2, 0.5]
#filter condition
filtercond = "0.0 .<= filterdata .<= filterquantiles[50]"

#locations of the point data (fetch distance coordinates for x) [m]
p1 = [2, 5.0, 0.3]
p2 = [10.0, 5, 1.9]
p3 = [10.0, 5.0, 1.2]

#height to choose from forcing data (if ==-1, then use same as first point) [m]
forcingheight = -1

#flag for loading sonic measurement data as well
loadsonics = true
#which sonics
listsonics = ["kaijo", "t1irg", "t2irg", "t2lcsat", "t2ucsat", "tjk"]
allowedsonics = ["kaijo", "t1irg", "t2irg", "t2lcsat", "t2ucsat", "tjk"]
#period to load from sonics
sonicsstart = DateTime(2021,05,31,12,00,00)
sonicsend = sonicsstart + Minute(5)
drtime = Second(300)

#location of sonic-input files (turbulence data)
tjkinfile =     joinpath(duerrpath, "scripts", "src", "data", "tjktmp.nc")
t1irginfile =   joinpath(duerrpath, "scripts", "src", "data", "t1irgtmp.nc")
t2irginfile =   joinpath(duerrpath, "scripts", "src", "data", "t2irgtmp.nc")
t2lcsatinfile = joinpath(duerrpath, "scripts", "src", "data", "t2lcsattmp.nc")
t2ucsatinfile = joinpath(duerrpath, "scripts", "src", "data", "t2ucsattmp.nc")
kaijoinfile =   joinpath(duerrpath, "scripts", "src", "data", "kaijotmp.nc")

#################################################
#load additional functions if sonics should be included
if loadsonics
    if !@isdefined turb_loaded
        turb_loaded = false
    end
    if !turb_loaded
        include("/home/haugened/Documents/ibl_patch_snow/code/src/turb_data.jl")
        import .turb    
        turb_loaded = true
    end
end
#################################################

if loadsonics
    for i in listsonics
        if !(i in allowedsonics)
            println(string(i, " is not in ", allowedsonics, ". Please set variable 'listsonics' to allowed value. Exiting."))
            exit()
        end
    end
    println("Loading sonic data from ", listsonics)
end

#################################################
#read in data
#=
if !@isdefined forcingdata1_read
    forcingdata1_read = false
    forcingdata1_read_sucessfully = false
end
if !forcingdata1_read
    println("Reading forcing data 1 ...")
    if isfile(uforcingfile1)
        (forcing1_heights, forcing1_times, forcing1_u, forcing1_T) = func.readforcingfromnetcdf(uforcingfile1)
        forcingdata1_read_sucessfully = true
        @info("Forcing data 1 sucessfully read.")
    else
        @warn("Forcing data 1 could not be read.")
    end
    forcingdata1_read = true
end
=#
#=
if !@isdefined forcingdata2_read
    forcingdata2_read = false
    forcingdata2_read_sucessfully = false
end
if !forcingdata2_read
    println("Reading forcing data 2 ...")
    if isfile(uforcingfile2)
        (forcing2_heights, forcing2_times, forcing2_u, forcing2_T) = func.readforcingfromnetcdf(uforcingfile2)
        forcingdata2_read_sucessfully = true
        @info("Forcing data 2 sucessfully read.")
    else
        @warn("Forcing data 2 could not be read.")
    end
    forcingdata2_read = true
end
=#
if !@isdefined data1_read
    data1_read = false
end
if !data1_read
    println("Reading cross section data 1 ...")
    (x1, y1, z1, times1, T1, U1) = func.readaerialprofilefromnetcdf(file1)
    data1_read = true
    T1 .-= 273.15
end
#=
if !@isdefined data2_read
    data2_read = false
end
if !data2_read
    println("Reading cross section data 2 ...")
    (x2, y2, z2, times2, T2, U2) = func.readaerialprofilefromnetcdf(file2)
    data2_read = true
end
=#
#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx1
    transformx1 = false
end
if !transformx1
    x1 = x1 .- snowtransition1
    transformx1 = true
end
#=
if !@isdefined transformx2
    transformx2 = false
end
if !transformx2
    x2 = x2 .- snowtransition2
    transformx2 = true
end
=#
skiptoidx1 = findfirst(x->x==minimum(abs.(times1 .- skipsec1)), abs.(times1 .- skipsec1))
#skiptoidx2 = findfirst(x->x==minimum(abs.(times2 .- skipsec2)), abs.(times2 .- skipsec2))
println(string("Skipping ", skipsec1, " seconds at the beginning of 1."))
#println(string("Skipping ", skipsec2, " seconds at the beginning of 2."))

#deal with sonic data
if loadsonics
    if !@isdefined sonics_read
        sonics_read = false
    end
    if !sonics_read
        #sonic heights
        sonic_heights = zeros(Float64, length(listsonics))
        println("Reading sonic data...")
        if "kaijo" in listsonics
            htmp = 0.3
            datatmp = func.readtotalturbasnetcdf(kaijoinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            @warn("Careful!! Rotate Kaijo data due to measurement position. Check pictures of setup!")
            new_w = datatmp.u
            new_u = datatmp.w
            datatmp.u = new_u
            datatmp.w = new_w
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            kaijodata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="kaijo", listsonics)] = htmp
        end
        if "t1irg" in listsonics
            htmp = 1.2
            datatmp = func.readtotalturbasnetcdf(t1irginfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            #flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t1irgdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t1irg", listsonics)] = htmp
        end
        if "t2irg" in listsonics
            htmp = 0.9
            datatmp = func.readtotalturbasnetcdf(t2irginfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2irgdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2irg", listsonics)] = htmp
        end
        if "t2lcsat" in listsonics
            htmp = 1.9
            datatmp = func.readtotalturbasnetcdf(t2lcsatinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2lcsatdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2lcsat", listsonics)] = htmp
        end
        if "t2ucsat" in listsonics
            htmp = 2.85
            datatmp = func.readtotalturbasnetcdf(t2ucsatinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2ucsatdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2ucsat", listsonics)] = htmp
        end
        if "tjk" in listsonics
            htmp = 5.0
            datatmp = func.readtotalturbasnetcdf(tjkinfile)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            tjkdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="tjk", listsonics)] = htmp
        end
        sonics_read = true
    end
end

#################################################
#correct sonic temperature with humidity (from T2IRG)
humidity_data = t2irgdata.h2o
#humidity_data = tjkmeteodata[starttime .<= tjkmeteodata.time .<= endtime, ["time", humidity_var, temperature_var, "barom_pressure_absolute_ClimaVUE50"]]
#if size(humidity_data, 1) == 0 #take closest value
#    meantime = starttime + (endtime-starttime)/2
#    minvector = abs.(tjkmeteodata.time .- meantime)
#    minvalue = minimum(minvector)
#    closestidx = findfirst(x->x==minvalue, minvector)
#    humidity_data = tjkmeteodata[closestidx, ["time", humidity_var, temperature_var, "barom_pressure_absolute_ClimaVUE50"]] 
#end

#correction according to Schotanus et al. (1983) in the revision described in van Dijk et al. (2004, eq. 353)
#accessed from: https://www.licor.com/env/support/EddyPro/topics/calculate-flux-level-123.html (equ 6-88)

#following stuff from: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
#sat = 611 .* exp.((17.67 .*(humidity_data[:,3]))./(humidity_data[:, 3].-29.65.+273.15))
#ws = 0.622 .* (sat ./ (humidity_data[:,4] .*100))
#w = humidity_data[:,2]./100 .* ws
#q = w ./ (w .+ 1)
q = mean(humidity_data./1000)
#q=0.008
Tcorrfactor = 1/(1+0.51*q)
#turbdata.T = Tcorrfactor .* turbdata.T

if "kaijo" in listsonics
    kaijodata.T = ((kaijodata.T .+ 273.15).*Tcorrfactor) .-273.15
end
if "t1irg" in listsonics
    q = t1irgdata.h2o ./ 1000
    Tcorrfactor = 1 ./ (1 .+ 0.51 .* q)
    t1irgdata.T = ((t1irgdata.T .+ 273.15) .*Tcorrfactor) .-273.15
end
if "t2irg" in listsonics
    q = t2irgdata.h2o ./ 1000
    Tcorrfactor = 1 ./ (1 .+ 0.51 .* q)
    t2irgdata.T = ((t2irgdata.T .+ 273.15).*Tcorrfactor) .-273.15
end
if "t2lcsat" in listsonics
    t2lcsatdata.T = ((t2lcsatdata.T .+ 273.15).*Tcorrfactor) .-273.15
end
if "t2ucsat" in listsonics
    t2ucsatdata.T = ((t2ucsatdata.T .+ 273.15).*Tcorrfactor) .-273.15
end
if "tjk" in listsonics
    tjkdata.T = ((tjkdata.T .+ 273.15).*Tcorrfactor) .-273.15
end


#################################################
take1_1 = func.chooseval3D(p1, hcat(x1, y1, z1))
take2_1 = func.chooseval3D(p2, hcat(x1, y1, z1))
take3_1 = func.chooseval3D(p3, hcat(x1, y1, z1))
#=
take1_2 = func.chooseval3D(p1, hcat(x2, y2, z2))
take2_2 = func.chooseval3D(p2, hcat(x2, y2, z2))
take3_2 = func.chooseval3D(p3, hcat(x2, y2, z2))
=#

println(string("run 1, point1 (x, y, z) [m] = (", x1[take1_1], ", ", y1[take1_1], ", ", z1[take1_1], ")"))
println(string("run 1, point2 (x, y, z) [m] = (", x1[take2_1], ", ", y1[take2_1], ", ", z1[take2_1], ")"))
println(string("run 1, point3 (x, y, z) [m] = (", x1[take3_1], ", ", y1[take3_1], ", ", z1[take3_1], ")"))
#=
println(string("run 2, point1 (x, y, z) [m] = (", x2[take1_2], ", ", y2[take1_2], ", ", z2[take1_2], ")"))
println(string("run 2, point2 (x, y, z) [m] = (", x2[take2_2], ", ", y2[take2_2], ", ", z2[take2_2], ")"))
println(string("run 2, point3 (x, y, z) [m] = (", x2[take3_2], ", ", y2[take3_2], ", ", z2[take3_2], ")"))
=#
u1_1 = U1[take1_1, skiptoidx1:end, 1]
u2_1 = U1[take2_1, skiptoidx1:end, 1]
u3_1 = U1[take3_1, skiptoidx1:end, 1]
w1_1 = U1[take1_1, skiptoidx1:end, 3]
w2_1 = U1[take2_1, skiptoidx1:end, 3]
w3_1 = U1[take3_1, skiptoidx1:end, 3]
T1_1 = T1[take1_1, skiptoidx1:end]
T2_1 = T1[take2_1, skiptoidx1:end]
T3_1 = T1[take3_1, skiptoidx1:end]
times_points_1 = times1[skiptoidx1:end]

#heat flux data
timesteps1 = func.gettimesteps(times_points_1)
meantimesteps1 = round(mean(timesteps1), digits=3)
movavgidcs1 = round(Int, reyavgtime/Millisecond(meantimesteps1*1e3))

movavgT1_1 = func.movingaverage(T1_1, movavgidcs1)
movavgT2_1 = func.movingaverage(T2_1, movavgidcs1)
movavgT3_1 = func.movingaverage(T3_1, movavgidcs1)

movavgw1_1 = func.movingaverage(w1_1, movavgidcs1)
movavgw2_1 = func.movingaverage(w2_1, movavgidcs1)
movavgw3_1 = func.movingaverage(w3_1, movavgidcs1)

devT1_1 = T1_1 .- movavgT1_1
devT2_1 = T2_1 .- movavgT2_1
devT3_1 = T3_1 .- movavgT3_1

devw1_1 = w1_1 .- movavgw1_1
devw2_1 = w2_1 .- movavgw2_1
devw3_1 = w3_1 .- movavgw3_1

wT1_1 = devT1_1 .* devw1_1
wT2_1 = devT2_1 .* devw2_1
wT3_1 = devT3_1 .* devw3_1

#=
u1_2 = U2[take1_2, skiptoidx2:end, 1]
u2_2 = U2[take2_2, skiptoidx2:end, 1]
u3_2 = U2[take3_2, skiptoidx2:end, 1]
w1_2 = U2[take1_2, skiptoidx2:end, 3]
w2_2 = U2[take2_2, skiptoidx2:end, 3]
w3_2 = U2[take3_2, skiptoidx2:end, 3]
T1_2 = T2[take1_2, skiptoidx2:end]
T2_2 = T2[take2_2, skiptoidx2:end]
T3_2 = T2[take3_2, skiptoidx2:end]
times_points_2 = times2[skiptoidx2:end]
=#
#spectra

#evaluate time steps
timesteps1 = func.gettimesteps(times_points_1)
#timesteps2 = func.gettimesteps(times_points_2)

#choose new time step
timesteptotake1 = 0.008 #0.045 #[s]
#timesteptotake2 = 0.008 #0.045 #[s]

#plot histogram of time steps
fig1 = PyPlot.figure()
ax1 = fig1.add_subplot(111)
ax1.set_title("Time steps taken for following plot")
ax1.set_xlabel("time step [1e-3 s]")
ax1.set_ylabel("density")
ax1.hist(timesteps1 .* 1e3, bins=collect(4:0.05:8), density=true, label="distribution");
ax1.axvline(timesteptotake1 *1e3, 0, 1, color="red", label="time step taken")
ax1.grid(true)
ax1.legend()
#=
fig2 = PyPlot.figure()
ax2 = fig2.add_subplot(111)
ax2.set_title("Time steps for second run")
ax2.set_xlabel("time step [1e-3 s]")
ax2.set_ylabel("density")
ax2.hist(timesteps2 .* 1e3, bins=collect(0:0.05:100), density=true, label="distribution");
ax2.axvline(timesteptotake2 *1e3, 0, 1, color="red", label="time step taken")
ax2.grid(true)
ax2.legend()
=#
println(string("run 1: mean(timesteps) [s] = ", mean(timesteps1), "; median(timesteps) [s] = ", median(timesteps1), ", taken [s] = ", timesteptotake1))
#println(string("run 2: mean(timesteps) [s] = ", mean(timesteps2), "; median(timesteps) [s] = ", median(timesteps2), ", taken [s] = ", timesteptotake2))

#make new time series with equal distance between times
finaltimes_1 = func.createequaltimesteps(timesteptotake1, times_points_1)
u1_1_equal = func.linearmaptotimeseries(u1_1, times_points_1, finaltimes_1)
u2_1_equal = func.linearmaptotimeseries(u2_1, times_points_1, finaltimes_1)
u3_1_equal = func.linearmaptotimeseries(u3_1, times_points_1, finaltimes_1)
w1_1_equal = func.linearmaptotimeseries(w1_1, times_points_1, finaltimes_1)
w2_1_equal = func.linearmaptotimeseries(w2_1, times_points_1, finaltimes_1)
w3_1_equal = func.linearmaptotimeseries(w3_1, times_points_1, finaltimes_1)
T1_1_equal = func.linearmaptotimeseries(T1_1, times_points_1, finaltimes_1)
T2_1_equal = func.linearmaptotimeseries(T2_1, times_points_1, finaltimes_1)
T3_1_equal = func.linearmaptotimeseries(T3_1, times_points_1, finaltimes_1)
wT1_1_equal = func.linearmaptotimeseries(wT1_1, times_points_1, finaltimes_1)
wT2_1_equal = func.linearmaptotimeseries(wT2_1, times_points_1, finaltimes_1)
wT3_1_equal = func.linearmaptotimeseries(wT3_1, times_points_1, finaltimes_1)
#=
finaltimes_2 = func.createequaltimesteps(timesteptotake2, times_points_2)
u1_2_equal = func.linearmaptotimeseries(u1_2, times_points_2, finaltimes_2)
u2_2_equal = func.linearmaptotimeseries(u2_2, times_points_2, finaltimes_2)
u3_2_equal = func.linearmaptotimeseries(u3_2, times_points_2, finaltimes_2)
w1_2_equal = func.linearmaptotimeseries(w1_2, times_points_2, finaltimes_2)
w2_2_equal = func.linearmaptotimeseries(w2_2, times_points_2, finaltimes_2)
w3_2_equal = func.linearmaptotimeseries(w3_2, times_points_2, finaltimes_2)
T1_2_equal = func.linearmaptotimeseries(T1_2, times_points_2, finaltimes_2)
T2_2_equal = func.linearmaptotimeseries(T2_2, times_points_2, finaltimes_2)
T3_2_equal = func.linearmaptotimeseries(T3_2, times_points_2, finaltimes_2)
=#
#code copied from ibl_patch_snow/code/scripts/turb/spectra/fourier_decomp.jl @a107495a
using LsqFit, FFTW, Statistics, Dates
using PyPlot, PyCall
#computer name and adapt file structure to it
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    pathtofile = "/home/michi/Documents/slf/data/"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    pathtofile = "/home/haugened/Documents/data/ir/"
end
include(joinpath(importdir, "src", "fourier.jl"))
import .ft

#data to decompose

#duration in seconds
dur_1 = finaltimes_1[end] - finaltimes_1[1]
#dur_2 = finaltimes_2[end] - finaltimes_2[1]

println("Preprocessing...")
#detrend the vector (subtract linear fit)
u1_1_ftin = ft.belltaper(ft.detrend(u1_1_equal))
u2_1_ftin = ft.belltaper(ft.detrend(u2_1_equal))
u3_1_ftin = ft.belltaper(ft.detrend(u3_1_equal))
w1_1_ftin = ft.belltaper(ft.detrend(w1_1_equal))
w2_1_ftin = ft.belltaper(ft.detrend(w2_1_equal))
w3_1_ftin = ft.belltaper(ft.detrend(w3_1_equal))
T1_1_ftin = ft.belltaper(ft.detrend(T1_1_equal))
T2_1_ftin = ft.belltaper(ft.detrend(T2_1_equal))
T3_1_ftin = ft.belltaper(ft.detrend(T3_1_equal))
wT1_1_ftin = ft.belltaper(ft.detrend(wT1_1_equal))
wT2_1_ftin = ft.belltaper(ft.detrend(wT2_1_equal))
wT3_1_ftin = ft.belltaper(ft.detrend(wT3_1_equal))
#=
u1_2_ftin = ft.belltaper(ft.detrend(u1_2_equal))
u2_2_ftin = ft.belltaper(ft.detrend(u2_2_equal))
u3_2_ftin = ft.belltaper(ft.detrend(u3_2_equal))
w1_2_ftin = ft.belltaper(ft.detrend(w1_2_equal))
w2_2_ftin = ft.belltaper(ft.detrend(w2_2_equal))
w3_2_ftin = ft.belltaper(ft.detrend(w3_2_equal))
T1_2_ftin = ft.belltaper(ft.detrend(T1_2_equal))
T2_2_ftin = ft.belltaper(ft.detrend(T2_2_equal))
T3_2_ftin = ft.belltaper(ft.detrend(T3_2_equal))
=#
println("Applying Fourier-Trafo...")
(a , FTu1_1) = ft.fourierplus(u1_1_ftin)
(a , FTu2_1) = ft.fourierplus(u2_1_ftin)
(a , FTu3_1) = ft.fourierplus(u3_1_ftin)
(a , FTw1_1) = ft.fourierplus(w1_1_ftin)
(a , FTw2_1) = ft.fourierplus(w2_1_ftin)
(a , FTw3_1) = ft.fourierplus(w3_1_ftin)
(a , FTT1_1) = ft.fourierplus(T1_1_ftin)
(a , FTT2_1) = ft.fourierplus(T2_1_ftin)
(a , FTT3_1) = ft.fourierplus(T3_1_ftin)
(a , FTwT1_1) = ft.fourierplus(wT1_1_ftin)
(a , FTwT2_1) = ft.fourierplus(wT2_1_ftin)
(a , FTwT3_1) = ft.fourierplus(wT3_1_ftin)
(Soff_u1_1, frequ1_1) = ft.spectralenergydensity(FTu1_1, dur_1)
(Soff_u2_1, frequ2_1) = ft.spectralenergydensity(FTu2_1, dur_1)
(Soff_u3_1, frequ3_1) = ft.spectralenergydensity(FTu3_1, dur_1)
(Soff_w1_1, freqw1_1) = ft.spectralenergydensity(FTw1_1, dur_1)
(Soff_w2_1, freqw2_1) = ft.spectralenergydensity(FTw2_1, dur_1)
(Soff_w3_1, freqw3_1) = ft.spectralenergydensity(FTw3_1, dur_1)
(Soff_T1_1, freqT1_1) = ft.spectralenergydensity(FTT1_1, dur_1)
(Soff_T2_1, freqT2_1) = ft.spectralenergydensity(FTT2_1, dur_1)
(Soff_T3_1, freqT3_1) = ft.spectralenergydensity(FTT3_1, dur_1)
(Soff_wT1_1, freqwT1_1) = ft.spectralenergydensity(FTwT1_1, dur_1)
(Soff_wT2_1, freqwT2_1) = ft.spectralenergydensity(FTwT2_1, dur_1)
(Soff_wT3_1, freqwT3_1) = ft.spectralenergydensity(FTwT3_1, dur_1)
#=
(a , FTu1_2) = ft.fourierplus(u1_2_ftin)
(a , FTu2_2) = ft.fourierplus(u2_2_ftin)
(a , FTu3_2) = ft.fourierplus(u3_2_ftin)
(a , FTw1_2) = ft.fourierplus(w1_2_ftin)
(a , FTw2_2) = ft.fourierplus(w2_2_ftin)
(a , FTw3_2) = ft.fourierplus(w3_2_ftin)
(a , FTT1_2) = ft.fourierplus(T1_2_ftin)
(a , FTT2_2) = ft.fourierplus(T2_2_ftin)
(a , FTT3_2) = ft.fourierplus(T3_2_ftin)
(Soff_u1_2, frequ1_2) = ft.spectralenergydensity(FTu1_2, dur_2)
(Soff_u2_2, frequ2_2) = ft.spectralenergydensity(FTu2_2, dur_2)
(Soff_u3_2, frequ3_2) = ft.spectralenergydensity(FTu3_2, dur_2)
(Soff_w1_2, freqw1_2) = ft.spectralenergydensity(FTw1_2, dur_2)
(Soff_w2_2, freqw2_2) = ft.spectralenergydensity(FTw2_2, dur_2)
(Soff_w3_2, freqw3_2) = ft.spectralenergydensity(FTw3_2, dur_2)
(Soff_T1_2, freqT1_2) = ft.spectralenergydensity(FTT1_2, dur_2)
(Soff_T2_2, freqT2_2) = ft.spectralenergydensity(FTT2_2, dur_2)
(Soff_T3_2, freqT3_2) = ft.spectralenergydensity(FTT3_2, dur_2)
=#
freqtimesSoff_u1_1 = frequ1_1.*Soff_u1_1
freqtimesSoff_u2_1 = frequ2_1.*Soff_u2_1
freqtimesSoff_u3_1 = frequ3_1.*Soff_u3_1
freqtimesSoff_w1_1 = freqw1_1.*Soff_w1_1
freqtimesSoff_w2_1 = freqw2_1.*Soff_w2_1
freqtimesSoff_w3_1 = freqw3_1.*Soff_w3_1
freqtimesSoff_T1_1 = freqT1_1.*Soff_T1_1
freqtimesSoff_T2_1 = freqT2_1.*Soff_T2_1
freqtimesSoff_T3_1 = freqT3_1.*Soff_T3_1
freqtimesSoff_wT1_1 = freqwT1_1 .* Soff_wT1_1
freqtimesSoff_wT2_1 = freqwT2_1 .* Soff_wT2_1
freqtimesSoff_wT3_1 = freqwT3_1 .* Soff_wT3_1
(freq_u1_1, freqtimesSoff_u1_1) = ft.logavg(frequ1_1, freqtimesSoff_u1_1, 0.05)
(freq_u2_1, freqtimesSoff_u2_1) = ft.logavg(frequ2_1, freqtimesSoff_u2_1, 0.05)
(freq_u3_1, freqtimesSoff_u3_1) = ft.logavg(frequ3_1, freqtimesSoff_u3_1, 0.05)
(freq_w1_1, freqtimesSoff_w1_1) = ft.logavg(freqw1_1, freqtimesSoff_w1_1, 0.05)
(freq_w2_1, freqtimesSoff_w2_1) = ft.logavg(freqw2_1, freqtimesSoff_w2_1, 0.05)
(freq_w3_1, freqtimesSoff_w3_1) = ft.logavg(freqw3_1, freqtimesSoff_w3_1, 0.05)
(freq_T1_1, freqtimesSoff_T1_1) = ft.logavg(freqT1_1, freqtimesSoff_T1_1, 0.05)
(freq_T2_1, freqtimesSoff_T2_1) = ft.logavg(freqT2_1, freqtimesSoff_T2_1, 0.05)
(freq_T3_1, freqtimesSoff_T3_1) = ft.logavg(freqT3_1, freqtimesSoff_T3_1, 0.05)
(freq_wT1_1, freqtimesSoff_wT1_1) = ft.logavg(freqwT1_1, freqtimesSoff_wT1_1, 0.05)
(freq_wT2_1, freqtimesSoff_wT2_1) = ft.logavg(freqwT2_1, freqtimesSoff_wT2_1, 0.05)
(freq_wT3_1, freqtimesSoff_wT3_1) = ft.logavg(freqwT3_1, freqtimesSoff_wT3_1, 0.05)
#=
freqtimesSoff_u1_2 = frequ1_2.*Soff_u1_2
freqtimesSoff_u2_2 = frequ2_2.*Soff_u2_2
freqtimesSoff_u3_2 = frequ3_2.*Soff_u3_2
freqtimesSoff_w1_2 = freqw1_2.*Soff_w1_2
freqtimesSoff_w2_2 = freqw2_2.*Soff_w2_2
freqtimesSoff_w3_2 = freqw3_2.*Soff_w3_2
freqtimesSoff_T1_2 = freqT1_2.*Soff_T1_2
freqtimesSoff_T2_2 = freqT2_2.*Soff_T2_2
freqtimesSoff_T3_2 = freqT3_2.*Soff_T3_2
(freq_u1_2, freqtimesSoff_u1_2) = ft.logavg(frequ1_2, freqtimesSoff_u1_2, 0.05)
(freq_u2_2, freqtimesSoff_u2_2) = ft.logavg(frequ2_2, freqtimesSoff_u2_2, 0.05)
(freq_u3_2, freqtimesSoff_u3_2) = ft.logavg(frequ3_2, freqtimesSoff_u3_2, 0.05)
(freq_w1_2, freqtimesSoff_w1_2) = ft.logavg(freqw1_2, freqtimesSoff_w1_2, 0.05)
(freq_w2_2, freqtimesSoff_w2_2) = ft.logavg(freqw2_2, freqtimesSoff_w2_2, 0.05)
(freq_w3_2, freqtimesSoff_w3_2) = ft.logavg(freqw3_2, freqtimesSoff_w3_2, 0.05)
(freq_T1_2, freqtimesSoff_T1_2) = ft.logavg(freqT1_2, freqtimesSoff_T1_2, 0.05)
(freq_T2_2, freqtimesSoff_T2_2) = ft.logavg(freqT2_2, freqtimesSoff_T2_2, 0.05)
(freq_T3_2, freqtimesSoff_T3_2) = ft.logavg(freqT3_2, freqtimesSoff_T3_2, 0.05)
=#
plot_u1_1 = func.movingaverage(freqtimesSoff_u1_1, 8)
plot_u2_1 = func.movingaverage(freqtimesSoff_u2_1, 8)
plot_u3_1 = func.movingaverage(freqtimesSoff_u3_1, 8)
plot_w1_1 = func.movingaverage(freqtimesSoff_w1_1, 8)
plot_w2_1 = func.movingaverage(freqtimesSoff_w2_1, 8)
plot_w3_1 = func.movingaverage(freqtimesSoff_w3_1, 8)
plot_T1_1 = func.movingaverage(freqtimesSoff_T1_1, 8)
plot_T2_1 = func.movingaverage(freqtimesSoff_T2_1, 8)
plot_T3_1 = func.movingaverage(freqtimesSoff_T3_1, 8)
plot_wT1_1 = func.movingaverage(freqtimesSoff_wT1_1, 8)
plot_wT2_1 = func.movingaverage(freqtimesSoff_wT2_1, 8)
plot_wT3_1 = func.movingaverage(freqtimesSoff_wT3_1, 8)
#=
plot_u1_2 = func.movingaverage(freqtimesSoff_u1_2, 8)
plot_u2_2 = func.movingaverage(freqtimesSoff_u2_2, 8)
plot_u3_2 = func.movingaverage(freqtimesSoff_u3_2, 8)
plot_w1_2 = func.movingaverage(freqtimesSoff_w1_2, 8)
plot_w2_2 = func.movingaverage(freqtimesSoff_w2_2, 8)
plot_w3_2 = func.movingaverage(freqtimesSoff_w3_2, 8)
plot_T1_2 = func.movingaverage(freqtimesSoff_T1_2, 8)
plot_T2_2 = func.movingaverage(freqtimesSoff_T2_2, 8)
plot_T3_2 = func.movingaverage(freqtimesSoff_T3_2, 8)
=#
#sonics
dur_kaijo = Dates.value(Second(kaijodata.time[end] - kaijodata.time[1]))
dur_t2irg = Dates.value(Second(t2irgdata.time[end] - t2irgdata.time[1]))
dur_t1irg = Dates.value(Second(t1irgdata.time[end] - t1irgdata.time[1]))
dur_t2lcsat = Dates.value(Second(t2lcsatdata.time[end] - t2lcsatdata.time[1]))
dur_t2ucsat = Dates.value(Second(t2ucsatdata.time[end] - t2ucsatdata.time[1]))
dur_tjk = Dates.value(Second(tjkdata.time[end] - tjkdata.time[1]))

kaijou_ftin = ft.belltaper(ft.detrend(kaijodata.u))
kaijow_ftin = ft.belltaper(ft.detrend(kaijodata.w))
kaijoT_ftin = ft.belltaper(ft.detrend(kaijodata.T))
kaijowT_ftin = ft.belltaper(ft.detrend(kaijodata.wT))
t2irgu_ftin = ft.belltaper(ft.detrend(t2irgdata.u))
t2irgw_ftin = ft.belltaper(ft.detrend(t2irgdata.w))
t2irgT_ftin = ft.belltaper(ft.detrend(t2irgdata.T))
t2irgwT_ftin = ft.belltaper(ft.detrend(t2irgdata.wT))
t1irgu_ftin = ft.belltaper(ft.detrend(t1irgdata.u))
t1irgw_ftin = ft.belltaper(ft.detrend(t1irgdata.w))
t1irgT_ftin = ft.belltaper(ft.detrend(t1irgdata.T))
t1irgwT_ftin = ft.belltaper(ft.detrend(t1irgdata.wT))
t2lcsatu_ftin = ft.belltaper(ft.detrend(t2lcsatdata.u))
t2lcsatw_ftin = ft.belltaper(ft.detrend(t2lcsatdata.w))
t2lcsatT_ftin = ft.belltaper(ft.detrend(t2lcsatdata.T))
t2lcsatwT_ftin = ft.belltaper(ft.detrend(t2lcsatdata.wT))
t2ucsatu_ftin = ft.belltaper(ft.detrend(t2ucsatdata.u))
t2ucsatw_ftin = ft.belltaper(ft.detrend(t2ucsatdata.w))
t2ucsatT_ftin = ft.belltaper(ft.detrend(t2ucsatdata.T))
t2ucsatwT_ftin = ft.belltaper(ft.detrend(t2ucsatdata.wT))
tjku_ftin = ft.belltaper(ft.detrend(tjkdata.u))
tjkw_ftin = ft.belltaper(ft.detrend(tjkdata.w))
tjkT_ftin = ft.belltaper(ft.detrend(tjkdata.T))
tjkwT_ftin = ft.belltaper(ft.detrend(tjkdata.wT))

(a, FTvec_kaijo_u) = ft.fourierplus(kaijou_ftin)
(a, FTvec_kaijo_w) = ft.fourierplus(kaijow_ftin)
(a, FTvec_kaijo_T) = ft.fourierplus(kaijoT_ftin)
(a, FTvec_kaijo_wT) = ft.fourierplus(kaijowT_ftin)
(a, FTvec_t2irg_u) = ft.fourierplus(t2irgu_ftin)
(a, FTvec_t2irg_w) = ft.fourierplus(t2irgw_ftin)
(a, FTvec_t2irg_T) = ft.fourierplus(t2irgT_ftin)
(a, FTvec_t2irg_wT) = ft.fourierplus(t2irgwT_ftin)
(a, FTvec_t1irg_u) = ft.fourierplus(t1irgu_ftin)
(a, FTvec_t1irg_w) = ft.fourierplus(t1irgw_ftin)
(a, FTvec_t1irg_T) = ft.fourierplus(t1irgT_ftin)
(a, FTvec_t1irg_wT) = ft.fourierplus(t1irgwT_ftin)
(a, FTvec_t2lcsat_u) = ft.fourierplus(t2lcsatu_ftin)
(a, FTvec_t2lcsat_w) = ft.fourierplus(t2lcsatw_ftin)
(a, FTvec_t2lcsat_T) = ft.fourierplus(t2lcsatT_ftin)
(a, FTvec_t2lcsat_wT) = ft.fourierplus(t2lcsatwT_ftin)
(a, FTvec_t2ucsat_u) = ft.fourierplus(t2ucsatu_ftin)
(a, FTvec_t2ucsat_w) = ft.fourierplus(t2ucsatw_ftin)
(a, FTvec_t2ucsat_T) = ft.fourierplus(t2ucsatT_ftin)
(a, FTvec_t2ucsat_wT) = ft.fourierplus(t2ucsatwT_ftin)
(a, FTvec_tjk_u) = ft.fourierplus(tjku_ftin)
(a, FTvec_tjk_w) = ft.fourierplus(tjkw_ftin)
(a, FTvec_tjk_T) = ft.fourierplus(tjkT_ftin)
(a, FTvec_tjk_wT) = ft.fourierplus(tjkwT_ftin)

(Soff_kaijo_u, freq_kaijo_u) = ft.spectralenergydensity(FTvec_kaijo_u, dur_kaijo)
(Soff_kaijo_w, freq_kaijo_w) = ft.spectralenergydensity(FTvec_kaijo_w, dur_kaijo)
(Soff_kaijo_T, freq_kaijo_T) = ft.spectralenergydensity(FTvec_kaijo_T, dur_kaijo)
(Soff_kaijo_wT, freq_kaijo_wT) = ft.spectralenergydensity(FTvec_kaijo_wT, dur_kaijo)
(Soff_t2irg_u, freq_t2irg_u) = ft.spectralenergydensity(FTvec_t2irg_u, dur_t2irg)
(Soff_t2irg_w, freq_t2irg_w) = ft.spectralenergydensity(FTvec_t2irg_w, dur_t2irg)
(Soff_t2irg_T, freq_t2irg_T) = ft.spectralenergydensity(FTvec_t2irg_T, dur_t2irg)
(Soff_t2irg_wT, freq_t2irg_wT) = ft.spectralenergydensity(FTvec_t2irg_wT, dur_t2irg)
(Soff_t1irg_u, freq_t1irg_u) = ft.spectralenergydensity(FTvec_t1irg_u, dur_t1irg)
(Soff_t1irg_w, freq_t1irg_w) = ft.spectralenergydensity(FTvec_t1irg_w, dur_t1irg)
(Soff_t1irg_T, freq_t1irg_T) = ft.spectralenergydensity(FTvec_t1irg_T, dur_t1irg)
(Soff_t1irg_wT, freq_t1irg_wT) = ft.spectralenergydensity(FTvec_t1irg_wT, dur_t1irg)
(Soff_t2lcsat_u, freq_t2lcsat_u) = ft.spectralenergydensity(FTvec_t2lcsat_u, dur_t2lcsat)
(Soff_t2lcsat_w, freq_t2lcsat_w) = ft.spectralenergydensity(FTvec_t2lcsat_w, dur_t2lcsat)
(Soff_t2lcsat_T, freq_t2lcsat_T) = ft.spectralenergydensity(FTvec_t2lcsat_T, dur_t2lcsat)
(Soff_t2lcsat_wT, freq_t2lcsat_wT) = ft.spectralenergydensity(FTvec_t2lcsat_wT, dur_t2lcsat)
(Soff_t2ucsat_u, freq_t2ucsat_u) = ft.spectralenergydensity(FTvec_t2ucsat_u, dur_t2ucsat)
(Soff_t2ucsat_w, freq_t2ucsat_w) = ft.spectralenergydensity(FTvec_t2ucsat_w, dur_t2ucsat)
(Soff_t2ucsat_T, freq_t2ucsat_T) = ft.spectralenergydensity(FTvec_t2ucsat_T, dur_t2ucsat)
(Soff_t2ucsat_wT, freq_t2ucsat_wT) = ft.spectralenergydensity(FTvec_t2ucsat_wT, dur_t2ucsat)
(Soff_tjk_u, freq_tjk_u) = ft.spectralenergydensity(FTvec_tjk_u, dur_tjk)
(Soff_tjk_w, freq_tjk_w) = ft.spectralenergydensity(FTvec_tjk_w, dur_tjk)
(Soff_tjk_T, freq_tjk_T) = ft.spectralenergydensity(FTvec_tjk_T, dur_tjk)
(Soff_tjk_wT, freq_tjk_wT) = ft.spectralenergydensity(FTvec_tjk_wT, dur_tjk)

freqtimesSoff_kaijo_u = freq_kaijo_u .* Soff_kaijo_u
freqtimesSoff_kaijo_w = freq_kaijo_w .* Soff_kaijo_w
freqtimesSoff_kaijo_T = freq_kaijo_T .* Soff_kaijo_T
freqtimesSoff_kaijo_wT = freq_kaijo_wT .* Soff_kaijo_wT
freqtimesSoff_t2irg_u = freq_t2irg_u .* Soff_t2irg_u
freqtimesSoff_t2irg_w = freq_t2irg_w .* Soff_t2irg_w
freqtimesSoff_t2irg_T = freq_t2irg_T .* Soff_t2irg_T
freqtimesSoff_t2irg_wT = freq_t2irg_wT .* Soff_t2irg_wT
freqtimesSoff_t1irg_u = freq_t1irg_u .* Soff_t1irg_u
freqtimesSoff_t1irg_w = freq_t1irg_w .* Soff_t1irg_w
freqtimesSoff_t1irg_T = freq_t1irg_T .* Soff_t1irg_T
freqtimesSoff_t1irg_wT = freq_t1irg_wT .* Soff_t1irg_wT
freqtimesSoff_t2lcsat_u = freq_t2lcsat_u .* Soff_t2lcsat_u
freqtimesSoff_t2lcsat_w = freq_t2lcsat_w .* Soff_t2lcsat_w
freqtimesSoff_t2lcsat_T = freq_t2lcsat_T .* Soff_t2lcsat_T
freqtimesSoff_t2lcsat_wT = freq_t2lcsat_wT .* Soff_t2lcsat_wT
freqtimesSoff_t2ucsat_u = freq_t2ucsat_u .* Soff_t2ucsat_u
freqtimesSoff_t2ucsat_w = freq_t2ucsat_w .* Soff_t2ucsat_w
freqtimesSoff_t2ucsat_T = freq_t2ucsat_T .* Soff_t2ucsat_T
freqtimesSoff_t2ucsat_wT = freq_t2ucsat_wT .* Soff_t2ucsat_wT
freqtimesSoff_tjk_u = freq_tjk_u .* Soff_tjk_u
freqtimesSoff_tjk_w = freq_tjk_w .* Soff_tjk_w
freqtimesSoff_tjk_T = freq_tjk_T .* Soff_tjk_T
freqtimesSoff_tjk_wT = freq_tjk_wT .* Soff_tjk_wT

(freq_kaijo_u, freqtimesSoff_kaijo_u) = ft.logavg(freq_kaijo_u, freqtimesSoff_kaijo_u, 0.05)
(freq_kaijo_w, freqtimesSoff_kaijo_w) = ft.logavg(freq_kaijo_w, freqtimesSoff_kaijo_w, 0.05)
(freq_kaijo_T, freqtimesSoff_kaijo_T) = ft.logavg(freq_kaijo_T, freqtimesSoff_kaijo_T, 0.05)
(freq_kaijo_wT, freqtimesSoff_kaijo_wT) = ft.logavg(freq_kaijo_wT, freqtimesSoff_kaijo_wT, 0.05)
(freq_t2irg_u, freqtimesSoff_t2irg_u) = ft.logavg(freq_t2irg_u, freqtimesSoff_t2irg_u, 0.05)
(freq_t2irg_w, freqtimesSoff_t2irg_w) = ft.logavg(freq_t2irg_w, freqtimesSoff_t2irg_w, 0.05)
(freq_t2irg_T, freqtimesSoff_t2irg_T) = ft.logavg(freq_t2irg_T, freqtimesSoff_t2irg_T, 0.05)
(freq_t2irg_wT, freqtimesSoff_t2irg_wT) = ft.logavg(freq_t2irg_wT, freqtimesSoff_t2irg_wT, 0.05)
(freq_t1irg_u, freqtimesSoff_t1irg_u) = ft.logavg(freq_t1irg_u, freqtimesSoff_t1irg_u, 0.05)
(freq_t1irg_w, freqtimesSoff_t1irg_w) = ft.logavg(freq_t1irg_w, freqtimesSoff_t1irg_w, 0.05)
(freq_t1irg_T, freqtimesSoff_t1irg_T) = ft.logavg(freq_t1irg_T, freqtimesSoff_t1irg_T, 0.05)
(freq_t1irg_wT, freqtimesSoff_t1irg_wT) = ft.logavg(freq_t1irg_wT, freqtimesSoff_t1irg_wT, 0.05)
(freq_t2lcsat_u, freqtimesSoff_t2lcsat_u) = ft.logavg(freq_t2lcsat_u, freqtimesSoff_t2lcsat_u, 0.05)
(freq_t2lcsat_w, freqtimesSoff_t2lcsat_w) = ft.logavg(freq_t2lcsat_w, freqtimesSoff_t2lcsat_w, 0.05)
(freq_t2lcsat_T, freqtimesSoff_t2lcsat_T) = ft.logavg(freq_t2lcsat_T, freqtimesSoff_t2lcsat_T, 0.05)
(freq_t2lcsat_wT, freqtimesSoff_t2lcsat_wT) = ft.logavg(freq_t2lcsat_wT, freqtimesSoff_t2lcsat_wT, 0.05)
(freq_t2ucsat_u, freqtimesSoff_t2ucsat_u) = ft.logavg(freq_t2ucsat_u, freqtimesSoff_t2ucsat_u, 0.05)
(freq_t2ucsat_w, freqtimesSoff_t2ucsat_w) = ft.logavg(freq_t2ucsat_w, freqtimesSoff_t2ucsat_w, 0.05)
(freq_t2ucsat_T, freqtimesSoff_t2ucsat_T) = ft.logavg(freq_t2ucsat_T, freqtimesSoff_t2ucsat_T, 0.05)
(freq_t2ucsat_wT, freqtimesSoff_t2ucsat_wT) = ft.logavg(freq_t2ucsat_wT, freqtimesSoff_t2ucsat_wT, 0.05)
(freq_tjk_u, freqtimesSoff_tjk_u) = ft.logavg(freq_tjk_u, freqtimesSoff_tjk_u, 0.05)
(freq_tjk_w, freqtimesSoff_tjk_w) = ft.logavg(freq_tjk_w, freqtimesSoff_tjk_w, 0.05)
(freq_tjk_T, freqtimesSoff_tjk_T) = ft.logavg(freq_tjk_T, freqtimesSoff_tjk_T, 0.05)
(freq_tjk_wT, freqtimesSoff_tjk_wT) = ft.logavg(freq_tjk_wT, freqtimesSoff_tjk_wT, 0.05)

plot_kaijo_u = func.movingaverage(freqtimesSoff_kaijo_u, 8)
plot_kaijo_w = func.movingaverage(freqtimesSoff_kaijo_w, 8)
plot_kaijo_T = func.movingaverage(freqtimesSoff_kaijo_T, 8)
plot_kaijo_wT = func.movingaverage(freqtimesSoff_kaijo_wT, 8)
plot_t2irg_u = func.movingaverage(freqtimesSoff_t2irg_u, 8)
plot_t2irg_w = func.movingaverage(freqtimesSoff_t2irg_w, 8)
plot_t2irg_T = func.movingaverage(freqtimesSoff_t2irg_T, 8)
plot_t2irg_wT = func.movingaverage(freqtimesSoff_t2irg_wT, 8)
plot_t1irg_u = func.movingaverage(freqtimesSoff_t1irg_u, 8)
plot_t1irg_w = func.movingaverage(freqtimesSoff_t1irg_w, 8)
plot_t1irg_T = func.movingaverage(freqtimesSoff_t1irg_T, 8)
plot_t1irg_wT = func.movingaverage(freqtimesSoff_t1irg_wT, 8)
plot_t2lcsat_u = func.movingaverage(freqtimesSoff_t2lcsat_u, 8)
plot_t2lcsat_w = func.movingaverage(freqtimesSoff_t2lcsat_w, 8)
plot_t2lcsat_T = func.movingaverage(freqtimesSoff_t2lcsat_T, 8)
plot_t2lcsat_wT = func.movingaverage(freqtimesSoff_t2lcsat_wT, 8)
plot_t2ucsat_u = func.movingaverage(freqtimesSoff_t2ucsat_u, 8)
plot_t2ucsat_w = func.movingaverage(freqtimesSoff_t2ucsat_w, 8)
plot_t2ucsat_T = func.movingaverage(freqtimesSoff_t2ucsat_T, 8)
plot_t2ucsat_wT = func.movingaverage(freqtimesSoff_t2ucsat_wT, 8)
plot_tjk_u = func.movingaverage(freqtimesSoff_tjk_u, 8)
plot_tjk_w = func.movingaverage(freqtimesSoff_tjk_w, 8)
plot_tjk_T = func.movingaverage(freqtimesSoff_tjk_T, 8)
plot_tjk_wT = func.movingaverage(freqtimesSoff_tjk_wT, 8)

##
#PyPlot.pygui(true)
fig2 = PyPlot.figure(figsize=(12,11))
gs = gridspec.GridSpec(4,2, height_ratios=[10,10,10,2])
for i in 1:gs.nrows-1
    #time series
    ax1 = fig2.add_subplot(gs[i,1])
    plotx = times_points_1
    if i == 1
        ploty1 = u1_1
        ploty2 = u2_1
        ploty3 = u3_1
        sonic_symbol = :u
        textstr = "a)"
        ax1.text(0.41, 1.1, "time series", size="large", transform=ax1.transAxes)
    elseif i == 2
        ploty1 = w1_1
        ploty2 = w2_1
        ploty3 = w3_1
        sonic_symbol = :w
        textstr = "c)"
    elseif i==3
        ploty1 = T1_1
        ploty2 = T2_1
        ploty3 = T3_1
        sonic_symbol = :T
        textstr = "e)"
    elseif i==4
        ploty1 = wT1_1
        ploty2 = wT2_1
        ploty3 = wT3_1
        textstr = "g)"
        sonic_symbol = :wT
    end

    #sonics
    #ax1.plot(1/1000 .* Dates.value.(Millisecond.(kaijodata.time .- sonicsstart)), kaijodata[:,sonic_symbol])
    #ax1.plot(1/1000 .* Dates.value.(Millisecond.(t2irgdata.time .- sonicsstart)), t2irgdata[:,sonic_symbol])
    ax1.plot(1/1000 .* Dates.value.(Millisecond.(t1irgdata.time .- sonicsstart)), t1irgdata[:,sonic_symbol])
    #ax1.plot(1/1000 .* Dates.value.(Millisecond.(t2lcsatdata.time .- sonicsstart)), t2lcsatdata[:,sonic_symbol])
    #ax1.plot(1/1000 .* Dates.value.(Millisecond.(t2ucsatdata.time .- sonicsstart)), t2ucsatdata[:,sonic_symbol])
    #ax1.plot(1/1000 .* Dates.value.(Millisecond.(tjkdata.time .- sonicsstart)), tjkdata[:,sonic_symbol])

    #ax1.plot(plotx, ploty1, "--")#, alpha=0.5)
    #ax1.plot(plotx, ploty2, "--")#, alpha=0.5)
    ax1.plot(plotx, ploty3, "--")#, alpha=0.5)
    
    ax1.text(-0.1, 1.01, textstr, size="large", transform=ax1.transAxes)
    ax1.grid(true)
    ax1.set_xlim(0,300)

    if i==1
        ax1.set_ylabel(L"u~\mathrm{[m~s^{-1}]}")
        ax1.set_ylim(0,6)
    elseif i==2
        ax1.set_ylabel(L"w~\mathrm{[m~s^{-1}]}")
        ax1.set_ylim(-1.4, 1.4)
    elseif i==3
        ax1.set_ylabel(L"T~\mathrm{[^\circ C]}")
        ax1.set_ylim(7, 11.5)
    elseif i==4
        ax1.set_ylabel(L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}")
        ax1.set_ylim(-0.6, 0.4)
    end
    if i==gs.nrows-1
        ax1.set_xlabel(L"t~\mathrm{[s]}")
    else
        #ax1.tick_params(axis="x", bottom=false)
        ax1.tick_params(axis="x", labelbottom=false)
    end
end
for i in 1:gs.nrows-1
    #spectra
    ax = fig2.add_subplot(gs[i,2])
    if i == 1
        plotx1 = freq_u1_1
        plotx2 = freq_u2_1
        plotx3 = freq_u3_1
        ploty1 = plot_u1_1
        ploty2 = plot_u2_1
        ploty3 = plot_u3_1
        kaijox = freq_kaijo_u
        kaijoy = plot_kaijo_u
        t2irgx = freq_t2irg_u
        t2irgy = plot_t2irg_u
        t1irgx = freq_t1irg_u
        t1irgy = plot_t1irg_u
        t2lcsatx = freq_t2lcsat_u
        t2lcsaty = plot_t2lcsat_u
        t2ucsatx = freq_t2ucsat_u
        t2ucsaty = plot_t2ucsat_u
        tjkx = freq_tjk_u
        tjky = plot_tjk_u
        textstr = "b)"
        ax.text(0.43, 1.1, "spectra", size="large", transform=ax.transAxes)
    elseif i == 2
        plotx1 = freq_w1_1
        plotx2 = freq_w2_1
        plotx3 = freq_w3_1
        ploty1 = plot_w1_1
        ploty2 = plot_w2_1
        ploty3 = plot_w3_1
        kaijox = freq_kaijo_w
        kaijoy = plot_kaijo_w
        t2irgx = freq_t2irg_w
        t2irgy = plot_t2irg_w
        t1irgx = freq_t1irg_w
        t1irgy = plot_t1irg_w
        t2lcsatx = freq_t2lcsat_w
        t2lcsaty = plot_t2lcsat_w
        t2ucsatx = freq_t2ucsat_w
        t2ucsaty = plot_t2ucsat_w
        tjkx = freq_tjk_w
        tjky = plot_tjk_w
        textstr = "d)"
    elseif i==3
        plotx1 = freq_T1_1
        plotx2 = freq_T2_1
        plotx3 = freq_T3_1
        ploty1 = plot_T1_1
        ploty2 = plot_T2_1
        ploty3 = plot_T3_1
        kaijox = freq_kaijo_T
        kaijoy = plot_kaijo_T
        t2irgx = freq_t2irg_T
        t2irgy = plot_t2irg_T
        t1irgx = freq_t1irg_T
        t1irgy = plot_t1irg_T
        t2lcsatx = freq_t2lcsat_T
        t2lcsaty = plot_t2lcsat_T
        t2ucsatx = freq_t2ucsat_T
        t2ucsaty = plot_t2ucsat_T
        tjkx = freq_tjk_T
        tjky = plot_tjk_T
        textstr = "f)"
    elseif i==4
        plotx1 = freq_wT1_1
        plotx2 = freq_wT2_1
        plotx3 = freq_wT3_1
        ploty1 = plot_wT1_1
        ploty2 = plot_wT2_1
        ploty3 = plot_wT3_1
        kaijox = freq_kaijo_wT
        kaijoy = plot_kaijo_wT
        t2irgx = freq_t2irg_wT
        t2irgy = plot_t2irg_wT
        t1irgx = freq_t1irg_wT
        t1irgy = plot_t1irg_wT
        t2lcsatx = freq_t2lcsat_wT
        t2lcsaty = plot_t2lcsat_wT
        t2ucsatx = freq_t2ucsat_wT
        t2ucsaty = plot_t2ucsat_wT
        tjkx = freq_tjk_wT
        tjky = plot_tjk_wT
        textstr = "h)"
    end

    #sonics
    #ax.plot(kaijox, kaijoy, label="Kaijo")
    #ax.plot(t2irgx, t2irgy, label="T2IRG")
    ax.plot(t1irgx, t1irgy, label="measurement @ 1.2m")
    #ax.plot(t2lcsatx, t2lcsaty, label="measurement @ h=1.9m")
    #ax.plot(t2ucsatx, t2ucsaty, label="T2UCSAT")
    #ax.plot(tjkx, tjky, label="TJK")

    #ax.plot(plotx1, ploty1, "--", label=string("x=", round(x1[take1_1], digits=2), "m , h=", round(z1[take1_1], digits=2), "m"))
    #ax.plot(plotx2, ploty2, "--",label=string("x=", round(x1[take2_1], digits=2), "m , h=", round(z1[take2_1], digits=2), "m"))
    ax.plot(plotx3, ploty3, "--", label=string("x=", round(x1[take3_1], digits=2), "m , h=", round(z1[take3_1], digits=1), "m"))

    ax.text(-0.13, 1.01, textstr, size="large", transform=ax.transAxes)

    #xtheory=collect(2:0.05:8)
    #ytheory=exp.(-2*log.(xtheory)/3)./15
    #ax.plot(xtheory,ytheory, color="black", label="e^{-2/3}")
    PyPlot.xscale("log")
    PyPlot.yscale("log")
    ax.set_xlim([0.002, 10])
    ax.grid()
    #ax.legend()

    if i == 1
        ax.set_ylim(0.001, 0.3)
    elseif i == 2
        ax.set_ylim(1e-8, 5e-2)
    elseif i == 3
        ax.set_ylim(1e-6, 1e-1)
    elseif i == 4
        ax.set_ylim(1e-7, 5e-3)
    end

    if i == 1 || i == 2
        ax.set_ylabel(L"fS(f)~[\mathrm{m^2~s^{-2}}]")
    elseif i == 3
        ax.set_ylabel(L"fS(f)~[\mathrm{K^2}]")
    elseif i == 4
        ax.set_ylabel(L"fS(f)~[\mathrm{K^2~m^2~s^{-2}}]")
    end
    if i==gs.nrows-1
        ax.set_xlabel(L"f~\mathrm{[Hz]}")
        global (h, l) = ax.get_legend_handles_labels()
    else
        ax.tick_params(axis="x", labelbottom=false)
    end
end
ax_leg = fig2.add_subplot(py"$(gs)[3,0:2]")
ax_leg.legend(h, l, loc="upper center")#,borderaxespad=0)
ax_leg.axis("off")
#PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/cmp_meas_les_3rd.jpg")
##
=#
####################################################################################
####################################################################################
###  Fig. 4:          SNAPSHOTS SHOWING INTERMITTENT ADVECTION                   ###
####################################################################################
#####################################################################################figure comparing intermittent advection in IR-screen and simulation
#=
using PyCall, Statistics, LaTeXStrings, NCDatasets, DataFrames, StatsBase, PyCall, Dates, ProgressMeter, CSV
import PyPlot
gridspec = pyimport("matplotlib.gridspec")
mpimg = pyimport("matplotlib.image")
animation = pyimport("matplotlib.animation")
numpy = pyimport("numpy")
cramericm = pyimport("cmcrameri.cm")
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    pathtofile = "/home/michi/Documents/slf/data/"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    pathtofile = "/home/haugened/Documents/data/ir/"
end
include(joinpath(importdir, "src", "ir_evaluation.jl"))
import .irev

#outer directory of case
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/third_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#folder containing the time steps
readfromfolder = joinpath(casedir, "postProcessing/surfaces/yNormal")

#netCDF4-file containing the profile (length 3); use same time steps!!
file = joinpath.(readfromfolder, "yNormal_time5th.nc")

#value of x where bare to snow transition is
snowtransition = 8.0 #m

#quantity for profiles 
quantity_to_plot = "T"

#read in data
if !@isdefined data_read
    data_read = false
end
if !data_read
    println("Reading cross section data...")
    (x, y, z, times, T, U) = func.readaerialprofilefromnetcdf(file)
    data_read = true
end

#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx
    transformx = false
end
if !transformx
    x = x .- snowtransition
    transformx = true
end

#create bottom line for plot
xtopo = func.createbottom(x, z)

#################################################
#take and/or quantity to plot
vals = T
vals .-= 273.15 #K to C
valmin = 9
valmax = 13
colmap = "viridis"
cblabel = L"T~\mathrm{[^\circ C]}"

#functions
function loadfromNetCDF4(file::String, varname::String)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    dattot = dat[:,:,2900:end-300]
    close(netfile)
    return dattot
end

#definitions
irfile = "/home/haugened/Documents/data/ir/210531_115756/irdata_03_to_05.nc"
cmperpxl = 0.6
cmperpxlu = 0.606
cmperpxlw = 0.549
framespersec = 29.9

IRdata = loadfromNetCDF4(irfile, "irdata")
#IRdata = IRdata[:,:,2000:end]

#cut
IRdata[:,443:540,:] .= NaN
IRdata = IRdata[11:end,19:end,:]
(IRdata, snowsurf) = irev.setsnowsurfacetonan(IRdata, -0.2, framewise=true)

#load setup image
img_setup = mpimg.imread("/home/haugened/Documents/process_paper_23/figures/figures_new/ir_excerpt.png")

##
fntsze = 16
lblsze = fntsze-2
xis0pxl = 5
PyPlot.pygui(true)
pubfig = PyPlot.figure(figsize=(11.5,8))
pubfig.text(0.07, 0.89, "a)", size="large", )
pubfig.text(0.07, 0.65, "b)", size="large", )
pubfig.text(0.5, 0.65, "c)", size="large", )
pubfig.text(0.07, 0.39, "d)", size="large", )
#pubfig.text(0.26, 0.415, L"\Delta t_{IR} = 7.1\mathrm{s}", size="large")
#pubfig.text(0.66, 0.415, L"\Delta t_{LES} = 3.2\mathrm{s}", size="large")
pubfig.text(0.5, 0.39, "e)", size="large", )
gs = gridspec.GridSpec(5, 3, width_ratios=[30,1,30], height_ratios=[8,10,10,0.1,0.8], wspace=0.1)
ax0 = pubfig.add_subplot(py"$(gs)[0,0]")
PyPlot.imshow(img_setup)
PyPlot.axis("off")
ax1 = pubfig.add_subplot(gs[2,1])
#ax1.set_title(L"t_{IR, 1}=12:01:35")
xleft = -(xis0pxl-1)*cmperpxlu/100
xright = size(IRdata, 2)*cmperpxlu/100+xleft
ybottom = 0
ytop = size(IRdata, 1) * cmperpxlw/100
hm1 = ax1.imshow(IRdata[:,:,2966-2900], extent = [xleft,xright,ybottom,ytop], cmap=cramericm.batlow, vmin=4, vmax=8)
ax1.set_ylabel(L"h~\mathrm{[m]}")#, fontsize=fntsze)
#ax1.tick_params(labelsize=lblsze)
ax1.tick_params(axis="x", labelbottom=false)
ax1.set_ylim(0,2)
IRwindarrow = ax1.annotate("", xy=(0.8,2.15), xytext=(0,2.15), annotation_clip=false, arrowprops=Dict("lw" => 2, "arrowstyle" => "->", "facecolor" => "black"))
IRwindarrowtxt = ax1.text(0, 2.21, "wind")
ax3 = pubfig.add_subplot(gs[5,1])
cbar = pubfig.colorbar(hm1, cax=ax3, orientation="horizontal")
cbar.set_label(L"T~\mathrm{[^\circ C]}")#, fontsize=fntsze)
#cbar.ax.tick_params(labelsize=lblsze)
ax2 = pubfig.add_subplot(py"$(gs)[2,0]", sharex=ax1)
#ax2.set_title(L"t_{IR, 1} + 7.1\mathrm{s}")
hm2 = ax2.imshow(IRdata[:,:,3180-2900], extent = [xleft,xright,ybottom,ytop], cmap=cramericm.batlow, vmin=4, vmax=8)
ax2.set_xlabel(L"x~\mathrm{[m]}")#, fontsize=fntsze)
ax2.set_ylabel(L"h~\mathrm{[m]}")#, fontsize=fntsze)
#ax2.tick_params(labelsize=lblsze)
IRwindarrow2 = ax2.annotate("", xy=(0.8,2.15), xytext=(0,2.15), annotation_clip=false, arrowprops=Dict("lw" => 2, "arrowstyle" => "->", "facecolor" => "black"))
IRwindarrowtxt2 = ax2.text(0, 2.21, "wind")
ax4 = pubfig.add_subplot(gs[2,3], sharey=ax1)
PyPlot.gca().set_aspect("equal")
#ax4.set_title(L"t_{LES, 1}=292.2~\mathrm{s}")
hm_ani1 = ax4.tripcolor(x, z, vals[:, 722], cmap=colmap, vmin=valmin, vmax=valmax, shading="gouraud")
topo = ax4.plot(xtopo.x, xtopo.h, color="white", lw=1)
ax4.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
#ax4.set_xlim(minimum(x), maximum(x))
#ax4.set_ylim(0, maximum(z))
#ax4.set_xlabel(L"x~[\mathrm{m}]")
#ax4.set_ylabel(L"h~[\mathrm{m}]")
#ax4.set_xlim(-0.5, 6)
#ax4.set_ylim(0, 2)
ax4.tick_params(axis="x", labelbottom=false)
#hm_cb = fig.colorbar(hm, ax=ax1)
#hm_cb.set_label(cblabel)
aniwindarrow = ax4.annotate("", xy=(0.4,2.15), xytext=(-0.4,2.15), annotation_clip=false, arrowprops=Dict("lw" => 2, "arrowstyle" => "->", "facecolor" => "black"))
aniwindarrowtxt = ax4.text(-0.4, 2.21, "wind")
ax5 = pubfig.add_subplot(gs[3,3], sharey=ax2, sharex=ax4)
#ax5.set_title(L"t_{LES, 1} - 12.24~\mathrm{s}")
hm_ani2 = ax5.tripcolor(x, z, vals[:, 642], cmap=colmap, vmin=valmin, vmax=valmax, shading="gouraud")
topo = ax5.plot(xtopo.x, xtopo.h, color="white", lw=1)
ax5.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
ax5.set_xlim(minimum(x), maximum(x))
ax5.set_ylim(0, maximum(z))
ax5.set_xlabel(L"x~[\mathrm{m}]")
#ax5.set_ylabel(L"h~[\mathrm{m}]")
ax5.set_xlim(-0.5, 6)
ax5.set_ylim(0, 2)
#hm_cb = fig.colorbar(hm, ax=ax1)
#hm_cb.set_label(cblabel)
PyPlot.gca().set_aspect("equal")
aniwindarrow2 = ax5.annotate("", xy=(0.4,2.15), xytext=(-0.4,2.15), annotation_clip=false, arrowprops=Dict("lw" => 2, "arrowstyle" => "->", "facecolor" => "black"))
aniwindarrowtxt2 = ax5.text(-0.4, 2.21, "wind")
ax6 = pubfig.add_subplot(gs[5,3])
cbar = pubfig.colorbar(hm_ani2, cax=ax6, orientation="horizontal")
cbar.set_label(cblabel)#, fontsize=fntsze)
#PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/cmp_IR_LES_3rd.jpg")
##
=#

####################################################################################
####################################################################################
###  Fig. 5:          MEDIAN BUOYANCY FLUXES AND SIBL GROWTH                     ###
####################################################################################
#####################################################################################SIBL depth growth
#=
using NCDatasets, DataFrames, StatsBase, LaTeXStrings, PyCall, Dates, LsqFit, CSV
import PyPlot
PyPlot.pygui(true)
GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")
animation = pyimport("matplotlib.animation")
numpy = pyimport("numpy")
cramericm = pyimport("cmcrameri.cm")

#outer directory of case
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/third_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#folder containing the time steps
readfromfolder = joinpath(casedir, "postProcessing/surfaces/yNormal")

#netCDF4-file containing the profile (length 3); use same time steps!!
file = joinpath.(readfromfolder, "yNormal_time5th.nc")

#value of x where bare to snow transition is
snowtransition = 8.0 #m

reyavgtime = Second(50)

#seconds to skip (spinuptime)
skipsec = 10

#flag for computing and plotting SIBL depth
calcsibl = true
#distance between two neighboring profiles calculating SIBL
calcevery = 0.1 #m
calchalfwidth = 0.03 #m

#quantity used to diagnose SIBL depth
diagnosefromquantity = "wT"
allowedquantitiestodiagnosefrom = ["wT", "z/L"]

#flag for filtering of profiles (see below in code to specify filter condition)
postfilter = false
postfilterquantity = "u" #u, v, w, or T
#location of point to evaluate for post-filtering
postfilterloc = [5, 2, 0.5]
#filter condition
filtercond = "filterdata .>= filterquantiles[80]"

#Interval [ms] for plotting.
#intervalmsec = 30

#calculated for video
#framespersec = 1000/intervalmsec

#file name of video
#anistring = string(joinpath(casedir, "video.mp4"))
#################################################
#input check
if !(diagnosefromquantity in allowedquantitiestodiagnosefrom)
    println(string(diagnosefromquantity), " can not be used to diagnose SIBL depth. Possible values are ", allowedquantities, ". Exiting.")
    exit()
else
    println("Diagnosing SIBL depth from ", diagnosefromquantity)
end

#################################################
#read in data
if !@isdefined data_read
    data_read = false
end
if !data_read
    println("Reading cross section data...")
    (x, y, z, times, T, U) = func.readaerialprofilefromnetcdf(file)
    data_read = true
end

#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx
    transformx = false
end
if !transformx
    x = x .- snowtransition
    transformx = true
end

#create bottom line for plot
xtopo = func.createbottom(x, z)

#################################################
if diagnosefromquantity == "wT"
    plotlim = 10
    cblabel = L"\overline{w'T'}~\mathrm{[mK~m~s^{-1}]}"
    if postfilter
        figtitle = string("Median w'T' skipping first ", skipsec, "s; ", filtercond)
    else
        figtitle = string("Median w'T' skipping first ", skipsec, "s; No filtering")
    end

    #get stuff for evaluating wT
    vals_w = U[:,:,3] .* 1000
    vals_T = T

    #calculate fluxes
    quant_crosssection = similar(vals_T)
    fill!(quant_crosssection, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))

    for jrow in 1:size(vals_w, 1)
            movavg_w = func.movingaverage(vals_w[jrow, :], movavgidcs)
            movavg_T = func.movingaverage(vals_T[jrow, :], movavgidcs)
            dev1 = vals_w[jrow, :] .- movavg_w
            dev2 = vals_T[jrow, :] .- movavg_T
            quant_crosssection[jrow, :] = dev1 .* dev2
    end
elseif diagnosefromquantity == "z/L"
    plotlim = 0.04
    cblabel = L"\frac{z}{L}~\mathrm{[1]}"
    if postfilter
        figtitle = string("Median z/L skipping first ", skipsec, "s; ", filtercond)
    else
        figtitle = string("Median z/L skipping first ", skipsec, "s; No filtering")
    end


    valsu = U[:,:,1]
    valsv = U[:,:,2]
    valsw = U[:,:,3]
    valsT = T

    #calculate fluxes
    quant_crosssection = similar(valsT)
    fill!(quant_crosssection, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))

    for jrow in 1:size(valsw, 1)
            movavg_u = func.movingaverage(valsu[jrow, :], movavgidcs)
            movavg_v = func.movingaverage(valsv[jrow, :], movavgidcs)
            movavg_w = func.movingaverage(valsw[jrow, :], movavgidcs)
            movavg_T = func.movingaverage(valsT[jrow, :], movavgidcs)
            dev_u = valsu[jrow, :] .- movavg_u
            dev_v = valsv[jrow, :] .- movavg_v
            dev_w = valsw[jrow, :] .- movavg_w
            dev_T = valsT[jrow, :] .- movavg_T
            ustar = ((dev_u .* dev_w) .^2 .+ (dev_v .* dev_w) .^2) .^(1/4)
            wT_tmp = dev_w .* dev_T
            L_tmp = - (ustar .^ 3 .* movavg_T) ./ (0.4*9.81 .* wT_tmp)
            quant_crosssection[jrow, :] = z[jrow] ./ L_tmp
    end
end

#timesround = round.(times, digits=3)
#animidcs = func.getanimidcs(timesround, intervalmsec)

#get rid of NaN column
x_anim = x

quant_median = fill(NaN, size(quant_crosssection, 1))

if postfilter
    println("-------------------------------------")
    println("Filtering of profiles activated.")
    println("No filtering of sonic measurement data!")
    println("Filter condition: ", filtercond)

    #get data at filter location
    post_take = func.chooseval3D(postfilterloc, hcat(x, y, z))

    println("location for filtering:")
    println(string("x = ", x[post_take], "m"))
    println(string("y = ", y[post_take], "m"))
    println(string("z = ", z[post_take], "m"))
    println("-------------------------------------")


    if postfilterquantity == "u"
        filterdata = U[post_take,:,1]
    elseif postfilterquantity == "v"
        filterdata = U[post_take,:,2]
    elseif postfilterquantity == "w"
        filterdata = U[post_take,:,3]
    elseif postfilterquantity == "T"
        filterdata = T[post_take,:]
    else
        @error("Wrong filter quantity! Please check.")
    end

    #calculate 1%-quantiles
    quantiles = collect(0.01:0.01:1)
    filterquantiles = quantile(filter(!isnan, filterdata), quantiles)

    #post filter condition
    postfiltercond = eval(Meta.parse(filtercond))
else #no postfiltering
    postfiltercond = fill(true, size(quant_crosssection, 2))
    println("No post-filtering applied.")
end

quant_crosssection[:, .!(postfiltercond)] .= NaN

skiptoidx = findfirst(x->x==minimum(abs.(times .- skipsec)), abs.(times .- skipsec))
println(string("Skipping ", skipsec, " seconds at the beginning for calculating median."))

for jrow in 1:size(quant_median, 1)
    quant_median[jrow] = median(filter(!isnan, quant_crosssection[jrow, skiptoidx:end]))
end

if calcsibl
    #diagnose SIBL depth according to w'T'
    println("Diagnosing SIBL depth from ", diagnosefromquantity)

    sibllocs = collect(minimum(x):calcevery:maximum(x))
    sibldepth_idx = zeros(Int, size(sibllocs, 1))
    sibldepth = zeros(Float64, size(sibllocs, 1))
    sibldepth2 = zeros(Float64, size(sibllocs, 1))

    for iloc in 1:size(sibllocs, 1)
        sortperms = fill(-1, size(x, 1))
        profilelines = func.chooselistval("z", sibllocs[iloc]-calchalfwidth, sibllocs[iloc]+calchalfwidth, hcat(x, y, z))
        sortperms = func.sortprofileentries("z", hcat(x, y, z), profilelines)

        #extract profiles
        tmp = quant_median[profilelines]
        profile = tmp[sortperms]
        dist_along_profile = fill(NaN, size(sortperms, 1))
        tmpdist = z[profilelines]
        dist_along_profile = tmpdist[sortperms]

        sibldepth_idx[iloc] = func.diagnosesibldepth(profile, diagnosefromquantity)#findfirst(x->x>0, func.movingaverage(quant_median[irow, 5:end], 10))
        idcs = sibllocs[iloc]-calchalfwidth .<= xtopo.x .<= sibllocs[iloc]+calchalfwidth
        if sibldepth_idx[iloc] == 0
            if count(idcs )== 0
                tmpmini = abs.(xtopo.x .- minimum(sibllocs[iloc]))
                mintmpmini = minimum(tmpmini)
                sibldepth[iloc] = xtopo.h[findfirst(x->x==mintmpmini, tmpmini)]
            else
                sibldepth[iloc] = 0 #mean(xtopo.h[idcs])
            end
        else
            sibldepth[iloc] = dist_along_profile[sibldepth_idx[iloc]] - mean(xtopo.h[idcs])
        end
        sibldepth2[iloc] = sibldepth[iloc] + mean(xtopo.h[idcs])
    end

    #get rid of NaN values
    sibllocs3 = copy(sibllocs)
    sibldepth3 = copy(sibldepth)
    sibllocs[1.1 .<= sibllocs .<= 1.3] .= NaN
    sibllocs[2.3 .<= sibllocs .<= 2.4] .= NaN
    nans1 = isnan.(sibllocs)
    nans2 = isnan.(sibldepth)
    nans_total = nans1 .|| nans2
    sibllocs4a = sibllocs3[(1.0 .<= sibllocs3 .<= 1.4)]
    sibllocs4b = sibllocs3[(2.2 .<= sibllocs3 .<= 2.5)]
    sibldepth4a = sibldepth3[(1.0 .<= sibllocs3 .<=1.4)]
    sibldepth4b = sibldepth3[(2.2 .<= sibllocs3 .<= 2.5)]

    sibllocs2 = sibllocs
    xbiggerequal0_2 = sibllocs2 .>= 0
    sibllocs = sibllocs3[.!nans_total]
    sibldepth = sibldepth[.!nans_total]

    #fit Brutsaert1982 function
    brutsaert82(x, p) = p[1] .* x.^(p[2])# .+ 4.7 * (z / p[3])) #p[1]=u⋆, p[2]=z₀
    p0 = [0.334, 0.77]#, 0.1] #initial parameters

    xbiggerequal0 = sibllocs .>= 0

    brutsaert82_fit = curve_fit(brutsaert82, sibllocs[xbiggerequal0], sibldepth[xbiggerequal0], p0)#; lower=[lower_fric,lower_z0,lower_L], upper=[upper_fric, upper_z0, upper_L])
    brutsaert82_param = brutsaert82_fit.param
    println("----------------------------------------------------------")
    println("Parameters for SIBL depth fit c*x^b (according to Brutsaert et al. 1982):")
    println(string("c = ", brutsaert82_param[1]))
    println(string("b = ", brutsaert82_param[2]))
    println("----------------------------------------------------------")
    #=
    #plot SIBL depth
    fig = PyPlot.figure(figsize=(8,6))
    ax2 = fig.add_subplot(111)
    #ax2.set_title(string("SIBL depth from ", diagnosefromquantity," & Brutsaert82 fit (h=c*x^b); c=", round(brutsaert82_param[1]*100, digits=2), "e-2, b=", round(brutsaert82_param[2], digits=2)))
    ax2.plot(sibllocs3, sibldepth3, label="LES")
    ax2.plot(sibllocs4a, sibldepth4a, color="white", alpha=0.6)
    ax2.plot(sibllocs4b, sibldepth4b, color="white", alpha=0.6)
    ax2.plot(sibllocs[xbiggerequal0], brutsaert82(sibllocs[xbiggerequal0], brutsaert82_param), "--", label="Brutsaert82 fit")
    ax2.set_xlim(-0, 12)
    ax2.set_ylim(0, 0.50)
    ax2.set_xlabel(L"x~\mathrm{[m]}")
    ax2.set_ylabel(L"h_{SIBL} [m]")
    ax2.grid()
    #ax2.legend()
    =#
end

##
fig = PyPlot.figure(figsize=(15,4.5))
gs = GridSpec.GridSpec(2, 2, width_ratios=(50,1), height_ratios=(2, 1))
ax1 = fig.add_subplot(gs[1, 1])
hm = ax1.tripcolor(x, z, quant_median, cmap=cramericm.vik, vmin=-plotlim, vmax=plotlim, shading="gouraud")
topo = ax1.plot(xtopo.x, xtopo.h, color="black", lw=3)
ax1.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
ax1.set_ylabel(L"h~[\mathrm{m}]")
if calcsibl
    #ax1.plot(sibllocs2[xbiggerequal0_2], sibldepth2[xbiggerequal0_2], lw=3, c="green")
end
ax1.set_ylim(0,2)
#ax1.set_xlim(-0.5, 12)
ax1.tick_params(axis="x", labelbottom=false)
ax1.text(-0.06, 1.05, "a)", size="large", transform=ax1.transAxes)
PyPlot.gca().set_aspect("equal")
ax1c = fig.add_subplot(gs[1, 2])
hm_cb = fig.colorbar(hm, cax=ax1c)
hm_cb.set_label(cblabel)

ax2 = fig.add_subplot(gs[2,1], sharex=ax1)
ax2.plot(sibllocs3, sibldepth3, label="diagnosed SIBL depth")
ax2.plot(sibllocs4a, sibldepth4a, color="white", alpha=0.6)
ax2.plot(sibllocs4b, sibldepth4b, color="white", alpha=0.6)
ax2.plot(sibllocs[xbiggerequal0], brutsaert82(sibllocs[xbiggerequal0], brutsaert82_param), "--", label="power law fit")
ax2.set_xlim(-0.5, 11)
ax2.set_ylim(0, 0.50)
ax2.set_xlabel(L"x~\mathrm{[m]}")
ax2.set_ylabel(L"h_{SIBL} [m]")
ax2.grid()
ax2.text(-0.06, 1.05, "b)", size="large", transform=ax2.transAxes)

(h, l) = ax2.get_legend_handles_labels()
ax_leg = fig.add_subplot(gs[2, 2])
l1 = ax_leg.legend(h, l, loc="best", bbox_to_anchor=(6.5, 0.75))#,borderaxespad=0)
ax_leg.axis("off")
PyPlot.gca().add_artist(l1)

##
#PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/SIBL_depth.jpg")
=#

####################################################################################
####################################################################################
###  Fig. 6:          PROFILES FOR HIGH AND LOW WIND SPEEDS                      ###
####################################################################################
#####################################################################################profiles for low and high wind speeds
##
#=
using NCDatasets, DataFrames, StatsBase, LaTeXStrings, Dates, PyCall
import PyPlot
PyPlot.pygui(true)
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
cramericm = pyimport("cmcrameri.cm")
lines = pyimport("matplotlib.lines")
#outer directory of case
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/third_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#netCDF4-file containing the profile (length 3)
file = joinpath.(joinpath(casedir, "postProcessing/surfaces/yNormal/"), "yNormal_time5th.nc")
#location of .nc-file containing the forcing wind speed
uforcingfile = joinpath(casedir, "forcing.nc")

#define plottype
plottype = "profiles"
allowedplottypes = ["profiles", "points"]

#axis along which profile is taken ("x", "y", or "z")
profilealong = "z"

#value of x where bare to snow transition is
snowtransition = 8.0 #m

#profile location [m] in new fetch-distance coordinates for x (multiple profiles, -1 along profile axis)
profilelocations = [-0.2, 0.5, 2.0, 10.5]

#profile half width
profhalfwidth = 0.03

#quantity for profiles (T, u, v, w, wT, Tvar, tke, or z/L)
quantities_to_plot = ["T", "wT", "u"] #["u", "w", "T", "wT", "tke"]
allowedquantities = ["u", "v", "w", "T", "wT", "tke", "Tvar", "tke", "z/L"]

#reynolds averaging time for simulation output
reyavgtime = Second(50)

#seconds to skip (spinuptime)
skipsec = 10

#flag for filtering of profiles
postfilter = true
postfilterquantity = "u" #u, v, w, or T
#location of point to evaluate for post-filtering (fetch distance coordinates for x) [m]
postfilterloc = [0, 3, 1.5]
#filter condition (not working here!! set manually in code below; separately for sonics!!)
filtercond_list = ["a", "b"]#["0.0 .<= filterdata .<= filterquantiles[20]", "filterquantiles[80] .<= filterdata"]

#flag for loading sonic measurement data as well
loadsonics = true
#which sonics
listsonics = ["kaijo", "t1irg", "t2lcsat"]#, "t2ucsat", "tjk"]
allowedsonics = ["kaijo", "t1irg", "t2irg", "t2lcsat", "t2ucsat", "tjk"]
#period to load from sonics
sonicsstart = DateTime(2021,05,31,12,00,00)
sonicsend = sonicsstart + Minute(10)
drtime = Second(600)

#location of sonic-input files (turbulence data)
tjkinfile =     joinpath(duerrpath, "scripts", "src", "data", "tjktmp.nc")
t1irginfile =   joinpath(duerrpath, "scripts", "src", "data", "t1irgtmp.nc")
t2irginfile =   joinpath(duerrpath, "scripts", "src", "data", "t2irgtmp.nc")
t2lcsatinfile = joinpath(duerrpath, "scripts", "src", "data", "t2lcsattmp.nc")
t2ucsatinfile = joinpath(duerrpath, "scripts", "src", "data", "t2ucsattmp.nc")
kaijoinfile =   joinpath(duerrpath, "scripts", "src", "data", "kaijotmp.nc")

#location of slow data input file
tjkslowfile = joinpath("/home/haugened/Documents/openfoam/duerr_les/", "scripts", "src", "data", "tjk_data.csv")
t2slowfile = joinpath("/home/haugened/Documents/openfoam/duerr_les/", "scripts", "src", "data", "t2slow.nc")
#################################################
#load additional functions if sonics should be included
if loadsonics
    if !@isdefined turb_loaded
        turb_loaded = false
    end
    if !turb_loaded
        include("/home/haugened/Documents/ibl_patch_snow/code/src/turb_data.jl")
        import .turb    
        turb_loaded = true
    end
end
#################################################
#input check
if !(plottype in allowedplottypes)
    println(string(plottype, " is not in ", allowedplottypes, ". Please set variable 'plottype' to allowed value. Exiting."))
    exit()
else
    println("Plot type: ", plottype)
end

if loadsonics
    for i in listsonics
        if !(i in allowedsonics)
            println(string(i, " is not in ", allowedsonics, ". Please set variable 'listsonics' to allowed value. Exiting."))
            exit()
        end
    end
    println("Loading sonic data from ", listsonics)
end

#################################################
#read in data
if !@isdefined forcingdata_read
    forcingdata_read = false
    forcingdata_read_sucessfully = false
end
if !forcingdata_read
    println("Reading forcing data...")
    if isfile(uforcingfile)
        (forcing_heights, forcing_times, forcing_u, forcing_T) = func.readforcingfromnetcdf(uforcingfile)
        forcingdata_read_sucessfully = true
        @info("Forcing data sucessfully read.")
    else
        @warn("Forcing data could not be read.")
    end
    forcingdata_read = true
end

if !@isdefined data_read
    data_read = false
end
if !data_read
    println("Reading cross section data...")
    (x, y, z, times, T, U) = func.readaerialprofilefromnetcdf(file)
    data_read = true
end

#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx
    transformx = false
end
if !transformx
    x = x .- snowtransition
    transformx = true
end

#create bottom line for plot
xtopo = func.createbottom(x, z)

#deal with sonic data
if loadsonics
    if !@isdefined sonics_read
        sonics_read = false
    end
    if !sonics_read
        #sonic heights
        sonic_heights = zeros(Float64, length(listsonics))
        println("Reading sonic data...")
        if "kaijo" in listsonics
            htmp = 0.3
            datatmp = func.readtotalturbasnetcdf(kaijoinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            @warn("Careful!! Rotate Kaijo data due to measurement position. Check pictures of setup!")
            new_w = datatmp.u
            new_u = datatmp.w
            datatmp.u = new_u
            datatmp.w = new_w
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            kaijodata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="kaijo", listsonics)] = htmp
        end
        if "t1irg" in listsonics
            htmp = 1.2
            datatmp = func.readtotalturbasnetcdf(t1irginfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t1irgdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t1irg", listsonics)] = htmp
        end
        if "t2irg" in listsonics
            htmp = 0.9
            datatmp = func.readtotalturbasnetcdf(t2irginfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2irgdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2irg", listsonics)] = htmp
        end
        if "t2lcsat" in listsonics
            htmp = 1.9
            datatmp = func.readtotalturbasnetcdf(t2lcsatinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2lcsatdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2lcsat", listsonics)] = htmp
        end
        if "t2ucsat" in listsonics
            htmp = 2.85
            datatmp = func.readtotalturbasnetcdf(t2ucsatinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2ucsatdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2ucsat", listsonics)] = htmp
        end
        if "tjk" in listsonics
            htmp = 5.0
            datatmp = func.readtotalturbasnetcdf(tjkinfile)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            tjkdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="tjk", listsonics)] = htmp
        end
        sonics_read = true
    end
end

fig = PyPlot.figure(figsize=(10,10))
gs = gridspec.GridSpec(7, 2, width_ratios=[1, 1], height_ratios=[5, 1.5, 5, 1.5, 5, 1, 2])

subfig_letters = ["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)", "i)", "j)", "k)", "l)", "m)"]

for q in 1:size(quantities_to_plot, 1)
    quantity_to_plot = quantities_to_plot[q]
    #take quantity to plot
    if quantity_to_plot == "u"
        vals = U[:,:,1]
        x_label_text = L"u~\mathrm{[m~s^{-1}]}"
        quantity_symbol = :U0
    elseif quantity_to_plot == "v"
        vals = U[:,:,2]
        x_label_text = L"v~\mathrm{[m~s^{-1}]}"
        quantity_symbol = :U1
    elseif quantity_to_plot == "w"
        vals = U[:,:,3]
        x_label_text = L"w~\mathrm{[m~s^{-1}]}"
        quantity_symbol = :U2
    elseif quantity_to_plot == "T"
        vals = T .- 273.15
        x_label_text = L"T~\mathrm{[^\circ C]}"
        #x_label_text = L"T~\mathrm{[K]}"
        quantity_symbol = :T
    elseif quantity_to_plot == "wT"
        vals1 = T
        vals2 = U[:,:,3] #w
        x_label_text = L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}"
        quantity_symbol = :wT
    elseif quantity_to_plot == "Tvar"
        vals = T
        x_label_text = L"\overline{\left( T' \right)^2}~\mathrm{[K^2]}"
        quantity_symbol = :Tvar
    elseif quantity_to_plot == "tke"
        vals1 = U[:,:,1] #u
        vals2 = U[:,:,2] #v
        vals3 = U[:,:,3] #w
        x_label_text = L"e~\mathrm{[m^2~s^{-2}]}"
        quantity_symbol = :tke
    elseif quantity_to_plot == "z/L"
        valsu = U[:,:,1] #u
        valsv = U[:,:,2] #v
        valsw = U[:,:,3] #w
        valsT = T
        x_label_text = L"\frac{z}{L}~\mathrm{[1]}"
        quantity_symbol = :zoverL
    end

    skiptoidx = findfirst(x->x==minimum(abs.(times .- skipsec)), abs.(times .- skipsec))
    println(string("Skipping ", skipsec, " seconds at the beginning."))

    if loadsonics
        if quantity_to_plot == "u"
            sonic_symbol = :u
        elseif quantity_to_plot == "v"
            sonic_symbol = :v
        elseif quantity_to_plot == "w"
            sonic_symbol = :w
        elseif quantity_to_plot == "T"
            sonic_symbol = :T
        elseif quantity_to_plot == "wT"
            sonic_symbol = :wT
        elseif quantity_to_plot == "Tvar"
            sonic_symbol = :TT
        elseif quantity_to_plot == "tke"
            sonic_symbol = :tke
        elseif quantity_to_plot == "z/L"
            sonic_symbol = :zoverL
        end
    end

    tjkmeteodata = func.csvtodataframe(tjkslowfile)
    t2vent = turb.readturbasnetcdf(t2slowfile, sonicsstart, sonicsend + Hour(2)) #t2 is in LT, not in UTC+1
    t2vent.time .-= Hour(1)
    height_t2vent = 1.00 #[m]
    tjkmeteodata_excerpt = tjkmeteodata[sonicsstart .<= tjkmeteodata.time .<= sonicsend, :]
    t2vent_excerpt = t2vent[sonicsstart .<= t2vent.time .<= sonicsend, :]

    meanTair_tjk = mean(tjkmeteodata_excerpt.tair_HygroVUE10)
    meanTair_t2vent = mean(t2vent_excerpt.vent_air_temp)

    scaledTair = t2vent_excerpt.vent_air_temp .* (meanTair_tjk/meanTair_t2vent)

    #choose profile lines
    profilelines = fill(false, size(x, 1), size(profilelocations, 1))
    sortperms = fill(-1, size(x, 1), size(profilelocations, 1))
    for i in 1:length(profilelocations)
        profilelines[:, i] = func.chooselistval(profilealong, profilelocations[i]-profhalfwidth, profilelocations[i]+profhalfwidth, hcat(x, y, z))
        tmp = func.sortprofileentries(profilealong, hcat(x, y, z), profilelines[:,i])
        sortperms[1:length(tmp), i] = tmp
    end

    firstall0idx = 0
    for i in 1:size(sortperms, 2)
        firstall0idx = maximum([firstall0idx, findlast(x->x!=-1, sortperms[:, i])+1])
    end
    sortperms = sortperms[1:firstall0idx-1, :]

    #extract profiles

    profilesu = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
    profilesv = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
    profilesw = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
    profilesT = fill(NaN, size(sortperms,1), length(profilelocations), size(T, 2))
    for i in 1:length(profilelocations)
        tmpu = U[profilelines[:, i],:,1]
        profilesu[1:size(tmpu,1),i,:] = tmpu[sortperms[1:size(tmpu, 1), i], :]
        tmpv = U[profilelines[:, i],:,2]
        profilesv[1:size(tmpu,1),i,:] = tmpv[sortperms[1:size(tmpv, 1), i], :]
        tmpw = U[profilelines[:, i],:, 3]
        profilesw[1:size(tmpu,1),i,:] = tmpw[sortperms[1:size(tmpw, 1), i], :]
        tmpT = T[profilelines[:, i],:]
        profilesT[1:size(tmpu,1),i,:] = tmpT[sortperms[1:size(tmpT, 1), i], :]
    end
    dist_along_profile = fill(NaN, size(sortperms, 1), length(profilelocations))
    if profilealong == "x"
        for i in 1:length(profilelocations)
            tmpdist = x[profilelines[:, i]]
            dist_along_profile[1:length(tmpdist), i] = tmpdist[sortperms[1:size(tmpdist,1), i]]
        end
    elseif profilealong == "z"
        for i in 1:length(profilelocations)
            tmpdist = z[profilelines[:, i]]
            dist_along_profile[1:length(tmpdist), i] = tmpdist[sortperms[1:size(tmpdist,1), i]]
        end
    end

    for i in 1:length(profilelocations)
        (profilesu[:,i,:], ) = func.reducedoubleentries(profilesu[:, i, :], dist_along_profile[:, i], profhalfwidth)
        (profilesv[:,i,:], ) = func.reducedoubleentries(profilesv[:, i, :], dist_along_profile[:, i], profhalfwidth)
        (profilesw[:,i,:], ) = func.reducedoubleentries(profilesw[:, i, :], dist_along_profile[:, i], profhalfwidth)
        (profilesT[:,i,:], dist_along_profile[:, i]) = func.reducedoubleentries(profilesT[:, i, :], dist_along_profile[:, i], profhalfwidth)    
        firstall0idx = 0
    end

    firstall0idx = 0
    for i in 1:size(profilelocations, 2)
        firstall0idx = maximum([firstall0idx, findfirst(x->isnan(x), dist_along_profile[:, i])])
    end
    profilesu = profilesu[1:firstall0idx, :, :]
    profilesv = profilesv[1:firstall0idx, :, :]
    profilesw = profilesw[1:firstall0idx, :, :]
    profilesT = profilesT[1:firstall0idx, :, :]
    dist_along_profile = dist_along_profile[1:firstall0idx, :]

    #get terrain height at profile locations
    xidcs = zeros(Int64, size(profilelocations, 1))
    meanheights = zeros(Float64, size(profilelocations, 1))
    for l in 1:length(xidcs)
        idcstmp = findall(x->profilelocations[l]-profhalfwidth .<= x .<= profilelocations[l]+profhalfwidth, xtopo.x)
        mintmp = minimum(abs.(xtopo.x .- profilelocations[l]))
        xidcs[l] = findfirst(x->x==mintmp, abs.(xtopo.x .- profilelocations[l]))
        meanheights[l] = mean(xtopo.h[idcstmp])
        dist_along_profile[:,l] .-= meanheights[l]
    end

    if quantity_to_plot == "z/L"
        #calculate fluxes
        profiles = zeros(Float64, size(profilesw, 1), size(profilesw, 2), size(profilesw, 3))

        timesteps = func.gettimesteps(times)
        meantimesteps = round(mean(timesteps), digits=3)
        movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
        
        for icol in 1:size(profilesw, 2)
            for jrow in 1:size(profilesw, 1)
                movavgprofu = func.movingaverage(profilesu[jrow, icol, :], movavgidcs)
                movavgprofv = func.movingaverage(profilesv[jrow, icol, :], movavgidcs)
                movavgprofw = func.movingaverage(profilesw[jrow, icol, :], movavgidcs)
                movavgprofT = func.movingaverage(profilesT[jrow, icol, :], movavgidcs)
                devu = profilesu[jrow, icol, :] .- movavgprofu
                devv = profilesv[jrow, icol, :] .- movavgprofv
                devw = profilesw[jrow, icol, :] .- movavgprofw
                devT = profilesT[jrow, icol, :] .- movavgprofT
                ustar = ((devu .* devw) .^2 .+ (devv .* devw) .^2) .^(1/4)
                wT_tmp = devw .* devT
                L_tmp = - (ustar .^ 3 .* movavgprofT) ./ (0.4*9.81 .* wT_tmp)
                if profilealong == "x" .|| profilealong == "y"
                    profiles[jrow, icol, :] = prof_locs_final[icol] ./ L_tmp
                else #profilealong == "z"
                    profiles[jrow, icol, :] = dist_along_profile[jrow] ./ L_tmp
                end
            end
        end

    elseif quantity_to_plot == "Tvar"
        #calculate fluxes
        profiles = zeros(Float64, size(profilesT, 1), size(profilesT, 2), size(profilesT, 3))

        timesteps = func.gettimesteps(times)
        meantimesteps = round(mean(timesteps), digits=3)
        movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
        
        for icol in 1:size(profilesT, 2)
            for jrow in 1:size(profilesT, 1)
                movavgprof1 = func.movingaverage(profilesT[jrow, icol, :], movavgidcs)
                dev1 = profilesT[jrow, icol, :] .- movavgprof1
                profiles[jrow, icol, :] = dev1 .* dev1
            end
        end
    elseif quantity_to_plot == "tke"
        #calculate fluxes
        profiles = zeros(Float64, size(profilesu, 1), size(profilesu, 2), size(profilesu, 3))

        timesteps = func.gettimesteps(times)
        meantimesteps = round(mean(timesteps), digits=3)
        movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
        
        for icol in 1:size(profilesu, 2)
            for jrow in 1:size(profilesu, 1)
                movavgprof1 = func.movingaverage(profilesu[jrow, icol, :], movavgidcs)
                movavgprof2 = func.movingaverage(profilesv[jrow, icol, :], movavgidcs)
                movavgprof3 = func.movingaverage(profilesw[jrow, icol, :], movavgidcs)
                dev1 = profilesu[jrow, icol, :] .- movavgprof1
                dev2 = profilesv[jrow, icol, :] .- movavgprof2
                dev3 = profilesw[jrow, icol, :] .- movavgprof3
                profiles[jrow, icol, :] = 0.5 .*(dev1 .^2 .+ dev2 .^2 .+ dev3 .^2)
            end
        end
    elseif quantity_to_plot == "wT" #fluxes
        #calculate fluxes
        profiles = zeros(Float64, size(profilesw, 1), size(profilesw, 2), size(profilesw, 3))

        timesteps = func.gettimesteps(times)
        meantimesteps = round(mean(timesteps), digits=3)
        movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
        
        for icol in 1:size(profilesw, 2)
            for jrow in 1:size(profilesw, 1)
                movavgprof1 = func.movingaverage(profilesw[jrow, icol, :], movavgidcs)
                movavgprof2 = func.movingaverage(profilesT[jrow, icol, :], movavgidcs)
                dev1 = profilesw[jrow, icol, :] .- movavgprof1
                dev2 = profilesT[jrow, icol, :] .- movavgprof2
                profiles[jrow, icol, :] = dev1 .* dev2
            end
        end
    elseif quantity_to_plot == "u"
        profiles = profilesu
    elseif quantity_to_plot == "v"
        profiles = profilesv
    elseif quantity_to_plot == "w"
        profiles = profilesw
    elseif quantity_to_plot == "T"
        profiles = profilesT .- 273.15
    end
    
    profiles_before_filter = copy(profiles)

    for f in 1:size(filtercond_list, 1)
        filtercond = filtercond_list[f]
        if postfilter
            if q==1
                println("-------------------------------------")
                println("Filtering of profiles activated.")
                #println("No filtering of sonic measurement data!")
                println("Filter condition: ", filtercond)

                #get data at filter location
                global post_take = func.chooseval3D(postfilterloc, hcat(x, y, z))

                println("location for filtering:")
                println(string("x = ", x[post_take], "m"))
                println(string("y = ", y[post_take], "m"))
                println(string("z = ", z[post_take], "m"))
                println("-------------------------------------")
            end

            if postfilterquantity == "u"
                filterdata = U[post_take,:,1]
            elseif postfilterquantity == "v"
                filterdata = U[post_take,:,2]
            elseif postfilterquantity == "w"
                filterdata = U[post_take,:,3]
            elseif postfilterquantity == "T"
                filterdata = T[post_take,:]
            else
                @error("Wrong filter quantity! Please check.")
            end

            #calculate 1%-quantiles
            quantiles = collect(0.01:0.01:1)
            filterquantiles = quantile(filter(!isnan, filterdata), quantiles)

            #post filter condition
            if f == 1
                #postfiltercond = filterquantiles[20] .<= filterdata .<= filterquantiles[40]
                postfiltercond = 0.0 .<= filterdata .<= filterquantiles[20]
            elseif f==2
                #postfiltercond = filterquantiles[60] .<= filterdata .<= filterquantiles[80]
                postfiltercond = filterquantiles[80] .<= filterdata
            end
            #postfiltercond = eval(Meta.parse(filtercond))
        else #no postfiltering
            postfiltercond = fill(true, size(profiles, 3))
            println("No post-filtering applied.")
        end

        profiles = copy(profiles_before_filter)
        profiles[:,:, .!(postfiltercond)] .= NaN
        #println(count(postfiltercond)/size(postfiltercond, 1), " times remaining after filtering")

        medianprof = fill(NaN, size(profiles, 1), size(profiles, 2), 3)

        for j in 1:size(medianprof, 2)
            for i in 1:size(medianprof, 1)
                try
                    (medianprof[i,j,1], medianprof[i,j,2], medianprof[i,j,3]) = quantile(filter(!isnan, profiles[i, j, skiptoidx:end]), [0.25, 0.5, 0.75])
                catch e
                end
            end
        end

        pltrow = (q-1)*2+1
        pltcol = f

        ax1 = fig.add_subplot(gs[pltrow, pltcol])
        #=if postfilter
            ax1.set_title(filtercond)
        else
            ax1.set_title("No filtering")
        end=#
        for k in 1:size(profilelocations, 1)
            ax1.plot(medianprof[:,k,2], dist_along_profile[:, k], label=string("x=", round(profilelocations[k], digits=1), "m"))#-profhalfwidth, " - ", profilelocations[k]+profhalfwidth, "m"))
            ax1.fill_betweenx(dist_along_profile[:, k], medianprof[:,k,1], medianprof[:,k,3], alpha=0.3)
        end
        #ax1.plot(medianprof[:,2,2], dist_along_profile[:, 2], label=string("x=", profilelocations[2]-profhalfwidth, " - ", profilelocations[2]+profhalfwidth, "m"))
        #ax1.fill_betweenx(dist_along_profile[:, 2], medianprof[:,2,1], medianprof[:,2,3], alpha=0.3)
        #ax1.plot(medianprof[:,3,2], dist_along_profile[:, 3], label=string("x=", profilelocations[3]-profhalfwidth, " - ", profilelocations[3]+profhalfwidth, "m"))
        #ax1.fill_betweenx(dist_along_profile[:, 3], medianprof[:,3,1], medianprof[:,3,3], alpha=0.3)
        #calculate mean and IQR for sonics
        if loadsonics .&& quantity_to_plot in ["u", "w", "tke"]
            #post filter condition
            if "kaijo" in listsonics
                sonicdata_filt = copy(kaijodata[:, postfilterquantity])
                sonicdata_col = copy(kaijodata[:, sonic_symbol])
                sonicfilterquantiles = quantile(filter(!isnan, sonicdata_filt), collect(0.01:0.01:1))
                if f == 1
                    sonicpostfiltercond = sonicdata_filt .<= sonicfilterquantiles[20]
                elseif f==2
                    sonicpostfiltercond = sonicfilterquantiles[80] .<= sonicdata_filt
                end
                sonicdata_col[.!(sonicpostfiltercond)] .= NaN
                sonicmed = quantile(filter(!isnan, sonicdata_col), [0.25,0.5,0.75])
                ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="kaijo", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="short-path ultrasonic")
            end
            if "t1irg" in listsonics
                sonicdata_filt = copy(t1irgdata[:, postfilterquantity])
                sonicdata_col = copy(t1irgdata[:, sonic_symbol])
                sonicfilterquantiles = quantile(filter(!isnan, sonicdata_filt), collect(0.01:0.01:1))
                if f == 1
                    sonicpostfiltercond = sonicdata_filt .<= sonicfilterquantiles[20]
                elseif f==2
                    sonicpostfiltercond = sonicfilterquantiles[80] .<= sonicdata_filt
                end
                sonicdata_col[.!(sonicpostfiltercond)] .= NaN
                sonicmed = quantile(filter(!isnan, sonicdata_col), [0.25,0.5,0.75])
                ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t1irg", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="ECT1, h=1.2m")
            end
            if "t2irg" in listsonics
                sonicdata_filt = copy(t2irgdata[:, postfilterquantity])
                sonicdata_col = copy(t2irgdata[:, sonic_symbol])
                sonicfilterquantiles = quantile(filter(!isnan, sonicdata_filt), collect(0.01:0.01:1))
                if f == 1
                    sonicpostfiltercond = sonicdata_filt .<= sonicfilterquantiles[20]
                elseif f==2
                    sonicpostfiltercond = sonicfilterquantiles[80] .<= sonicdata_filt
                end
                sonicdata_col[.!(sonicpostfiltercond)] .= NaN

                sonicmed = quantile(filter(!isnan, sonicdata_col), [0.25,0.5,0.75])
                ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t2irg", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="ECT2, h=0.9m")
            end
            if "t2lcsat" in listsonics
                sonicdata_filt = copy(t2lcsatdata[:, postfilterquantity])
                sonicdata_col = copy(t2lcsatdata[:, sonic_symbol])
                sonicfilterquantiles = quantile(filter(!isnan, sonicdata_filt), collect(0.01:0.01:1))
                if f == 1
                    sonicpostfiltercond = sonicdata_filt .<= sonicfilterquantiles[20]
                elseif f==2
                    sonicpostfiltercond = sonicfilterquantiles[80] .<= sonicdata_filt
                end
                sonicdata_col[.!(sonicpostfiltercond)] .= NaN
                sonicmed = quantile(filter(!isnan, sonicdata_col), [0.25,0.5,0.75])
                ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t2lcsat", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="ECT2, h=1.9m")
            end
            if "t2ucsat" in listsonics
                sonicdata_filt = copy(t2ucsatdata[:, postfilterquantity])
                sonicdata_col = copy(t2ucsatdata[:, sonic_symbol])
                sonicfilterquantiles = quantile(filter(!isnan, sonicdata_filt), collect(0.01:0.01:1))
                if f == 1
                    sonicpostfiltercond = sonicdata_filt .<= sonicfilterquantiles[20]
                elseif f==2
                    sonicpostfiltercond = sonicfilterquantiles[80] .<= sonicdata_filt
                end
                sonicdata_col[.!(sonicpostfiltercond)] .= NaN
                sonicmed = quantile(filter(!isnan, sonicdata_col), [0.25,0.5,0.75])
                ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t2ucsat", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="ECT2, h=2.9m")
            end
            if "tjk" in listsonics
                sonicdata_filt = copy(tjkdata[:, postfilterquantity])
                sonicdata_col = copy(tjkdata[:, sonic_symbol])
                sonicfilterquantiles = quantile(filter(!isnan, sonicdata_filt), collect(0.01:0.01:1))
                if f == 1
                    sonicpostfiltercond = sonicdata_filt .<= sonicfilterquantiles[20]
                elseif f==2
                    sonicpostfiltercond = sonicfilterquantiles[80] .<= sonicdata_filt
                end
                sonicdata_col[.!(sonicpostfiltercond)] .= NaN
                sonicmed = quantile(filter(!isnan, sonicdata_col), [0.25,0.5,0.75])
                ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="tjk", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="TJK (0.25-0.75)")
            end
        end
        #if quantity_symbol == :T
            #t2vent_excerpt[:, :vent_air_temp]
            #ventmed = quantile(filter(!isnan, scaledTair), [0.25,0.5,0.75])
            #ax1.errorbar(ventmed[2], 1.0, xerr=[[abs(ventmed[1]-ventmed[2])], [abs(ventmed[3]-ventmed[2])]], fmt="o", label=L"T_{vent}~\mathrm{(scaled)}")
        #end
        ax1.set_xlabel(x_label_text)
        ax1.set_ylim(0,2.0)

        if q == 1
            if f==1
                texttmp = L"u \leq u_{0.2}"
                ax1.text(0.4, 1.2, texttmp, size="large", transform=ax1.transAxes)
            elseif f==2
                texttmp = L"u \geq u_{0.8}"
                ax1.text(0.37, 1.2, texttmp, size="large", transform=ax1.transAxes)
            end
        end

        if quantity_to_plot == "u"
            ax1.set_xlim(-0.2, 4.3)
        elseif quantity_to_plot == "w"
            ax1.set_xlim(-0.11, 0.11)
        elseif quantity_to_plot == "T"
            ax1.set_xlim(0, 15)
        elseif quantity_to_plot =="wT"
            ax1.set_xlim(-0.05, 0.05)
        elseif quantity_to_plot == "tke"
            ax1.set_xlim(0,0.7)
        end
        ax1.grid()
        if f == 2
            ax1.tick_params(axis="y", labelleft=false)
            ax1.text(-0.1, 1.03, subfig_letters[(q-1)*2+f], size="large", transform=ax1.transAxes)
        else
            ax1.set_ylabel(L"h~\mathrm{[m]}")
            ax1.text(-0.15, 1.03, subfig_letters[(q-1)*2+f], size="large", transform=ax1.transAxes)
        end
        if quantity_to_plot == "u" && f==1 #q==size(quantities_to_plot, 1)
            (h, l) = ax1.get_legend_handles_labels()
            println(h)
            println(l)
            ax_leg = fig.add_subplot(py"$(gs)[6, :]")
            l1 = ax_leg.legend(h[1:size(profilelocations, 1)], l[1:size(profilelocations, 1)], loc="best", bbox_to_anchor=(0.5, 0.6))#,borderaxespad=0)
            l2 = ax_leg.legend(h[size(profilelocations, 1)+1:end], l[size(profilelocations, 1)+1:end], loc="best", bbox_to_anchor=(0.77, 0.6))#,borderaxespad=0)
            ax_leg.axis("off")
            PyPlot.gca().add_artist(l1)
            ax_leg.text(0.4, 0.62, "LES", size="medium", transform=ax_leg.transAxes)
            ax_leg.text(0.56, 0.62, "measurements", size="medium", transform=ax_leg.transAxes)
        end
    end #postfilter loop
end #quantity_to_plot loop
#PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/quant_profiles_3rd_new.jpg")
=#
##
####################################################################################
####################################################################################
###  Fig. 7:   BUOYANCY FLUXES AND SIBL GROWTH FOR HIGH AND LOW WIND SPEEDS      ###
####################################################################################
#####################################################################################SIBL depth growth for low and high wind speeds
#=
using NCDatasets, DataFrames, StatsBase, LaTeXStrings, PyCall, Dates, LsqFit, CSV
import PyPlot
PyPlot.pygui(true)
GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")
animation = pyimport("matplotlib.animation")
numpy = pyimport("numpy")
cramericm = pyimport("cmcrameri.cm")

#outer directory of case
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/third_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#folder containing the time steps
readfromfolder = joinpath(casedir, "postProcessing/surfaces/yNormal")

#netCDF4-file containing the profile (length 3); use same time steps!!
file = joinpath.(readfromfolder, "yNormal_time5th.nc")

#value of x where bare to snow transition is
snowtransition = 8.0 #m

reyavgtime = Second(50)

#seconds to skip (spinuptime)
skipsec = 10

#flag for computing and plotting SIBL depth
calcsibl = true
#distance between two neighboring profiles calculating SIBL
calcevery = 0.1 #m
calchalfwidth = 0.03 #m

#quantity used to diagnose SIBL depth
diagnosefromquantity = "wT"
allowedquantitiestodiagnosefrom = ["wT", "z/L"]

#flag for filtering of profiles (see below in code to specify filter condition)
postfilter = true
postfilterquantity = "u" #u, v, w, or T
#location of point to evaluate for post-filtering
postfilterloc = [0, 3, 1.5]
#filter condition
filtercond_list = ["filterdata .<= filterquantiles[20]", "filterdata .>= filterquantiles[80]"] #not read! (just number of entries) change in code further down!!

#Interval [ms] for plotting.
#intervalmsec = 30

#calculated for video
#framespersec = 1000/intervalmsec

#file name of video
#anistring = string(joinpath(casedir, "video.mp4"))
#################################################
#input check
if !(diagnosefromquantity in allowedquantitiestodiagnosefrom)
    println(string(diagnosefromquantity), " can not be used to diagnose SIBL depth. Possible values are ", allowedquantities, ". Exiting.")
    exit()
else
    println("Diagnosing SIBL depth from ", diagnosefromquantity)
end

#################################################
#read in data
if !@isdefined data_read
    data_read = false
end
if !data_read
    println("Reading cross section data...")
    (x, y, z, times, T, U) = func.readaerialprofilefromnetcdf(file)
    data_read = true
end

#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx
    transformx = false
end
if !transformx
    x = x .- snowtransition
    transformx = true
end

#create bottom line for plot
xtopo = func.createbottom(x, z)

#################################################
if diagnosefromquantity == "wT"
    plotlim = 10
    cblabel = L"\overline{w'T'}~\mathrm{[mK~m~s^{-1}]}"
    if postfilter
        figtitle = string("Median w'T' skipping first ", skipsec, "s; ", filtercond_list)
    else
        figtitle = string("Median w'T' skipping first ", skipsec, "s; No filtering")
    end

    #get stuff for evaluating wT
    vals_w = U[:,:,3] .* 1000
    vals_T = T

    #calculate fluxes
    quant_crosssection = similar(vals_T)
    fill!(quant_crosssection, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))

    for jrow in 1:size(vals_w, 1)
            movavg_w = func.movingaverage(vals_w[jrow, :], movavgidcs)
            movavg_T = func.movingaverage(vals_T[jrow, :], movavgidcs)
            dev1 = vals_w[jrow, :] .- movavg_w
            dev2 = vals_T[jrow, :] .- movavg_T
            quant_crosssection[jrow, :] = dev1 .* dev2
    end
elseif diagnosefromquantity == "z/L"
    plotlim = 0.04
    cblabel = L"\frac{z}{L}~\mathrm{[1]}"
    if postfilter
        figtitle = string("Median z/L skipping first ", skipsec, "s; ", filtercond)
    else
        figtitle = string("Median z/L skipping first ", skipsec, "s; No filtering")
    end


    valsu = U[:,:,1]
    valsv = U[:,:,2]
    valsw = U[:,:,3]
    valsT = T

    #calculate fluxes
    quant_crosssection = similar(valsT)
    fill!(quant_crosssection, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))

    for jrow in 1:size(valsw, 1)
            movavg_u = func.movingaverage(valsu[jrow, :], movavgidcs)
            movavg_v = func.movingaverage(valsv[jrow, :], movavgidcs)
            movavg_w = func.movingaverage(valsw[jrow, :], movavgidcs)
            movavg_T = func.movingaverage(valsT[jrow, :], movavgidcs)
            dev_u = valsu[jrow, :] .- movavg_u
            dev_v = valsv[jrow, :] .- movavg_v
            dev_w = valsw[jrow, :] .- movavg_w
            dev_T = valsT[jrow, :] .- movavg_T
            ustar = ((dev_u .* dev_w) .^2 .+ (dev_v .* dev_w) .^2) .^(1/4)
            wT_tmp = dev_w .* dev_T
            L_tmp = - (ustar .^ 3 .* movavg_T) ./ (0.4*9.81 .* wT_tmp)
            quant_crosssection[jrow, :] = z[jrow] ./ L_tmp
    end
end

#timesround = round.(times, digits=3)
#animidcs = func.getanimidcs(timesround, intervalmsec)

quant_crosssection_unchanged = copy(quant_crosssection)

#get rid of NaN column
x_anim = x

quant_median = fill(NaN, size(quant_crosssection, 1))

##
fig = PyPlot.figure(figsize=(15,10))
gs = GridSpec.GridSpec(5, 2, width_ratios=(50,1), height_ratios=(2, 1, 0.5, 2, 1))
for fi in 1:size(filtercond_list, 1)
    filtercond = filtercond_list[fi]
    if postfilter
        println("-------------------------------------")
        println("Filtering of profiles activated.")
        println("No filtering of sonic measurement data!")
        println("Filter condition: ", filtercond)

        #get data at filter location
        post_take = func.chooseval3D(postfilterloc, hcat(x, y, z))

        println("location for filtering:")
        println(string("x = ", x[post_take], "m"))
        println(string("y = ", y[post_take], "m"))
        println(string("z = ", z[post_take], "m"))
        println("-------------------------------------")

        filterdata = fill(NaN, 0)
        if postfilterquantity == "u"
            filterdata = U[post_take,:,1]
        elseif postfilterquantity == "v"
            filterdata = U[post_take,:,2]
        elseif postfilterquantity == "w"
            filterdata = U[post_take,:,3]
        elseif postfilterquantity == "T"
            filterdata = T[post_take,:]
        else
            @error("Wrong filter quantity! Please check.")
        end
        
        #calculate 1%-quantiles
        quantiles = collect(0.01:0.01:1)
        filterquantiles = quantile(filter(!isnan, filterdata), quantiles)

        println("s")

        #post filter condition
        if fi == 1
            postfiltercond = filterdata .<= filterquantiles[20]
        elseif fi == 2
            postfiltercond = filterdata .>= filterquantiles[80]
        end

    else #no postfiltering
        postfiltercond = fill(true, size(quant_crosssection, 3))
        println("No post-filtering applied.")
    end

    quant_crosssection = copy(quant_crosssection_unchanged)
    quant_crosssection[:, .!(postfiltercond)] .= NaN

    skiptoidx = findfirst(x->x==minimum(abs.(times .- skipsec)), abs.(times .- skipsec))
    println(string("Skipping ", skipsec, " seconds at the beginning for calculating median."))

    for jrow in 1:size(quant_median, 1)
        quant_median[jrow] = median(filter(!isnan, quant_crosssection[jrow, skiptoidx:end]))
    end

    if calcsibl
        #diagnose SIBL depth according to w'T'
        println("Diagnosing SIBL depth from ", diagnosefromquantity)

        sibllocs = collect(minimum(x):calcevery:maximum(x))
        sibldepth_idx = zeros(Int, size(sibllocs, 1))
        sibldepth = zeros(Float64, size(sibllocs, 1))
        sibldepth2 = zeros(Float64, size(sibllocs, 1))

        for iloc in 1:size(sibllocs, 1)
            sortperms = fill(-1, size(x, 1))
            profilelines = func.chooselistval("z", sibllocs[iloc]-calchalfwidth, sibllocs[iloc]+calchalfwidth, hcat(x, y, z))
            sortperms = func.sortprofileentries("z", hcat(x, y, z), profilelines)

            #extract profiles
            tmp = quant_median[profilelines]
            profile = tmp[sortperms]
            dist_along_profile = fill(NaN, size(sortperms, 1))
            tmpdist = z[profilelines]
            dist_along_profile = tmpdist[sortperms]

            sibldepth_idx[iloc] = func.diagnosesibldepth(profile, diagnosefromquantity)#findfirst(x->x>0, func.movingaverage(quant_median[irow, 5:end], 10))
            idcs = sibllocs[iloc]-calchalfwidth .<= xtopo.x .<= sibllocs[iloc]+calchalfwidth
            if sibldepth_idx[iloc] == 0
                if count(idcs )== 0
                    tmpmini = abs.(xtopo.x .- minimum(sibllocs[iloc]))
                    mintmpmini = minimum(tmpmini)
                    sibldepth[iloc] = xtopo.h[findfirst(x->x==mintmpmini, tmpmini)]
                else
                    sibldepth[iloc] = 0 #mean(xtopo.h[idcs])
                end
            else
                sibldepth[iloc] = dist_along_profile[sibldepth_idx[iloc]] - mean(xtopo.h[idcs])
            end
            sibldepth2[iloc] = sibldepth[iloc] + mean(xtopo.h[idcs])
        end

        #get rid of NaN values
        sibllocs3 = copy(sibllocs)
        sibldepth3 = copy(sibldepth)
        if fi == 1
            sibllocs[1.1 .<= sibllocs .<= 1.2] .= NaN
            sibllocs4a = sibllocs3[(1.0 .<= sibllocs3 .<= 1.3)]
            sibldepth4a = sibldepth3[(1.0 .<= sibllocs3 .<=1.3)]
            sibllocs[2.1 .<= sibllocs .<= 2.5] .= NaN
            sibllocs4b = sibllocs3[(2.0 .<= sibllocs3 .<= 2.6)]
            sibldepth4b = sibldepth3[(2.0 .<= sibllocs3 .<= 2.6)]
            sibllocs[3.9 .<= sibllocs .<= 4.2] .= NaN
            sibllocs4c = sibllocs3[(3.8 .<= sibllocs3 .<= 4.3)]
            sibldepth4c = sibldepth3[(3.8 .<= sibllocs3 .<= 4.3)]
            sibllocs[7.0 .<= sibllocs .<= 7.1] .= NaN
            sibllocs4d = sibllocs3[(6.9 .<= sibllocs3 .<= 7.2)]
            sibldepth4d = sibldepth3[(6.9 .<= sibllocs3 .<= 7.2)]
            sibllocs[7.45 .<= sibllocs .<= 7.55] .= NaN
            sibllocs4e = sibllocs3[(7.4 .<= sibllocs3 .<= 7.6)]
            sibldepth4e = sibldepth3[(7.4 .<= sibllocs3 .<= 7.6)]
            sibllocs[8.0 .<= sibllocs .<= 8.3] .= NaN
            sibllocs4f = sibllocs3[(7.9 .<= sibllocs3 .<= 8.4)]
            sibldepth4f = sibldepth3[(7.9 .<= sibllocs3 .<= 8.4)]
            sibllocs[8.65 .<= sibllocs .<= 8.75] .= NaN
            sibllocs4g = sibllocs3[(8.6 .<= sibllocs3 .<= 8.8)]
            sibldepth4g = sibldepth3[(8.6 .<= sibllocs3 .<= 8.8)]
        elseif fi == 2
            sibllocs[1.0 .<= sibllocs .<= 1.3] .= NaN
            sibllocs4a = sibllocs3[(0.9 .<= sibllocs3 .<= 1.4)]
            sibldepth4a = sibldepth3[(0.9 .<= sibllocs3 .<= 1.4)]
            sibllocs[2.0 .<= sibllocs .<= 2.5] .= NaN
            sibllocs4b = sibllocs3[(1.9 .<= sibllocs3 .<= 2.6)]
            sibldepth4b = sibldepth3[(1.9 .<= sibllocs3 .<= 2.6)]
        end
        nans1 = isnan.(sibllocs)
        nans2 = isnan.(sibldepth)
        nans_total = nans1 .|| nans2

        sibllocs2 = sibllocs
        xbiggerequal0_2 = sibllocs2 .>= 0
        sibllocs = sibllocs3[.!nans_total]
        sibldepth = sibldepth[.!nans_total]

        #fit Brutsaert1982 function
        brutsaert82(x, p) = p[1] .* x.^(p[2])# .+ 4.7 * (z / p[3])) #p[1]=u⋆, p[2]=z₀
        p0 = [0.334, 0.77]#, 0.1] #initial parameters

        xbiggerequal0 = sibllocs .>= 0

        brutsaert82_fit = curve_fit(brutsaert82, sibllocs[xbiggerequal0], sibldepth[xbiggerequal0], p0)#; lower=[lower_fric,lower_z0,lower_L], upper=[upper_fric, upper_z0, upper_L])
        brutsaert82_param = brutsaert82_fit.param
        println("----------------------------------------------------------")
        println("Parameters for SIBL depth fit c*x^b (according to Brutsaert et al. 1982):")
        println(string("c = ", brutsaert82_param[1]))
        println(string("b = ", brutsaert82_param[2]))
        println("----------------------------------------------------------")


        if fi == 1
            uprow = 1
            subfig_text = ["a)", "b)"]
            figtitle = L"u \leq u_{0.20}"
        elseif fi == 2
            uprow = 4
            subfig_text = ["c)", "d)"]
            figtitle = L"u \geq u_{0.80}"
        end
        ax1 = fig.add_subplot(gs[uprow, 1])
        hm = ax1.tripcolor(x, z, quant_median, cmap=cramericm.vik, vmin=-plotlim, vmax=plotlim, shading="gouraud")
        topo = ax1.plot(xtopo.x, xtopo.h, color="black", lw=3)
        ax1.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
        ax1.set_ylabel(L"h~[\mathrm{m}]")
        if calcsibl
            #ax1.plot(sibllocs2[xbiggerequal0_2], sibldepth2[xbiggerequal0_2], lw=3, c="green")
        end
        ax1.set_ylim(0,2)
        ax1.tick_params(axis="x", labelbottom=false)
        ax1.text(-0.06, 1.05, subfig_text[1], size="large", transform=ax1.transAxes)
        ax1.set_title(figtitle, fontweight="bold", fontsize=16)
        PyPlot.gca().set_aspect("equal")
        ax1c = fig.add_subplot(gs[uprow, 2])
        hm_cb = fig.colorbar(hm, cax=ax1c)
        hm_cb.set_label(cblabel)

        ax2 = fig.add_subplot(gs[uprow+1,1], sharex=ax1)
        ax2.plot(sibllocs3, sibldepth3, label="diagnosed SIBL depth")
        ax2.plot(sibllocs4a, sibldepth4a, color="white", alpha=0.6)
        ax2.plot(sibllocs4b, sibldepth4b, color="white", alpha=0.6)
        if fi == 1
            ax2.plot(sibllocs4c, sibldepth4c, color="white", alpha=0.6)
            ax2.plot(sibllocs4d, sibldepth4d, color="white", alpha=0.6)
            ax2.plot(sibllocs4e, sibldepth4e, color="white", alpha=0.6)
            ax2.plot(sibllocs4f, sibldepth4f, color="white", alpha=0.6)
            ax2.plot(sibllocs4g, sibldepth4g, color="white", alpha=0.6)
        end
        ax2.plot(sibllocs[xbiggerequal0], brutsaert82(sibllocs[xbiggerequal0], brutsaert82_param), "--", label="power law fit")
        ax2.set_xlim(-0.5, 11)
        ax2.set_ylim(0, 0.60)
        ax2.set_xlabel(L"x~\mathrm{[m]}")
        ax2.set_ylabel(L"h_{SIBL} [m]")
        ax2.grid()
        ax2.text(-0.06, 1.05, subfig_text[2], size="large", transform=ax2.transAxes)

        (h, l) = ax2.get_legend_handles_labels()
        ax_leg = fig.add_subplot(gs[uprow+1, 2])
        l1 = ax_leg.legend(h, l, loc="best", bbox_to_anchor=(6.5, 0.75))#,borderaxespad=0)
        ax_leg.axis("off")
        PyPlot.gca().add_artist(l1)
    end
end
PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/SIBL_depth_quant.jpg")
##
=#
####################################################################################
####################################################################################
###  SUP XXX:                         VIDEO                                      ###
####################################################################################
#####################################################################################video of the LES simulation
#=
using NCDatasets, DataFrames, StatsBase, LaTeXStrings, PyCall, Dates, ProgressMeter, CSV
import PyPlot
PyPlot.pygui(true)
GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")
animation = pyimport("matplotlib.animation")
numpy = pyimport("numpy")
cramericm = pyimport("cmcrameri.cm")

#outer directory of case
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/third_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#folder containing the time steps
readfromfolder = joinpath(casedir, "postProcessing/surfaces/yNormal")

#netCDF4-file containing the profile (length 3); use same time steps!!
file = joinpath.(readfromfolder, "yNormal_time5th.nc")

#value of x where bare to snow transition is
snowtransition = 8.0 #m

#quantity for profiles 
quantity_to_plot = "T"
allowedquantities = ["u", "v", "w", "T", "wT", "tke", "Tvar", "tke", "z/L"]

#only for flux-type variables
reyavgtime = Second(50)

#flag for filtering of profiles (only for plotstyle == "median")
postfilter = false
postfilterquantity = "u" #u, v, w, or T
#location of point to evaluate for post-filtering (fetch distance coordinates for x) [m]
postfilterloc = [-1.0, 2, 0.5]
#filter condition
filtercond = "0.0 .<= filterdata .<= filterquantiles[50]"

#seconds to skip (spinuptime); only for plotstyle == "median"
skipsec = 10

#Interval [ms] for plotting (only for plotstyle == "animation").
intervalmsec = 30

#calculated for video (only for plotstyle == "animation")
framespersec = 1000/intervalmsec

#flag if animation should be saved
anisave = false #true

#file name of video (only for plotstyle == "animation")
anistring = string(joinpath("/home/haugened/Documents/les_paper_24/figures/", "video_T_pub.mp4"))

#################################################
#input check
if !(quantity_to_plot in allowedquantities)
    println(string(quantity_to_plot), " is not in ", allowedquantities, ". Please set variable 'quantity_to_plot' to allowed value. Exiting.")
    exit()
else
    println("Plotting cross section of ", quantity_to_plot)
end
#################################################
#read in data
if !@isdefined data_read
    data_read = false
end
if !data_read
    println("Reading cross section data...")
    (x, y, z, times, T, U) = func.readaerialprofilefromnetcdf(file)
    data_read = true
end

#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx
    transformx = false
end
if !transformx
    x = x .- snowtransition
    transformx = true
end

#create bottom line for plot
xtopo = func.createbottom(x, z)

#################################################
#take and/or quantity to plot
if quantity_to_plot == "u"
    vals = U[:,:,1]
    valmin = -2
    valmax = 2
    colmap = "viridis"
    cblabel = L"u~\mathrm{[m~s^{-1}]}"
    figtitle = string("u")
elseif quantity_to_plot == "v"
    vals = U[:,:,2]
    valmin = -1
    valmax = 2
    colmap = "viridis"
    cblabel = L"v~\mathrm{[m~s^{-1}]}"
    figtitle = string("v")
elseif quantity_to_plot == "w"
    vals = U[:,:,3]
    valmin = -1
    valmax = 1
    colmap = "viridis"
    cblabel = L"w~\mathrm{[m~s^{-1}]}"
    figtitle = string("w")
elseif quantity_to_plot == "T"
    vals = T .-273.15
    valmin = 9
    valmax = 13
    colmap = "viridis"
    cblabel = L"T~\mathrm{[^\circ C]}"
    figtitle = string("T")
elseif quantity_to_plot == "wT"
    vals1 = T
    vals2 = U[:,:,3] #w
    valmin = -0.1
    valmax = -valmin
    colmap = cramericm.vik
    cblabel = L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}"
    figtitle = string("w'T'")

    #calculate fluxes
    vals = similar(vals1)
    fill!(vals, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
    
    @showprogress "Calculating fluxes" for jrow in 1:size(vals1, 1)
        movavg1 = func.movingaverage(vals1[jrow, :], movavgidcs)
        movavg2 = func.movingaverage(vals2[jrow, :], movavgidcs)
        dev1 = vals1[jrow, :] .- movavg1
        dev2 = vals2[jrow, :] .- movavg2
        vals[jrow, :] = dev1 .* dev2
    end
elseif quantity_to_plot == "Tvar"
    vals1 = T
    valmin = 0
    valmax = 2.0
    colmap = "viridis"
    cblabel = L"\overline{\left( T' \right)^2}~\mathrm{[K^2]}"
    figtitle = string("Tvar")

    #calculate fluxes
    vals = similar(vals1)
    fill!(vals, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
    
    @showprogress "Calculating fluxes" for jrow in 1:size(vals1, 1)
            movavg1 = func.movingaverage(vals1[jrow, :], movavgidcs)
            dev1 = vals1[jrow, :] .- movavg1
            vals[jrow, :] = dev1 .* dev1
    end
elseif quantity_to_plot == "tke"
    vals1 = U[:,:,1] #u
    vals2 = U[:,:,2] #v
    vals3 = U[:,:,3] #w
    vmin = 0.0
    vmax = 2.0
    colmap = "viridis"
    cblabel = L"e~\mathrm{[m^2~s^{-2}]}"
    figtitle = string("TKE")

    #calculate fluxes
    vals = similar(vals1)
    fill!(vals, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
    
    @showprogress "Calculating fluxes" for jrow in 1:size(vals1, 1)
            movavg1 = func.movingaverage(vals1[jrow, :], movavgidcs)
            movavg2 = func.movingaverage(vals2[jrow, :], movavgidcs)
            movavg3 = func.movingaverage(vals3[jrow, :], movavgidcs)
            dev1 = vals1[jrow, :] .- movavg1
            dev2 = vals2[jrow, :] .- movavg2
            dev3 = vals3[jrow, :] .- movavg3
            vals[jrow, :] = 0.5 .*(dev1 .^2 .+ dev2 .^2 .+ dev3 .^2)
    end
elseif quantity_to_plot == "z/L"
    valsu = U[:,:,1] #u
    valsv = U[:,:,2] #v
    valsw = U[:,:,3] #w
    valsT = T[:,:]
    valmin = -0.1
    valmax = -valmin
    colmap = cramericm.vik
    cblabel = L"\frac{z}{L}~\mathrm{[1]}"
    figtitle = string("z/L")

    #calculate fluxes
    vals = similar(valsu)
    fill!(vals, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
    
    @showprogress "Calculating fluxes" for jrow in 1:size(valsu, 1)
            movavgu = func.movingaverage(valsu[jrow, :], movavgidcs)
            movavgv = func.movingaverage(valsv[jrow, :], movavgidcs)
            movavgw = func.movingaverage(valsw[jrow, :], movavgidcs)
            movavgT = func.movingaverage(valsT[jrow, :], movavgidcs)
            devu = valsu[jrow, :] .- movavgu
            devv = valsv[jrow, :] .- movavgv
            devw = valsw[jrow, :] .- movavgw
            devT = valsT[jrow, :] .- movavgT
            ustar = ((devu .* devw) .^2 .+ (devv .* devw) .^2) .^(1/4)
            wT_tmp = devw .* devT
            L_tmp = - (ustar .^ 3 .* movavgT) ./ (0.4*9.81 .* wT_tmp)
            vals[jrow, :] = z[jrow] ./ L_tmp
    end
end

#timesround = round.(times, digits=3)
animidcs = func.getanimidcs(times, intervalmsec)

"""
"Update-function for the animation"
animupdate(frame::Int64)
"""
function animupdate(frame::Int64)
    #ax1.set_title(string(figtitle, "; LES t = ", round(times[frame], digits=1), " s"))
    newdata = vals[:,frame]
    #newdata = newdata[1:end-1]
    #hm.set_array(numpy.ravel(newdata))
    hm.set_array(newdata)
    if mod(frame,100) == 0
        println(string(frame, "/", animidcs[end]))
    end
    return hm
end

fig = PyPlot.figure(figsize=(23,4))
ax1 = fig.add_subplot(111)
ax1.set_title("Haugeneder et al. 2024")
hm = ax1.tripcolor(x, z, vals[:, 100], cmap=colmap, vmin=valmin, vmax=valmax, shading="gouraud")
topo = ax1.plot(xtopo.x, xtopo.h, color="black", lw=3)
ax1.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
ax1.set_xlim(minimum(x), maximum(x))
ax1.set_ylim(0, maximum(z))
ax1.set_xlabel(L"x~[\mathrm{m}]")
ax1.set_ylabel(L"z~[\mathrm{m}]")
ax1.set_xlim(-0.5, 11)
ax1.set_ylim(0, 2)
hm_cb = fig.colorbar(hm, ax=ax1)
hm_cb.set_label(cblabel)
PyPlot.gca().set_aspect("equal")
ani = animation.FuncAnimation(fig, animupdate, frames=animidcs[1:end], interval=intervalmsec, repeat=false)#, blit=false)
if anisave
    ani.save(anistring, fps=framespersec)
end
=#






















####################################################################################
####################################################################################
###                              NOT USED FIGURES                                ###
####################################################################################
####################################################################################

####################################################################################
####################################################################################
###  Fig. XX:                     MEDIAN PROFILES                                ###
####################################################################################
#####################################################################################median profiles
##
#=
using NCDatasets, DataFrames, StatsBase, LaTeXStrings, Dates, PyCall
import PyPlot
#PyPlot.pygui(true)
gridspec = pyimport("matplotlib.gridspec")
cm = pyimport("matplotlib.cm")
cramericm = pyimport("cmcrameri.cm")
lines = pyimport("matplotlib.lines")
#outer directory of case
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/third_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#netCDF4-file containing the profile (length 3)
file = joinpath.(joinpath(casedir, "postProcessing/surfaces/yNormal/"), "yNormal_time5th.nc")
#location of .nc-file containing the forcing wind speed
uforcingfile = joinpath(casedir, "forcing.nc")

#define plottype
plottype = "profiles"
allowedplottypes = ["profiles", "points"]

#axis along which profile is taken ("x", "y", or "z")
profilealong = "z"

#value of x where bare to snow transition is
snowtransition = 8.0 #m

#profile location [m] in new fetch-distance coordinates for x (multiple profiles, -1 along profile axis)
profilelocations = [-0.2, 0.5, 2.0, 10.5]

#profile half width
profhalfwidth = 0.03

#quantity for profiles (T, u, v, w, wT, Tvar, tke, or z/L)
quantities_to_plot = ["u", "w", "T", "wT", "tke"]
allowedquantities = ["u", "v", "w", "T", "wT", "tke", "Tvar", "tke", "z/L"]

#reynolds averaging time for simulation output
reyavgtime = Second(50)

#seconds to skip (spinuptime)
skipsec = 10

#flag for filtering of profiles
postfilter = false
postfilterquantity = "u" #u, v, w, or T
#location of point to evaluate for post-filtering (fetch distance coordinates for x) [m]
postfilterloc = [-1.0, 0.5, 2.0]
#filter condition
filtercond = "0.0 .<= filterdata .<= filterquantiles[50]"

#locations of the point data (fetch distance coordinates for x) [m]
p1 = [1.0, 5.0, 0.3]
p2 = [1.0, 5, 1.5]
p3 = [1.0, 5.0, 2.0]

#height to choose from forcing data (if ==-1, then use same as first point) [m]
forcingheight = -1

#flag for loading sonic measurement data as well
loadsonics = true
#which sonics
listsonics = ["kaijo", "t1irg", "t2lcsat"]#, "t2ucsat", "tjk"]
allowedsonics = ["kaijo", "t1irg", "t2irg", "t2lcsat", "t2ucsat", "tjk"]
#period to load from sonics
sonicsstart = DateTime(2021,05,31,12,00,00)
sonicsend = sonicsstart + Minute(10)
drtime = Second(600)

#location of sonic-input files (turbulence data)
tjkinfile =     joinpath(duerrpath, "scripts", "src", "data", "tjktmp.nc")
t1irginfile =   joinpath(duerrpath, "scripts", "src", "data", "t1irgtmp.nc")
t2irginfile =   joinpath(duerrpath, "scripts", "src", "data", "t2irgtmp.nc")
t2lcsatinfile = joinpath(duerrpath, "scripts", "src", "data", "t2lcsattmp.nc")
t2ucsatinfile = joinpath(duerrpath, "scripts", "src", "data", "t2ucsattmp.nc")
kaijoinfile =   joinpath(duerrpath, "scripts", "src", "data", "kaijotmp.nc")

#location of slow data input file
tjkslowfile = joinpath("/home/haugened/Documents/openfoam/duerr_les/", "scripts", "src", "data", "tjk_data.csv")
t2slowfile = joinpath("/home/haugened/Documents/openfoam/duerr_les/", "scripts", "src", "data", "t2slow.nc")
#################################################
#load additional functions if sonics should be included
if loadsonics
    if !@isdefined turb_loaded
        turb_loaded = false
    end
    if !turb_loaded
        include("/home/haugened/Documents/ibl_patch_snow/code/src/turb_data.jl")
        import .turb    
        turb_loaded = true
    end
end
#################################################
#input check
if !(plottype in allowedplottypes)
    println(string(plottype, " is not in ", allowedplottypes, ". Please set variable 'plottype' to allowed value. Exiting."))
    exit()
else
    println("Plot type: ", plottype)
end

if loadsonics
    for i in listsonics
        if !(i in allowedsonics)
            println(string(i, " is not in ", allowedsonics, ". Please set variable 'listsonics' to allowed value. Exiting."))
            exit()
        end
    end
    println("Loading sonic data from ", listsonics)
end

#################################################
#read in data
if !@isdefined forcingdata_read
    forcingdata_read = false
    forcingdata_read_sucessfully = false
end
if !forcingdata_read
    println("Reading forcing data...")
    if isfile(uforcingfile)
        (forcing_heights, forcing_times, forcing_u, forcing_T) = func.readforcingfromnetcdf(uforcingfile)
        forcingdata_read_sucessfully = true
        @info("Forcing data sucessfully read.")
    else
        @warn("Forcing data could not be read.")
    end
    forcingdata_read = true
end

if !@isdefined data_read
    data_read = false
end
if !data_read
    println("Reading cross section data...")
    (x, y, z, times, T, U) = func.readaerialprofilefromnetcdf(file)
    data_read = true
end

#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx
    transformx = false
end
if !transformx
    x = x .- snowtransition
    transformx = true
end

#create bottom line for plot
xtopo = func.createbottom(x, z)

#deal with sonic data
if loadsonics
    if !@isdefined sonics_read
        sonics_read = false
    end
    if !sonics_read
        #sonic heights
        sonic_heights = zeros(Float64, length(listsonics))
        println("Reading sonic data...")
        if "kaijo" in listsonics
            htmp = 0.3
            datatmp = func.readtotalturbasnetcdf(kaijoinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            @warn("Careful!! Rotate Kaijo data due to measurement position. Check pictures of setup!")
            new_w = datatmp.u
            new_u = datatmp.w
            datatmp.u = new_u
            datatmp.w = new_w
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            kaijodata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="kaijo", listsonics)] = htmp
        end
        if "t1irg" in listsonics
            htmp = 1.2
            datatmp = func.readtotalturbasnetcdf(t1irginfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t1irgdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t1irg", listsonics)] = htmp
        end
        if "t2irg" in listsonics
            htmp = 0.9
            datatmp = func.readtotalturbasnetcdf(t2irginfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2irgdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2irg", listsonics)] = htmp
        end
        if "t2lcsat" in listsonics
            htmp = 1.9
            datatmp = func.readtotalturbasnetcdf(t2lcsatinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2lcsatdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2lcsat", listsonics)] = htmp
        end
        if "t2ucsat" in listsonics
            htmp = 2.85
            datatmp = func.readtotalturbasnetcdf(t2ucsatinfile)
            datatmp.time .-= Hour(1)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2ucsatdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2ucsat", listsonics)] = htmp
        end
        if "tjk" in listsonics
            htmp = 5.0
            datatmp = func.readtotalturbasnetcdf(tjkinfile)
            #datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, reyavgtime)
            flxtmp = func.avgflux(flxtmp, reyavgtime)
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            tjkdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="tjk", listsonics)] = htmp
        end
        sonics_read = true
    end
end
##
fig = PyPlot.figure(figsize=(10,11))
gs = gridspec.GridSpec(7, 3, width_ratios=[5, 0.01, 5], height_ratios=[3,1,5,0.7,5,0.7,5])

subfig_letters = ["a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"]

for q in 1:size(quantities_to_plot, 1)
    quantity_to_plot = quantities_to_plot[q]
    #take quantity to plot
    if quantity_to_plot == "u"
        vals = U[:,:,1]
        x_label_text = L"u~\mathrm{[m~s^{-1}]}"
        quantity_symbol = :U0
    elseif quantity_to_plot == "v"
        vals = U[:,:,2]
        x_label_text = L"v~\mathrm{[m~s^{-1}]}"
        quantity_symbol = :U1
    elseif quantity_to_plot == "w"
        vals = U[:,:,3]
        x_label_text = L"w~\mathrm{[m~s^{-1}]}"
        quantity_symbol = :U2
    elseif quantity_to_plot == "T"
        vals = T .- 273.15
        x_label_text = L"T~\mathrm{[^\circ C]}"
        #x_label_text = L"T~\mathrm{[K]}"
        quantity_symbol = :T
    elseif quantity_to_plot == "wT"
        vals1 = T
        vals2 = U[:,:,3] #w
        x_label_text = L"\overline{w'T'}~\mathrm{[10^{-2}~K~m~s^{-1}]}"
        quantity_symbol = :wT
    elseif quantity_to_plot == "Tvar"
        vals = T
        x_label_text = L"\overline{\left( T' \right)^2}~\mathrm{[K^2]}"
        quantity_symbol = :Tvar
    elseif quantity_to_plot == "tke"
        vals1 = U[:,:,1] #u
        vals2 = U[:,:,2] #v
        vals3 = U[:,:,3] #w
        x_label_text = L"e~\mathrm{[m^2~s^{-2}]}"
        quantity_symbol = :tke
    elseif quantity_to_plot == "z/L"
        valsu = U[:,:,1] #u
        valsv = U[:,:,2] #v
        valsw = U[:,:,3] #w
        valsT = T
        x_label_text = L"\frac{z}{L}~\mathrm{[1]}"
        quantity_symbol = :zoverL
    end

    skiptoidx = findfirst(x->x==minimum(abs.(times .- skipsec)), abs.(times .- skipsec))
    println(string("Skipping ", skipsec, " seconds at the beginning."))

    if loadsonics
        if quantity_to_plot == "u"
            sonic_symbol = :u
        elseif quantity_to_plot == "v"
            sonic_symbol = :v
        elseif quantity_to_plot == "w"
            sonic_symbol = :w
        elseif quantity_to_plot == "T"
            sonic_symbol = :T
        elseif quantity_to_plot == "wT"
            sonic_symbol = :wT
        elseif quantity_to_plot == "Tvar"
            sonic_symbol = :TT
        elseif quantity_to_plot == "tke"
            sonic_symbol = :tke
        elseif quantity_to_plot == "z/L"
            sonic_symbol = :zoverL
        end
    end

    tjkmeteodata = func.csvtodataframe(tjkslowfile)
    t2vent = turb.readturbasnetcdf(t2slowfile, sonicsstart, sonicsend + Hour(2)) #t2 is in LT, not in UTC+1
    t2vent.time .-= Hour(1)
    height_t2vent = 1.00 #[m]
    tjkmeteodata_excerpt = tjkmeteodata[sonicsstart .<= tjkmeteodata.time .<= sonicsend, :]
    t2vent_excerpt = t2vent[sonicsstart .<= t2vent.time .<= sonicsend, :]

    meanTair_tjk = mean(tjkmeteodata_excerpt.tair_HygroVUE10)
    meanTair_t2vent = mean(t2vent_excerpt.vent_air_temp)

    scaledTair = t2vent_excerpt.vent_air_temp .* (meanTair_tjk/meanTair_t2vent)

    #choose profile lines
    profilelines = fill(false, size(x, 1), size(profilelocations, 1))
    sortperms = fill(-1, size(x, 1), size(profilelocations, 1))
    for i in 1:length(profilelocations)
        profilelines[:, i] = func.chooselistval(profilealong, profilelocations[i]-profhalfwidth, profilelocations[i]+profhalfwidth, hcat(x, y, z))
        tmp = func.sortprofileentries(profilealong, hcat(x, y, z), profilelines[:,i])
        sortperms[1:length(tmp), i] = tmp
    end

    firstall0idx = 0
    for i in 1:size(sortperms, 2)
        firstall0idx = maximum([firstall0idx, findlast(x->x!=-1, sortperms[:, i])+1])
    end
    sortperms = sortperms[1:firstall0idx-1, :]

    #extract profiles

    profilesu = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
    profilesv = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
    profilesw = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
    profilesT = fill(NaN, size(sortperms,1), length(profilelocations), size(T, 2))
    for i in 1:length(profilelocations)
        tmpu = U[profilelines[:, i],:,1]
        profilesu[1:size(tmpu,1),i,:] = tmpu[sortperms[1:size(tmpu, 1), i], :]
        tmpv = U[profilelines[:, i],:,2]
        profilesv[1:size(tmpu,1),i,:] = tmpv[sortperms[1:size(tmpv, 1), i], :]
        tmpw = U[profilelines[:, i],:, 3]
        profilesw[1:size(tmpu,1),i,:] = tmpw[sortperms[1:size(tmpw, 1), i], :]
        tmpT = T[profilelines[:, i],:]
        profilesT[1:size(tmpu,1),i,:] = tmpT[sortperms[1:size(tmpT, 1), i], :]
    end
    dist_along_profile = fill(NaN, size(sortperms, 1), length(profilelocations))
    if profilealong == "x"
        for i in 1:length(profilelocations)
            tmpdist = x[profilelines[:, i]]
            dist_along_profile[1:length(tmpdist), i] = tmpdist[sortperms[1:size(tmpdist,1), i]]
        end
    elseif profilealong == "z"
        for i in 1:length(profilelocations)
            tmpdist = z[profilelines[:, i]]
            dist_along_profile[1:length(tmpdist), i] = tmpdist[sortperms[1:size(tmpdist,1), i]]
        end
    end

    for i in 1:length(profilelocations)
        (profilesu[:,i,:], ) = func.reducedoubleentries(profilesu[:, i, :], dist_along_profile[:, i], profhalfwidth)
        (profilesv[:,i,:], ) = func.reducedoubleentries(profilesv[:, i, :], dist_along_profile[:, i], profhalfwidth)
        (profilesw[:,i,:], ) = func.reducedoubleentries(profilesw[:, i, :], dist_along_profile[:, i], profhalfwidth)
        (profilesT[:,i,:], dist_along_profile[:, i]) = func.reducedoubleentries(profilesT[:, i, :], dist_along_profile[:, i], profhalfwidth)    
        firstall0idx = 0
    end

    firstall0idx = 0
    for i in 1:size(profilelocations, 2)
        firstall0idx = maximum([firstall0idx, findfirst(x->isnan(x), dist_along_profile[:, i])])
    end
    profilesu = profilesu[1:firstall0idx, :, :]
    profilesv = profilesv[1:firstall0idx, :, :]
    profilesw = profilesw[1:firstall0idx, :, :]
    profilesT = profilesT[1:firstall0idx, :, :]
    dist_along_profile = dist_along_profile[1:firstall0idx, :]

    #get terrain height at profile locations
    xidcs = zeros(Int64, size(profilelocations, 1))
    meanheights = zeros(Float64, size(profilelocations, 1))
    for l in 1:length(xidcs)
        idcstmp = findall(x->profilelocations[l]-profhalfwidth .<= x .<= profilelocations[l]+profhalfwidth, xtopo.x)
        mintmp = minimum(abs.(xtopo.x .- profilelocations[l]))
        xidcs[l] = findfirst(x->x==mintmp, abs.(xtopo.x .- profilelocations[l]))
        meanheights[l] = mean(xtopo.h[idcstmp])
        dist_along_profile[:,l] .-= meanheights[l]
    end

    if quantity_to_plot == "z/L"
        #calculate fluxes
        profiles = zeros(Float64, size(profilesw, 1), size(profilesw, 2), size(profilesw, 3))

        timesteps = func.gettimesteps(times)
        meantimesteps = round(mean(timesteps), digits=3)
        movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
        
        for icol in 1:size(profilesw, 2)
            for jrow in 1:size(profilesw, 1)
                movavgprofu = func.movingaverage(profilesu[jrow, icol, :], movavgidcs)
                movavgprofv = func.movingaverage(profilesv[jrow, icol, :], movavgidcs)
                movavgprofw = func.movingaverage(profilesw[jrow, icol, :], movavgidcs)
                movavgprofT = func.movingaverage(profilesT[jrow, icol, :], movavgidcs)
                devu = profilesu[jrow, icol, :] .- movavgprofu
                devv = profilesv[jrow, icol, :] .- movavgprofv
                devw = profilesw[jrow, icol, :] .- movavgprofw
                devT = profilesT[jrow, icol, :] .- movavgprofT
                ustar = ((devu .* devw) .^2 .+ (devv .* devw) .^2) .^(1/4)
                wT_tmp = devw .* devT
                L_tmp = - (ustar .^ 3 .* movavgprofT) ./ (0.4*9.81 .* wT_tmp)
                if profilealong == "x" .|| profilealong == "y"
                    profiles[jrow, icol, :] = prof_locs_final[icol] ./ L_tmp
                else #profilealong == "z"
                    profiles[jrow, icol, :] = dist_along_profile[jrow] ./ L_tmp
                end
            end
        end

    elseif quantity_to_plot == "Tvar"
        #calculate fluxes
        profiles = zeros(Float64, size(profilesT, 1), size(profilesT, 2), size(profilesT, 3))

        timesteps = func.gettimesteps(times)
        meantimesteps = round(mean(timesteps), digits=3)
        movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
        
        for icol in 1:size(profilesT, 2)
            for jrow in 1:size(profilesT, 1)
                movavgprof1 = func.movingaverage(profilesT[jrow, icol, :], movavgidcs)
                dev1 = profilesT[jrow, icol, :] .- movavgprof1
                profiles[jrow, icol, :] = dev1 .* dev1
            end
        end
    elseif quantity_to_plot == "tke"
        #calculate fluxes
        profiles = zeros(Float64, size(profilesu, 1), size(profilesu, 2), size(profilesu, 3))

        timesteps = func.gettimesteps(times)
        meantimesteps = round(mean(timesteps), digits=3)
        movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
        
        for icol in 1:size(profilesu, 2)
            for jrow in 1:size(profilesu, 1)
                movavgprof1 = func.movingaverage(profilesu[jrow, icol, :], movavgidcs)
                movavgprof2 = func.movingaverage(profilesv[jrow, icol, :], movavgidcs)
                movavgprof3 = func.movingaverage(profilesw[jrow, icol, :], movavgidcs)
                dev1 = profilesu[jrow, icol, :] .- movavgprof1
                dev2 = profilesv[jrow, icol, :] .- movavgprof2
                dev3 = profilesw[jrow, icol, :] .- movavgprof3
                profiles[jrow, icol, :] = 0.5 .*(dev1 .^2 .+ dev2 .^2 .+ dev3 .^2)
            end
        end
    elseif quantity_to_plot == "wT" #fluxes
        #calculate fluxes
        profiles = zeros(Float64, size(profilesw, 1), size(profilesw, 2), size(profilesw, 3))

        timesteps = func.gettimesteps(times)
        meantimesteps = round(mean(timesteps), digits=3)
        movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
        
        for icol in 1:size(profilesw, 2)
            for jrow in 1:size(profilesw, 1)
                movavgprof1 = func.movingaverage(profilesw[jrow, icol, :], movavgidcs)
                movavgprof2 = func.movingaverage(profilesT[jrow, icol, :], movavgidcs)
                dev1 = profilesw[jrow, icol, :] .- movavgprof1
                dev2 = profilesT[jrow, icol, :] .- movavgprof2
                profiles[jrow, icol, :] = dev1 .* dev2
            end
        end
    elseif quantity_to_plot == "u"
        profiles = profilesu
    elseif quantity_to_plot == "v"
        profiles = profilesv
    elseif quantity_to_plot == "w"
        profiles = profilesw
    elseif quantity_to_plot == "T"
        profiles = profilesT .- 273.15
    end

    if postfilter
        println("-------------------------------------")
        println("Filtering of profiles activated.")
        println("No filtering of sonic measurement data!")
        println("Filter condition: ", filtercond)

        #get data at filter location
        post_take = func.chooseval3D(postfilterloc, hcat(x, y, z))

        println("location for filtering:")
        println(string("x = ", x[post_take], "m"))
        println(string("y = ", y[post_take], "m"))
        println(string("z = ", z[post_take], "m"))
        println("-------------------------------------")


        if postfilterquantity == "u"
            filterdata = U[post_take,:,1]
        elseif postfilterquantity == "v"
            filterdata = U[post_take,:,2]
        elseif postfilterquantity == "w"
            filterdata = U[post_take,:,3]
        elseif postfilterquantity == "T"
            filterdata = T[post_take,:]
        else
            @error("Wrong filter quantity! Please check.")
        end

        #calculate 1%-quantiles
        quantiles = collect(0.01:0.01:1)
        filterquantiles = quantile(filter(!isnan, filterdata), quantiles)

        #post filter condition
        postfiltercond = eval(Meta.parse(filtercond))
    else #no postfiltering
        postfiltercond = fill(true, size(profiles, 3))
        println("No post-filtering applied.")
    end

    profiles[:,:, .!(postfiltercond)] .= NaN

    medianprof = fill(NaN, size(profiles, 1), size(profiles, 2), 3)

    for j in 1:size(medianprof, 2)
        for i in 1:size(medianprof, 1)
            try
                (medianprof[i,j,1], medianprof[i,j,2], medianprof[i,j,3]) = quantile(filter(!isnan, profiles[i, j, skiptoidx:end]), [0.25, 0.5, 0.75])
            catch e
            end
        end
    end

    if quantity_to_plot == "T"
        #plot median and indicate profile locations
        medval = fill(NaN, size(vals, 1))
        for k in 1:size(vals, 1)
            medval[k] = median(filter(!isnan, vals[k, :]))
        end

        cmap = PyPlot.get_cmap("tab10")
        axa = fig.add_subplot(py"$(gs)[0,0:3]")
        #axa.set_title(string("Median T - Indicating profile locations"))
        hm = axa.tripcolor(x, z, medval, cmap="viridis", vmin=5, vmax=11)#maximum(medval)*0.9)
        axa.set_xlabel(L"x~[\mathrm{m}]")
        axa.set_ylabel(L"h~[\mathrm{m}]")
        axa.set_xlim(minimum(x), maximum(x))
        axa.set_xlim(-0.5, 12)
        axa.set_ylim(0,2)
        for k in 1:size(profilelocations, 1)
            l = lines.Line2D([profilelocations[k], profilelocations[k]], [meanheights[k], 5], color=cmap(k-1))
            axa.fill_betweenx(collect(0:0.2:5), profilelocations[k]-profhalfwidth, profilelocations[k]+profhalfwidth, color=cmap(k-1), alpha=0.5 )
        end
        topo = axa.plot(xtopo.x, xtopo.h, color="black", lw=3)
        axa.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
        PyPlot.gca().set_aspect("equal")
        #axacb = fig.add_subplot(gs[1,4])
        #PyPlot.axis("off")
        hm_cb = fig.colorbar(hm, ax=axa)
        hm_cb.set_label(L"T~[\mathrm{^\circ C}]")
        axa.text(-0.087, 1.01, subfig_letters[1], size="large", transform=axa.transAxes)
    end

    pltrow = 2*div(q-1, 2) + 3
    pltcol = 3 - 2*mod(q, 2)

    ax1 = fig.add_subplot(gs[pltrow, pltcol])
    #=if postfilter
        ax1.set_title(filtercond)
    else
        ax1.set_title("No filtering")
    end=#
    for k in 1:size(profilelocations, 1)
        ax1.plot(medianprof[:,k,2], dist_along_profile[:, k], label=string("x=", round(profilelocations[k], digits=1), "m"))#-profhalfwidth, " - ", profilelocations[k]+profhalfwidth, "m"))
        ax1.fill_betweenx(dist_along_profile[:, k], medianprof[:,k,1], medianprof[:,k,3], alpha=0.3)
    end
    #ax1.plot(medianprof[:,2,2], dist_along_profile[:, 2], label=string("x=", profilelocations[2]-profhalfwidth, " - ", profilelocations[2]+profhalfwidth, "m"))
    #ax1.fill_betweenx(dist_along_profile[:, 2], medianprof[:,2,1], medianprof[:,2,3], alpha=0.3)
    #ax1.plot(medianprof[:,3,2], dist_along_profile[:, 3], label=string("x=", profilelocations[3]-profhalfwidth, " - ", profilelocations[3]+profhalfwidth, "m"))
    #ax1.fill_betweenx(dist_along_profile[:, 3], medianprof[:,3,1], medianprof[:,3,3], alpha=0.3)
    #calculate mean and IQR for sonics
    if loadsonics .&& quantity_to_plot in ["u", "w", "tke"]
        if "kaijo" in listsonics
            sonicmed = quantile(filter(!isnan, kaijodata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="kaijo", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="short-path ultrasonic")
        end
        if "t1irg" in listsonics
            sonicmed = quantile(filter(!isnan, t1irgdata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t1irg", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="ECT1, h=1.2m")
        end
        if "t2irg" in listsonics
            sonicmed = quantile(filter(!isnan, t2irgdata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t2irg", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="ECT2, h=0.9m")
        end
        if "t2lcsat" in listsonics
            sonicmed = quantile(filter(!isnan, t2lcsatdata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t2lcsat", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="ECT2, h=1.9m")
        end
        if "t2ucsat" in listsonics
            sonicmed = quantile(filter(!isnan, t2ucsatdata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t2ucsat", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="T2UCSAT (0.25-0.75)")
        end
        if "tjk" in listsonics
            sonicmed = quantile(filter(!isnan, tjkdata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="tjk", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="TJK (0.25-0.75)")
        end
    end
    #if quantity_symbol == :T
    #    #t2vent_excerpt[:, :vent_air_temp]
    #    ventmed = quantile(filter(!isnan, scaledTair), [0.25,0.5,0.75])
    #    ax1.errorbar(ventmed[2], 1.0, xerr=[[abs(ventmed[1]-ventmed[2])], [abs(ventmed[3]-ventmed[2])]], fmt="o", label=L"T_{vent}~\mathrm{(scaled)}")
    #end
    ax1.set_xlabel(x_label_text)
    ax1.set_ylim(0,2.0)
    if quantity_to_plot == "u"
        ax1.set_xlim(-0.1, 3.5)
    elseif quantity_to_plot == "w"
        ax1.set_xlim(-0.15, 0.15)
    elseif quantity_to_plot == "T"
        ax1.set_xlim(0, 15)
    elseif quantity_to_plot =="wT"
        ax1.set_xlim(-0.04, 0.04)
    elseif quantity_to_plot == "tke"
        ax1.set_xlim(0,0.5)
    end
    ax1.grid()
    if q in [2, 4]
        ax1.tick_params(axis="y", labelleft=false)
        ax1.text(-0.1, 1.03, subfig_letters[q+1], size="large", transform=ax1.transAxes)
    else
        ax1.set_ylabel(L"h~\mathrm{[m]}")
        ax1.text(-0.15, 1.03, subfig_letters[q+1], size="large", transform=ax1.transAxes)
    end
    if quantity_to_plot == "u" #q==size(quantities_to_plot, 1)
        (h, l) = ax1.get_legend_handles_labels()
        println(h)
        println(l)
        ax_leg = fig.add_subplot(gs[7, 3])
        l1 = ax_leg.legend(h[1:size(profilelocations, 1)], l[1:size(profilelocations, 1)], loc="center left")#,borderaxespad=0)
        l2 = ax_leg.legend(h[size(profilelocations, 1)+1:end], l[size(profilelocations, 1)+1:end], loc="center right")#,borderaxespad=0)
        ax_leg.axis("off")
        PyPlot.gca().add_artist(l1)
        ax_leg.text(0.15, 0.78, "LES", size="medium", transform=ax_leg.transAxes)
        ax_leg.text(0.54, 0.78, "measurements", size="medium", transform=ax_leg.transAxes)
    end
end
#PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/med_profiles_3rd.jpg")
##
=#

####################################################################################
####################################################################################
###  Fig. XX:           MEDIAN BUOYANCY FLUXES AS HEAT MAP                       ###
####################################################################################
####################################################################################
#Figxx: Mean heat fluxes as heatmap
#todo: change to new data,
#      enable surface plot
#      change snowtransition to 8m
##
#=
using NCDatasets, DataFrames, StatsBase, LaTeXStrings, PyCall, Dates, ProgressMeter, CSV
import PyPlot
PyPlot.pygui(true)
GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")
animation = pyimport("matplotlib.animation")
numpy = pyimport("numpy")
cramericm = pyimport("cmcrameri.cm")

#outer directory of case
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/second_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#folder containing the time steps
readfromfolder = joinpath(casedir, "postProcessing/surfaces/yNormal")

#netCDF4-file containing the profile (length 3); use same time steps!!
file = joinpath.(readfromfolder, "yNormal_time8th.nc")

#value of x where bare to snow transition is
snowtransition = 8.0 #m

#plot median crosssection or animation
plotstyle = "median"
allowedplotstyle = ["median", "animation"]

#quantity for profiles 
quantity_to_plot = "wT"
allowedquantities = ["u", "v", "w", "T", "wT", "tke", "Tvar", "tke", "z/L"]

#only for flux-type variables
reyavgtime = Second(50)

#flag for filtering of profiles (only for plotstyle == "median")
postfilter = false
postfilterquantity = "u" #u, v, w, or T
#location of point to evaluate for post-filtering (fetch distance coordinates for x) [m]
postfilterloc = [-1.0, 2, 0.5]
#filter condition
filtercond = "0.0 .<= filterdata .<= filterquantiles[50]"

#seconds to skip (spinuptime); only for plotstyle == "median"
skipsec = 10

#Interval [ms] for plotting (only for plotstyle == "animation").
intervalmsec = 30

#calculated for video (only for plotstyle == "animation")
framespersec = 1000/intervalmsec

#flag if animation should be saved
anisave = true

#file name of video (only for plotstyle == "animation")
anistring = string(joinpath(casedir, "video_T_cut_to_screen.mp4"))

#################################################
#input check
if !(quantity_to_plot in allowedquantities)
    println(string(quantity_to_plot), " is not in ", allowedquantities, ". Please set variable 'quantity_to_plot' to allowed value. Exiting.")
    exit()
else
    println("Plotting cross section of ", quantity_to_plot)
end

if !(plotstyle in allowedplotstyle)
    println(string(plotstyle), " is not in ", allowedplotstyle, ". Please set variable 'plotstyle' to allowed value. Exiting.")
    exit()
else
    println("Plot type ", plotstyle)
end
#################################################
#read in data
if !@isdefined data_read
    data_read = false
end
if !data_read
    println("Reading cross section data...")
    (x, y, z, times, T, U) = func.readaerialprofilefromnetcdf(file)
    data_read = true
end

#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx
    transformx = false
end
if !transformx
    x = x .- snowtransition
    transformx = true
end

#create bottom line for plot
xtopo = func.createbottom(x, z)

#################################################
#take and/or quantity to plot
if quantity_to_plot == "u"
    vals = U[:,:,1]
    valmin = -2
    valmax = 2
    colmap = "viridis"
    cblabel = L"u~\mathrm{[m~s^{-1}]}"
    figtitle = string("u")
elseif quantity_to_plot == "v"
    vals = U[:,:,2]
    valmin = -1
    valmax = 2
    colmap = "viridis"
    cblabel = L"v~\mathrm{[m~s^{-1}]}"
    figtitle = string("v")
elseif quantity_to_plot == "w"
    vals = U[:,:,3]
    valmin = -1
    valmax = 1
    colmap = "viridis"
    cblabel = L"w~\mathrm{[m~s^{-1}]}"
    figtitle = string("w")
elseif quantity_to_plot == "T"
    vals = T
    valmin = 280
    valmax = 287
    colmap = "viridis"
    cblabel = L"T~\mathrm{[K]}"
    figtitle = string("T")
elseif quantity_to_plot == "wT"
    vals1 = T
    vals2 = U[:,:,3] #w
    valmin = -0.1
    valmax = -valmin
    colmap = cramericm.vik
    cblabel = L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}"
    figtitle = string("w'T'")

    #calculate fluxes
    vals = similar(vals1)
    fill!(vals, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
    
    @showprogress "Calculating fluxes" for jrow in 1:size(vals1, 1)
        movavg1 = func.movingaverage(vals1[jrow, :], movavgidcs)
        movavg2 = func.movingaverage(vals2[jrow, :], movavgidcs)
        dev1 = vals1[jrow, :] .- movavg1
        dev2 = vals2[jrow, :] .- movavg2
        vals[jrow, :] = dev1 .* dev2
    end
elseif quantity_to_plot == "Tvar"
    vals1 = T
    valmin = 0
    valmax = 2.0
    colmap = "viridis"
    cblabel = L"\overline{\left( T' \right)^2}~\mathrm{[K^2]}"
    figtitle = string("Tvar")

    #calculate fluxes
    vals = similar(vals1)
    fill!(vals, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
    
    @showprogress "Calculating fluxes" for jrow in 1:size(vals1, 1)
            movavg1 = func.movingaverage(vals1[jrow, :], movavgidcs)
            dev1 = vals1[jrow, :] .- movavg1
            vals[jrow, :] = dev1 .* dev1
    end
elseif quantity_to_plot == "tke"
    vals1 = U[:,:,1] #u
    vals2 = U[:,:,2] #v
    vals3 = U[:,:,3] #w
    vmin = 0.0
    vmax = 2.0
    colmap = "viridis"
    cblabel = L"e~\mathrm{[m^2~s^{-2}]}"
    figtitle = string("TKE")

    #calculate fluxes
    vals = similar(vals1)
    fill!(vals, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
    
    @showprogress "Calculating fluxes" for jrow in 1:size(vals1, 1)
            movavg1 = func.movingaverage(vals1[jrow, :], movavgidcs)
            movavg2 = func.movingaverage(vals2[jrow, :], movavgidcs)
            movavg3 = func.movingaverage(vals3[jrow, :], movavgidcs)
            dev1 = vals1[jrow, :] .- movavg1
            dev2 = vals2[jrow, :] .- movavg2
            dev3 = vals3[jrow, :] .- movavg3
            vals[jrow, :] = 0.5 .*(dev1 .^2 .+ dev2 .^2 .+ dev3 .^2)
    end
elseif quantity_to_plot == "z/L"
    valsu = U[:,:,1] #u
    valsv = U[:,:,2] #v
    valsw = U[:,:,3] #w
    valsT = T[:,:]
    valmin = -0.1
    valmax = -valmin
    colmap = cramericm.vik
    cblabel = L"\frac{z}{L}~\mathrm{[1]}"
    figtitle = string("z/L")

    #calculate fluxes
    vals = similar(valsu)
    fill!(vals, NaN)

    timesteps = func.gettimesteps(times)
    meantimesteps = round(mean(timesteps), digits=3)
    movavgidcs = round(Int, reyavgtime/Millisecond(meantimesteps*1e3))
    
    @showprogress "Calculating fluxes" for jrow in 1:size(valsu, 1)
            movavgu = func.movingaverage(valsu[jrow, :], movavgidcs)
            movavgv = func.movingaverage(valsv[jrow, :], movavgidcs)
            movavgw = func.movingaverage(valsw[jrow, :], movavgidcs)
            movavgT = func.movingaverage(valsT[jrow, :], movavgidcs)
            devu = valsu[jrow, :] .- movavgu
            devv = valsv[jrow, :] .- movavgv
            devw = valsw[jrow, :] .- movavgw
            devT = valsT[jrow, :] .- movavgT
            ustar = ((devu .* devw) .^2 .+ (devv .* devw) .^2) .^(1/4)
            wT_tmp = devw .* devT
            L_tmp = - (ustar .^ 3 .* movavgT) ./ (0.4*9.81 .* wT_tmp)
            vals[jrow, :] = z[jrow] ./ L_tmp
    end
end

#timesround = round.(times, digits=3)
animidcs = func.getanimidcs(times, intervalmsec)

if plotstyle == "median"
    skiptoidx = findfirst(x->x==minimum(abs.(times .- skipsec)), abs.(times .- skipsec))
    println(string("Skipping ", skipsec, " seconds at the beginning."))

    if postfilter
        println("-------------------------------------")
        println("Filtering of profiles activated.")
        println("Filter condition: ", filtercond)
    
        #get data at filter location
        take_x = func.chooseval(postfilterloc[1], x);
        take_y = func.chooseval(postfilterloc[2], y);
        take_z = func.chooseval(postfilterloc[3], z);
    
        println("location for filtering:")
        println(string("x = ", x[take_x], "m"))
        println(string("y = ", y[take_y], "m"))
        println(string("z = ", z[take_z], "m"))
        println("-------------------------------------")
    
    
        if postfilterquantity == "u"
            filterdata = U[take_x,take_y,take_z,:,1]
        elseif postfilterquantity == "v"
            filterdata = U[take_x,take_y,take_z,:,2]
        elseif postfilterquantity == "w"
            filterdata = U[take_x,take_y,take_z,:,3]
        elseif postfilterquantity == "T"
            filterdata = T[take_x,take_y,take_z,:]
        else
            @error("Wrong filter quantity! Please check.")
        end
    
        #calculate 1%-quantiles
        quantiles = collect(0.01:0.01:1)
        filterquantiles = quantile(filter(!isnan, filterdata), quantiles)
    
        #post filter condition
        postfiltercond = eval(Meta.parse(filtercond))
        figtitle = string(figtitle, "; ", filtercond)
    else #no postfiltering
        postfiltercond = fill(true, size(vals, 3))
        println("No post-filtering applied.")
    end
    
    vals[:,:, .!(postfiltercond)] .= NaN
    
    medianvals = fill(NaN, size(vals, 1))
    for jrow in 1:size(vals, 1)
        medianvals[jrow] = median(filter(!isnan, vals[jrow, skiptoidx:end]))
    end

    fig = PyPlot.figure(figsize=(13,5))
    ax1 = fig.add_subplot(111)
    ax1.set_title(string("Median ", figtitle, ", skipping ", skipsec, "s (spin up)"))
    hm = ax1.tripcolor(x, z, medianvals, cmap=colmap, vmin=valmin, vmax=valmax)
    topo = ax1.plot(xtopo.x, xtopo.h, color="black", lw=3)
    ax1.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
    ax1.set_xlabel(L"x~[\mathrm{m}]")
    ax1.set_ylabel(L"z~[\mathrm{m}]")
    ax1.set_xlim(minimum(x), maximum(x))
    ax1.set_xlim(-0.5,6)
    ax1.set_ylim(0,2)
    hm_cb = fig.colorbar(hm, ax=ax1)
    hm_cb.set_label(cblabel)
    PyPlot.gca().set_aspect("equal")
    
elseif plotstyle == "animation"
    "Update-function for the animation"
    """
    animupdate(frame::Int64)
    """
    function animupdate(frame::Int64)
        ax1.set_title(string(figtitle, "; LES t = ", round(times[frame], digits=1), " s"))
        newdata = vals[:,frame]
        #newdata = newdata[1:end-1]
        #hm.set_array(numpy.ravel(newdata))
        hm.set_array(newdata)
        if mod(frame,100) == 0
            println(string(frame, "/", animidcs[end]))
        end
        return hm
    end

    fig = PyPlot.figure(figsize=(23,7))
    ax1 = fig.add_subplot(111)
    hm = ax1.tripcolor(x, z, vals[:, 100], cmap=colmap, vmin=valmin, vmax=valmax, shading="gouraud")
    topo = ax1.plot(xtopo.x, xtopo.h, color="black", lw=3)
    ax1.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
    ax1.set_xlim(minimum(x), maximum(x))
    ax1.set_ylim(0, maximum(z))
    ax1.set_xlabel(L"x~[\mathrm{m}]")
    ax1.set_ylabel(L"z~[\mathrm{m}]")
    ax1.set_xlim(-0.5, 6)
    ax1.set_ylim(0, 2)
    hm_cb = fig.colorbar(hm, ax=ax1)
    hm_cb.set_label(cblabel)
    PyPlot.gca().set_aspect("equal")
    ani = animation.FuncAnimation(fig, animupdate, frames=animidcs[1:end], interval=intervalmsec, repeat=false)#, blit=false)
    if anisave
        ani.save(anistring, fps=framespersec)
    end
end
#PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/median_heatmap.jpg")
=#
##
####################################################################################
####################################################################################
###  Fig. XX:              COMPARING TEMPERATURE PROFILES                        ###
####################################################################################
####################################################################################
#=
using PyCall, Statistics, LaTeXStrings, NCDatasets, DataFrames, StatsBase, PyCall, Dates, ProgressMeter, CSV
import PyPlot
gridspec = pyimport("matplotlib.gridspec")
mpimg = pyimport("matplotlib.image")
animation = pyimport("matplotlib.animation")
numpy = pyimport("numpy")
cramericm = pyimport("cmcrameri.cm")
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
patches = pyimport("matplotlib.patches")
if gethostname() == "Michi-T450s"
    importdir = "/home/michi/Documents/slf/ibl_patch_snow/code/"
    pathtofile = "/home/michi/Documents/slf/data/"
elseif gethostname() == "LINUX24"
    importdir = "/home/haugened/Documents/ibl_patch_snow/code/"
    pathtofile = "/home/haugened/Documents/data/ir/"
end
include(joinpath(importdir, "src", "ir_evaluation.jl"))
import .irev
PyPlot.pygui(true)

#outer directory of case
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/third_run/"
#"/home/haugened/Documents/openfoam/duerr_les/"

#load functions from src/functions.jl
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#folder containing the time steps
readfromfolder = joinpath(casedir, "postProcessing/surfaces/yNormal")

#netCDF4-file containing the profile (length 3); use same time steps!!
file = joinpath.(readfromfolder, "yNormal_x5cm.nc")

#value of x where bare to snow transition is
snowtransition = 8.0 #m

#quantity for profiles 
quantity_to_plot = "T"

#profile locations (x-value) for comparison
profilelocations = [0.4, 2.0, 5.0]
profhalfwidth = 0.05 #m

#if false, read already produced profile (profilelocations, profilehalfwidht get overwritten!)
#if true, generate new ones (might need to execute line by line so it doesn't crash)
IRreadnew = false
#IRprofsloc = "/home/haugened/Documents/data/ir/210531_143457/les_paper_24_IR_profiles.nc"
IRfinalprofsloc = "/home/haugened/Documents/data/ir/210531_115756/les_paper_24_IR_final_profiles.nc"

#read in data
if !@isdefined data_read
    data_read = false
end
if !data_read
    println("Reading cross section data...")
    (x, y, z, times, T, U) = func.readaerialprofilefromnetcdf(file)
    T .-= 273.15
    data_read = true
end

#transform x-coordinates to new fetch distance coordinates (x=0 at transition from bare to snow)
if !@isdefined transformx
    transformx = false
end
if !transformx
    x = x .- snowtransition
    transformx = true
end

#create bottom line for plot
xtopo = func.createbottom(x, z)

#################################################
#choose profile lines for LES
profilelines = fill(false, size(x, 1), size(profilelocations, 1))
sortperms = fill(-1, size(x, 1), size(profilelocations, 1))
for i in 1:length(profilelocations)
    profilelines[:, i] = func.chooselistval("z", profilelocations[i]-profhalfwidth, profilelocations[i]+profhalfwidth, hcat(x, y, z))
    tmp = func.sortprofileentries("z", hcat(x, y, z), profilelines[:,i])
    sortperms[1:length(tmp), i] = tmp
end

firstall0idx = 0
for i in 1:size(sortperms, 2)
    firstall0idx = maximum([firstall0idx, findlast(x->x!=-1, sortperms[:, i])+1])
end
sortperms = sortperms[1:firstall0idx-1, :]

#extract profiles

#profilesu = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
#profilesv = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
#profilesw = fill(NaN, size(sortperms,1), length(profilelocations), size(U, 2))
profilesT = fill(NaN, size(sortperms,1), length(profilelocations), size(T, 2))
dist_along_profile = fill(NaN, size(sortperms, 1), length(profilelocations))
for i in 1:length(profilelocations)
    #tmpu = U[profilelines[:, i],:,1]
    #profilesu[1:size(tmpu,1),i,:] = tmpu[sortperms[1:size(tmpu, 1), i], :]
    #tmpv = U[profilelines[:, i],:,2]
    #profilesv[1:size(tmpu,1),i,:] = tmpv[sortperms[1:size(tmpv, 1), i], :]
    #tmpw = U[profilelines[:, i],:, 3]
    #profilesw[1:size(tmpu,1),i,:] = tmpw[sortperms[1:size(tmpw, 1), i], :]
    tmpT = T[profilelines[:, i],:]
    profilesT[1:size(tmpT,1),i,:] = tmpT[sortperms[1:size(tmpT, 1), i], :]
    tmpdist = z[profilelines[:, i]]
    dist_along_profile[1:length(tmpdist), i] = tmpdist[sortperms[1:size(tmpdist,1), i]]

end

for i in 1:length(profilelocations)
    #(profilesu[:,i,:], ) = func.reducedoubleentries(profilesu[:, i, :], dist_along_profile[:, i], profhalfwidth)
    #(profilesv[:,i,:], ) = func.reducedoubleentries(profilesv[:, i, :], dist_along_profile[:, i], profhalfwidth)
    #(profilesw[:,i,:], ) = func.reducedoubleentries(profilesw[:, i, :], dist_along_profile[:, i], profhalfwidth)
    (profilesT[:,i,:], dist_along_profile[:, i]) = func.reducedoubleentries(profilesT[:, i, :], dist_along_profile[:, i], profhalfwidth)    
    firstall0idx = 0
end

firstall0idx = 0
for i in 1:size(profilelocations, 2)
    firstall0idx = maximum([firstall0idx, findfirst(x->isnan(x), dist_along_profile[:, i])])
end
#profilesu = profilesu[1:firstall0idx, :, :]
#profilesv = profilesv[1:firstall0idx, :, :]
#profilesw = profilesw[1:firstall0idx, :, :]
profilesT = profilesT[1:firstall0idx, :, :]
dist_along_profile = dist_along_profile[1:firstall0idx, :]

#get terrain height at profile locations
xidcs = zeros(Int64, size(profilelocations, 1))
meanheights = zeros(Float64, size(profilelocations, 1))
for l in 1:length(xidcs)
    idcstmp = findall(x->profilelocations[l]-profhalfwidth .<= x .<= profilelocations[l]+profhalfwidth, xtopo.x)
    mintmp = minimum(abs.(xtopo.x .- profilelocations[l]))
    xidcs[l] = findfirst(x->x==mintmp, abs.(xtopo.x .- profilelocations[l]))
    meanheights[l] = mean(xtopo.h[idcstmp])
    dist_along_profile[:,l] .-= meanheights[l]
end

#################################################
#take and/or quantity to plot
vals = T
vals .-= 273.15 #K to C
valmin = 7
valmax = 15
colmap = "viridis"
cblabel = L"T~\mathrm{[^\circ C]}"

#functions
function loadfromNetCDF4(file::String, varname::String, startidx::Int, endidx::Int)
    netfile = NCDataset(file, "r")
    dat = netfile[varname]
    if endidx == -1
        dattot = dat[:,:,startidx:end]
    else
        dattot = dat[:,:,startidx:endidx]
    end
    close(netfile)
    return dattot
end

function cutIRdata(IRdatain::Array)
    IRdatain[:,443:540,:] .= NaN
    IRdatain = IRdatain[11:end,19:end,:]
    return IRdatain
end

#definitions
irfile =  "/home/haugened/Documents/data/ir/210531_115756/irdata_03_to_05.nc"
irfile2 = "/home/haugened/Documents/data/ir/210531_115756/irdata_06_to_08.nc"
irfile3 = "/home/haugened/Documents/data/ir/210531_115756/irdata_09_to_11.nc"
cmperpxl = 0.6
cmperpxlu = 0.606
cmperpxlw = 0.549
framespersec = 29.9

IRdata = loadfromNetCDF4(irfile, "irdata", 1000, 1001)
#IRdata = IRdata[:,:,2000:end]
#cut
IRdata = cutIRdata(IRdata)
(IRdata, snowsurf) = irev.setsnowsurfacetonan(IRdata, -0.2, framewise=true)
irsample = IRdata[:,:,1]

if IRreadnew
    IRdata1 = loadfromNetCDF4(irfile, "irdata", 1, -1)

    #cut
    IRdata1 = cutIRdata(IRdata1)
    (IRdata1, snowsurf1) = irev.setsnowsurfacetonan(IRdata1, -0.2, framewise=true)

    #######################################
    #find profile pxls
    xis0pxl = 5

    profsleftcol = zeros(Int64, length(profilelocations))
    profsrightcol = zeros(Int64, length(profilelocations))
    profhalfwidthpxl = ceil(Int, (profhalfwidth*100)/cmperpxlu)

    for i in 1:length(profilelocations)
        tmp = abs.(collect(0:(size(IRdata1, 2)).-xis0pxl)/cmperpxlu .- profilelocations[i]*100)
        tmpmin = minimum(tmp)
        profmiddle = findfirst(x->x==tmpmin, tmp)
        profsleftcol[i] = profmiddle - profhalfwidthpxl
        profsrightcol[i] = profmiddle + profhalfwidthpxl
    end
    ########################################
    maxprofwidth = maximum(profsrightcol .- profsleftcol)+1
    IR1prof = fill(NaN, size(IRdata1, 1), maxprofwidth, length(profilelocations), size(IRdata1 ,3))
    snowsurf1prof = fill(0, size(IR1prof, 2), length(profilelocations), size(IR1prof, 4))

    for i in 1:length(profilelocations)
        IR1prof[:, :,  i, :] = IRdata1[:, profsleftcol[i]:profsrightcol[i], :]
        snowsurf1prof[:,i,:] = snowsurf1[profsleftcol[i]:profsrightcol[i], :]
    end


    IRdata1=nothing

    IRdata2 = loadfromNetCDF4(irfile2, "irdata", 1, -1)
    IRdata2 = cutIRdata(IRdata2)
    (IRdata2, snowsurf2) = irev.setsnowsurfacetonan(IRdata2, -0.2, framewise=true)

    IR2prof = fill(NaN, size(IRdata2, 1), maxprofwidth, length(profilelocations), size(IRdata2 ,3))
    snowsurf2prof = fill(0, size(IR2prof, 2), length(profilelocations), size(IR2prof, 4))
    for i in 1:length(profilelocations)
        IR2prof[:, :,  i, :] = IRdata2[:, profsleftcol[i]:profsrightcol[i], :]
        snowsurf2prof[:,i,:] = snowsurf2[profsleftcol[i]:profsrightcol[i], :]
    end

    IRdata2=nothing

    IRdata3 = loadfromNetCDF4(irfile3, "irdata", 1, 1800)
    IRdata3 = cutIRdata(IRdata3)
    (IRdata3, snowsurf3) = irev.setsnowsurfacetonan(IRdata3, -0.2, framewise=true)

    IR3prof = fill(NaN, size(IRdata3, 1), maxprofwidth, length(profilelocations), size(IRdata3 ,3))
    snowsurf3prof = fill(0, size(IR3prof, 2), length(profilelocations), size(IR3prof, 4))

    for i in 1:length(profilelocations)
        IR3prof[:, :,  i, :] = IRdata3[:, profsleftcol[i]:profsrightcol[i], :]
        snowsurf3prof[:,i,:] = snowsurf3[profsleftcol[i]:profsrightcol[i], :]
    end

    IRdata3=nothing

    profIR = cat(IR1prof, IR2prof, IR3prof, dims = 4)
    snowsurfprof = cat(snowsurf1prof, snowsurf2prof, snowsurf3prof, dims=3)

    IR1prof = nothing
    IR2prof = nothing
    IR3prof = nothing
    snowsurf1prof = nothing
    snowsurf2prof = nothing
    snowsurf3prof = nothing

    for i in 1:length(profilelocations)
        for j in 1:size(profIR, 2)
            for k in 1:size(profIR, 4)
                profIR[snowsurfprof[j, i, k]:end, j, i, k] .= NaN
            end
        end
    end 

    profIR = profIR[1:firstallNaNrow-1, :,:,:]
    #=
    function saveIRprofilesasnetcdf(location::String, profiles::Array, centers::Vector, halfwidth::Number, deflatelvl::Int64=5)
        @info("Saving IR profiles to NetCDF4-file")
        ds = NCDataset(location, "c")
        defDim(ds, "rows", size(profiles, 1))
        defDim(ds, "cols", size(profiles, 2))
        defDim(ds, "nr profiles", size(profiles, 3))
        defDim(ds, "frames", size(profiles, 4))
        defDim(ds, "hwidth", 1)
        defVar(ds, "IRprof", profiles, ("rows", "cols", "nr profiles", "frames"); shuffle=true, deflatelevel=deflatelvl)
        defVar(ds, "center", centers, ("nr profiles",); shuffle=true, deflatelevel=deflatelvl)
        defVar(ds, "halfwidth", halfwidth, ("hwidth",); shuffle=true, deflatelevel=deflatelvl)
        close(ds)
    end
    =#
    #saveIRprofilesasnetcdf(IRprofsloc, profIR, profilelocations, profhalfwidth)

    #generate heights for IR profiles
    #first: collapse profiles along width
    IRsingleprofs = fill(NaN, size(profIR, 1), size(profIR, 3), size(profIR, 4))
    for k in 1:size(IRsingleprofs, 3)
        for j in 1:size(IRsingleprofs, 2)
            for i in 1:size(IRsingleprofs, 1)
                try
                    IRsingleprofs[i, j, k] = median(filter(!isnan, profIR[i, :, j, k]))
                catch e
                end
            end
        end
    end

    IRheightsraw = fill(NaN, size(IRsingleprofs, 1), size(IRsingleprofs, 2), size(IRsingleprofs, 3))

    for i in 1:size(IRheightsraw, 3)
        for j in 1:size(IRheightsraw, 2)
            snowsurftmp = findfirst(x->isnan(x), IRsingleprofs[:, j, i])
            if isnothing(snowsurftmp)
                snowsurftmp = size(IRsingleprofs, 1)
            end
            pxlheighttmp = snowsurftmp .- collect(1:snowsurftmp-1)
            IRheightsraw[1:length(pxlheighttmp), j, i] = pxlheighttmp .* (cmperpxlw/100)
        end
    end

    maxIRheightraw = maximum(filter(!isnan, IRheightsraw))
    IRheights = collect(0.01:0.01:maxIRheightraw)

    IRsingleprofs_final = fill(NaN, size(IRheights, 1), size(IRsingleprofs, 2), size(IRsingleprofs, 3))

    using ProgressMeter

    @showprogress for k in 1:size(IRsingleprofs_final, 3)
        for j in 1:size(IRsingleprofs_final, 2)
            for i in 1:size(IRheights, 1)
                if i == 1
                    hlow = 0.0
                    hhigh = IRheights[1] + ((IRheights[2] - IRheights[1]) / 2.0)
                elseif i==size(IRheights, 1)
                    hlow = IRheights[end-1] + ((IRheights[end] - IRheights[end-1]) / 2.0)
                    hhigh = IRheights[end]+1.0
                else
                    hlow = IRheights[i-1] + ((IRheights[i] - IRheights[i-1]) / 2.0)
                    hhigh = IRheights[i] + ((IRheights[i+1] - IRheights[i]) / 2.0)
                end
                cond = hlow .< IRheightsraw[:, j, i] .<= hhigh
                if count(cond) == 0
                    IRsingleprofs_final[i, j, k] = NaN
                else
                    vals = filter(!isnan, IRsingleprofs[cond, j, k])
                    if size(vals, 1) == 0
                        IRsingleprofs_final[i, j, k] = NaN
                    else
                        IRsingleprofs_final[i, j, k] = median(vals)
                    end
                end
            end
        end
    end

    function savefinalIRprofilesasnetcdf(location::String, profiles::Array, heights::Vector, centers::Vector, deflatelvl::Int64=5)
        @info("Saving IR profiles to NetCDF4-file")
        ds = NCDataset(location, "c")
        defDim(ds, "rows", size(profiles, 1))
        defDim(ds, "nr profiles", size(profiles, 2))
        defDim(ds, "frames", size(profiles, 3))
        defDim(ds, "hwidth", 1)
        defVar(ds, "IRprof", profiles, ("rows", "nr profiles", "frames"); shuffle=true, deflatelevel=deflatelvl)
        defVar(ds, "height", heights, ("rows",); shuffle=true, deflatelevel=deflatelvl)
        defVar(ds, "center", centers, ("nr profiles",); shuffle=true, deflatelevel=deflatelvl)
        close(ds)
    end

    savefinalIRprofilesasnetcdf(IRfinalprofsloc, IRsingleprofs_final, IRheights, profilelocations)

    IRread = true
else
    #=
    function readIRprofilesfromnetcdf(location::String)
        netfile = NCDataset(location, "r")
        IRtmp = netfile["IRprof"]
        centertmp = netfile["center"]
        hwtmp = netfile["halfwidth"]
        IRprof = IRtmp[:,:,:,:]
        center = centertmp[:]
        hw = hwtmp[:]
        return IRprof, center, hw
    end

    (profIR, profilelocations, profhalfwidth) = readIRprofilesfromnetcdf(IRprofsloc)
    =#
    function readfinalIRprofilesfromnetcdf(location::String)
        netfile = NCDataset(location, "r")
        IRtmp = netfile["IRprof"]
        centertmp = netfile["center"]
        htmp = netfile["height"]
        IRprof = IRtmp[:,:,:,:]
        center = centertmp[:]
        h = htmp[:]
        return IRprof, center, h
    end

    (IRsingleprofs_final, profilelocations, IRheights) = readfinalIRprofilesfromnetcdf(IRfinalprofsloc)
end

########################################
#median and IQR for profiles

les_profiles_med = fill(NaN, size(profilesT, 1), length(profilelocations), 3)
ir_profiles_med = fill(NaN, size(IRsingleprofs_final, 1), length(profilelocations), 3)

for i in 1:size(profilelocations, 1)
    for j in 1:size(les_profiles_med, 1)
        try
            les_profiles_med[j, i, :] = quantile(filter(!isnan, profilesT[j, i, :]), [0.25, 0.5, 0.75])
        catch e
        end
    end
    for j in 1:size(IRsingleprofs_final, 1)
        try
            ir_profiles_med[j, i, :] = quantile(filter(!isnan, IRsingleprofs_final[j, i, :]), [0.25, 0.5, 0.75])
        catch e
        end
    end
end

les_plot_profile_1 = les_profiles_med[:, 1, :]
les_plot_profile_2 = les_profiles_med[:, 2, :]
les_plot_profile_3 = les_profiles_med[:, 3, :]

ir_plot_profile_1 = IRsingleprofs_final[:,1, :] #ir_profiles_med[:, 1, :]
ir_plot_profile_2 = IRsingleprofs_final[:,2, :] #ir_profiles_med[:, 2, :]
ir_plot_profile_3 = IRsingleprofs_final[:,3, :] #ir_profiles_med[:, 3, :]

#=
firstNaNs = fill(size(ir_profiles_med, 2),3)
for i in 1:size(firstNaNs, 1)
    try
        firstNaNs[i] = findfirst(x->isnan(x), ir_profiles_med[:,i,2])
    catch e
        firstNaNs[i] = size(ir_profiles_med, 1)
    end
end

dist_along_IR = fill(NaN, size(profIR, 1), length(profilelocations))
for i in 1:length(profilelocations)
    dist_along_IR[1:firstNaNs[i], i] = (firstNaNs[i] .- collect(1:firstNaNs[i])) .* (cmperpxlw/100)
end
=#

##
cmap = PyPlot.get_cmap("tab10")
xis0pxl = 5

fig = PyPlot.figure(figsize=(10,7))
gs = gridspec.GridSpec(5, 3, height_ratios=[3,1,6,1,1])
ys = collect(0:0.01:2)

ax0 = fig.add_subplot(py"$(gs)[0,0:3]")
xleft = -(xis0pxl-1)*cmperpxlu/100
xright = size(IRdata, 2)*cmperpxlu/100+xleft
ybottom = 0
ytop = size(IRdata, 1) * cmperpxlw/100
hm1 = ax0.imshow(irsample, extent = [xleft,xright,ybottom,ytop], cmap=cramericm.batlow, vmin=4, vmax=8)
#ax0.set_ylabel(L"h~\mathrm{[m]}")#, fontsize=fntsze)
ax0.tick_params(axis="y", labelleft=false)
ax0.set_xlabel(L"x~\mathrm{[m]}")
#ax1.tick_params(labelsize=lblsze)
#ax0.tick_params(axis="x", labelbottom=false)
ax0.axvline(0.4, 0, 1, color=cmap(3), ls="-")#, label="time step taken")
ax0.fill_betweenx(ys, fill(-profhalfwidth[1] + profilelocations[1], length(ys)), fill(profhalfwidth[1] + profilelocations[1], length(ys)), color=cmap(3), alpha=0.5)
ax0.axvline(2.0, 0, 1, color=cmap(4), ls=":")#, label="time step taken")
ax0.fill_betweenx(ys, fill(-profhalfwidth[1] + profilelocations[2], length(ys)), fill(profhalfwidth[1] + profilelocations[2], length(ys)), color=cmap(4), alpha=0.5)
ax0.axvline(5.0, 0, 1, color=cmap(5), ls="-.")#, label="time step taken")
ax0.fill_betweenx(ys, fill(-profhalfwidth[1] + profilelocations[3], length(ys)), fill(profhalfwidth[1] + profilelocations[3], length(ys)), color=cmap(5), alpha=0.5)
ax0.set_ylim(0,2)
ax0.axhline()
IRwindarrow = ax0.annotate("", xy=(0.8,2.15), xytext=(0,2.15), annotation_clip=false, arrowprops=Dict("lw" => 2, "arrowstyle" => "->", "facecolor" => "black"))
IRwindarrowtxt = ax0.text(0, 2.21, "wind")

ax = fig.add_subplot(gs[3,1])
ax.set_title(L"x=0.4m")
fig.add_artist(patches.Rectangle((-0.19, -0.2), 1.255, 1.32, edgecolor=cmap(3), transform=ax.transAxes, linewidth=3, fill=false, alpha=0.8))
ax.plot(les_plot_profile_1[:, 2], dist_along_profile[:,1])
ax.fill_betweenx(dist_along_profile[:,1], les_plot_profile_1[:, 1], les_plot_profile_1[:, 3], alpha=0.5 )
ax.plot(ir_plot_profile_1[:,2], IRheights)
ax.fill_betweenx(IRheights, ir_plot_profile_1[:,1], ir_plot_profile_1[:,3], alpha=0.5 )
ax.set_xlabel(L"T~[\mathrm{^\circ C}]")
ax.set_ylabel(L"h~[\mathrm{m}]")
ax.grid()
ax.set_xlim(0,13)
ax.set_ylim(0,2)
#ax.tick_params(labelsize=lblsze)
#ax.tick_params(axis="x", labelbottom=false)

ax2 = fig.add_subplot(gs[3,2], sharey=ax)
ax2.set_title(L"x=2.0m")
fig.add_artist(patches.Rectangle((-0.1, -0.2), 1.16, 1.32, edgecolor=cmap(4), transform=ax2.transAxes, linewidth=3, ls=":", fill=false, alpha=0.8))
ax2.plot(les_plot_profile_2[:,2], dist_along_profile[:,2])
ax2.fill_betweenx(dist_along_profile[:,2], les_plot_profile_2[:,1], les_plot_profile_2[:,3], alpha=0.5 )
ax2.plot(ir_plot_profile_2[:,2], IRheights)
ax2.fill_betweenx(IRheights, ir_plot_profile_2[:,1], ir_plot_profile_2[:,3], alpha=0.5 )
ax2.set_xlabel(L"T~[\mathrm{^\circ C}]")
#ax2.set_ylabel(L"h~[\mathrm{m}]")
ax2.grid()
ax2.set_xlim(0,13)
ax2.tick_params(axis="y", labelleft=false)

ax3 = fig.add_subplot(gs[3,3], sharex=ax)
ax3.set_title(L"x=5.0m")
fig.add_artist(patches.Rectangle((-0.1, -0.2), 1.16, 1.32, edgecolor=cmap(5), transform=ax3.transAxes, linewidth=3, ls="-.", fill=false, alpha=0.8))
ax3.plot(les_plot_profile_3[:,2], dist_along_profile[:,3], label="LES")
ax3.fill_betweenx(dist_along_profile[:,3], les_plot_profile_3[:,1], les_plot_profile_3[:,3], alpha=0.5 )
ax3.plot(ir_plot_profile_3[:,2], IRheights, label="IR screen")
ax3.fill_betweenx(IRheights, ir_plot_profile_3[:,1], ir_plot_profile_3[:,3], alpha=0.5 )
ax3.set_xlabel(L"T~[\mathrm{^\circ C}]")
#ax3.set_ylabel(L"h~[\mathrm{m}]")
ax3.grid()
#ax3.set_xlim(0,13)
ax3.tick_params(axis="y", labelleft=false)
ax3.set_ylim(0,2)

(h, l) = ax3.get_legend_handles_labels()
ax_leg = fig.add_subplot(py"$(gs)[4,0:3]")
l1 = ax_leg.legend(h, l, loc="lower center")#,borderaxespad=0)
ax_leg.axis("off")
PyPlot.gca().add_artist(l1)
#PyPlot.savefig("/home/haugened/Documents/les_paper_24/figures/cmp_tprofiles_IR_LES.jpg")
=#
