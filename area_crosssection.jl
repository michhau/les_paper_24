#=
xxxdocumentation
script for postprocessing openFOAM-LES-runs. Create synthetic screen videos from crosssection
(e.g. plot) 
=#

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
plotstyle = "animation"
allowedplotstyle = ["median", "animation"]

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
anisave = true

#file name of video (only for plotstyle == "animation")
anistring = string(joinpath(casedir, "video_T.mp4"))

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

    fig = PyPlot.figure(figsize=(30,4))
    ax1 = fig.add_subplot(111)
    ax1.set_title(string("Median ", figtitle, ", skipping ", skipsec, "s (spin up)"))
    hm = ax1.tripcolor(x, z, medianvals, cmap=colmap, vmin=valmin, vmax=valmax)
    topo = ax1.plot(xtopo.x, xtopo.h, color="black", lw=3)
    ax1.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
    ax1.set_xlabel(L"x~[\mathrm{m}]")
    ax1.set_ylabel(L"z~[\mathrm{m}]")
    ax1.set_xlim(minimum(x), maximum(x))
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

#=
#sorting for concentation of different point files
x1 = copy(x)
y1 = copy(y)
z1 = copy(z)
times1 = copy(times)
T1 = copy(T)
U1 = copy(U)

a1 = DataFrame(x = x1, y = y1, z=z1)
a1sortperm = sortperm(a1)
a1sorted = a1[a1sortperm,:]

T1sorted = fill(NaN, size(T1, 1), size(T1, 2))
U1sorted = fill(NaN, size(U1, 1), size(U1, 2), size(U1, 3))
for i in 1:size(T, 2)
    T1sorted[:, i] = T1[a1sortperm, i]
    U1sorted[:, i, 1] = U1[a1sortperm, i, 1]
    U1sorted[:, i, 2] = U1[a1sortperm, i, 2]
    U1sorted[:, i, 3] = U1[a1sortperm, i, 3]
end

x1 = copy(a1sorted[:, 1])
y1 = copy(a1sorted[:, 2])
z1 = copy(a1sorted[:, 3])
T1 = copy(T1sorted)
U1 = copy(U1sorted)

=#