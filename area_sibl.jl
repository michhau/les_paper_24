#=
xxxdocumentation
script for postprocessing openFOAM-LES-runs. Evaluate SIBL growth
=#

using NCDatasets, DataFrames, StatsBase, LaTeXStrings, PyCall, Dates, LsqFit, CSV
import PyPlot
PyPlot.pygui(true)
GridSpec = pyimport("matplotlib.gridspec")
#mpwidgets = pyimport("matplotlib.widgets")
animation = pyimport("matplotlib.animation")
numpy = pyimport("numpy")
cramericm = pyimport("cmcrameri.cm")

#outer directory of case
casedir = "/home/haugened/Documents/openfoam/old_runs/240211_newnewMesh/"
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
file = joinpath.(readfromfolder, "yNormal.nc")

#value of x where bare to snow transition is
snowtransition = 6.0 #m

reyavgtime = Second(50)

#seconds to skip (spinuptime)
skipsec = 30

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
anistring = string(joinpath(casedir, "video.mp4"))
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
    plotlim = 0.01
    cblabel = L"\overline{w'T'}~\mathrm{[K~m~s^{-1}]}"
    if postfilter
        figtitle = string("Median w'T' skipping first ", skipsec, "s; ", filtercond)
    else
        figtitle = string("Median w'T' skipping first ", skipsec, "s; No filtering")
    end

    #get stuff for evaluating wT
    vals_w = U[:,:,3]
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
    postfiltercond = fill(true, size(quant_crosssection, 3))
    println("No post-filtering applied.")
end

quant_crosssection[:,:, .!(postfiltercond)] .= NaN

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
        if sibldepth_idx[iloc] == 0
            idcs = sibllocs[iloc]-calchalfwidth .<= xtopo.x .<= sibllocs[iloc]+calchalfwidth
            if count(idcs )== 0
                tmpmini = abs.(xtopo.x .- minimum(sibllocs[iloc]))
                mintmpmini = minimum(tmpmini)
                sibldepth[iloc] = xtopo.h[findfirst(x->x==mintmpmini, tmpmini)]
            else
                sibldepth[iloc] = mean(xtopo.h[idcs])
            end
        else
            sibldepth[iloc] = dist_along_profile[sibldepth_idx[iloc]]
        end
    end

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

    #plot SIBL depth
    fig = PyPlot.figure()
    ax2 = fig.add_subplot(111)
    ax2.set_title(string("SIBL depth from ", diagnosefromquantity," & Brutsaert82 fit (h=c*x^b); c=", round(brutsaert82_param[1]*100, digits=2), "e-2, b=", round(brutsaert82_param[2], digits=2)))
    ax2.plot(sibllocs, sibldepth, label="LES")
    ax2.plot(sibllocs[xbiggerequal0], brutsaert82(sibllocs[xbiggerequal0], brutsaert82_param), label="Brutsaert82 fit")
    ax2.set_xlim(-0, 4)
    ax2.set_ylim(0, 0.35)
    ax2.set_xlabel(L"x~\mathrm{[m]}")
    ax2.set_ylabel("SIBL depth [m]")
    ax2.grid()
    ax2.legend()
end

#plot quant_crosssection
fig = PyPlot.figure(figsize=(30,4))
ax1 = fig.add_subplot(111)
ax1.set_title(figtitle)
hm = ax1.tripcolor(x, z, quant_median, cmap=cramericm.vik, vmin=-plotlim, vmax=plotlim, shading="gouraud")
topo = ax1.plot(xtopo.x, xtopo.h, color="black", lw=3)
ax1.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
ax1.set_xlabel(L"x~[\mathrm{m}]")
ax1.set_ylabel(L"z~[\mathrm{m}]")
if calcsibl
    ax1.plot(sibllocs[xbiggerequal0], sibldepth[xbiggerequal0], lw=3, c="green")
end
ax1.set_ylim(0,2)
ax1.set_xlim(-6, 4)
hm_cb = fig.colorbar(hm, ax=ax1)
hm_cb.set_label(cblabel)
PyPlot.gca().set_aspect("equal")

#=
"Update-function for the animation"
function animupdate(frame::Int64)
    fig.suptitle(string("LES t = ", round(times[frame], digits=1), " s"))
    newdata = transpose(vals[:,2:end,frame])
    newdata = newdata[1:end-1, 1:end-1]
    hm.set_array(numpy.ravel(newdata))
    #if mod(frame,200) == 0
    #    println(string(frame, "/", length(animidcs)))
    #end
    return hm
end

fig = PyPlot.figure(figsize=(23,7))
ax1 = fig.add_subplot(111)
hm = ax1.pcolormesh(x_anim, z[2:end], transpose(vals[:,2:end,100]), vmin=valmin, vmax=valmax, shading="flat")
ax1.set_xlabel(L"x~[\mathrm{m}]")
ax1.set_ylabel(L"z~[\mathrm{m}]")
hm_cb = fig.colorbar(hm, ax=ax1)
hm_cb.set_label(label_text)
PyPlot.gca().set_aspect("equal")
ani = animation.FuncAnimation(fig, animupdate, frames=animidcs[1:end], interval=intervalmsec, repeat=false)#, blit=false)
ani.save(anistring, fps=framespersec)
=#