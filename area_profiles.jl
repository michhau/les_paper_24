#=
xxxdocumentation
script for postprocessing openFOAM-LES-runs. Deal with the profiles from an aerial cross section
(e.g. plot) 
=#

using NCDatasets, DataFrames, StatsBase, LaTeXStrings, Dates, PyCall
import PyPlot
PyPlot.pygui(true)
lines = pyimport("matplotlib.lines")
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

#netCDF4-file containing the profile (length 3)
file = joinpath.(joinpath(casedir, "postProcessing/surfaces/yNormal/"), "yNormal_time8th.nc")
#location of .nc-file containing the forcing wind speed
uforcingfile = joinpath(casedir, "forcing.nc")

#define plottype
plottype = "points"
allowedplottypes = ["profiles", "points"]

#axis along which profile is taken ("x", "y", or "z")
profilealong = "z"

#value of x where bare to snow transition is
snowtransition = 8.0 #m

#profile location [m] in new fetch-distance coordinates for x (multiple profiles, -1 along profile axis)
profilelocations = [-0.2, 0.5, 2.0, 10.3]

#profile half width
profhalfwidth = 0.03

#quantity for profiles (T, u, v, w, wT, Tvar, tke, or z/L)
quantity_to_plot = "w"
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
listsonics = ["kaijo", "t2irg", "t2lcsat", "t2ucsat", "tjk"]
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

if !(quantity_to_plot in allowedquantities)
    println(string(quantity_to_plot, " is not in ", allowedquantities, ". Please set variable 'quantity_to_plot' to allowed value. Exiting."))
    exit()
else
    println(string("Plotting ", plottype, " of ", quantity_to_plot))
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
    vals = T
    x_label_text = L"T~\mathrm{[K]}"
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
            datatmp.T .+= 273.15
            @warn("Careful!! Rotate Kaijo data due to measurement position. Check pictures of setup!")
            new_w = datatmp.u
            new_u = datatmp.w
            datatmp.u = new_u
            datatmp.w = new_w
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, Second(50))
            flxtmp = func.avgflux(flxtmp, Second(50))
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
            datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, Second(50))
            flxtmp = func.avgflux(flxtmp, Second(50))
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
            datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, Second(50))
            flxtmp = func.avgflux(flxtmp, Second(50))
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
            datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, Second(50))
            flxtmp = func.avgflux(flxtmp, Second(50))
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
            datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, Second(50))
            flxtmp = func.avgflux(flxtmp, Second(50))
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            t2ucsatdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="t2ucsat", listsonics)] = htmp
        end
        if "tjk" in listsonics
            htmp = 5.0
            datatmp = func.readtotalturbasnetcdf(tjkinfile)
            datatmp.T .+= 273.15
            disallowmissing!(datatmp)
            func.drdf!(datatmp; blockdur=drtime)
            datatmp = turb.despiking(datatmp)
            datatmp = turb.interpolatemissing(datatmp)
            flxtmp = func.turbflux(datatmp, Second(50))
            flxtmp = func.avgflux(flxtmp, Second(50))
            flxtmp[:, "zoverL"] = htmp ./ flxtmp.L_highfreq
            datatmp = datatmp[sonicsstart .<= datatmp.time .<= sonicsend, :]
            flxtmp = flxtmp[sonicsstart .<= flxtmp.time .<= sonicsend, :]
            tjkdata = hcat(datatmp, flxtmp[:,2:end])
            sonic_heights[findfirst(x->x=="tjk", listsonics)] = htmp
        end
        sonics_read = true
    end

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

#################################################
if plottype == "profiles"
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
        profiles = profilesT
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

    #plot median and indicate profile locations
    medval = fill(NaN, size(vals, 1))
    for k in 1:size(vals, 1)
        medval[k] = median(filter(!isnan, vals[k, :]))
    end

    cmap = PyPlot.get_cmap("tab10")
    fig2 = PyPlot.figure(figsize=(20,4))
    ax1 = fig2.add_subplot(111)
    ax1.set_title(string("Median T - Indicating profile locations"))
    hm = ax1.tripcolor(x, z, medval, cmap="viridis", vmin=minimum(medval)*1.05, vmax=maximum(medval)*0.95)
    ax1.set_xlabel(L"x~[\mathrm{m}]")
    ax1.set_ylabel(L"h~[\mathrm{m}]")
    ax1.set_xlim(minimum(x), maximum(x))
    ax1.set_xlim(-0.5, 12)
    ax1.set_ylim(0,2)
    for k in 1:size(profilelocations, 1)
        l = lines.Line2D([profilelocations[k], profilelocations[k]], [meanheights[k], 5], color=cmap(k-1))
        ax1.fill_betweenx(collect(0:0.2:5), profilelocations[k]-profhalfwidth, profilelocations[k]+profhalfwidth, color=cmap(k-1), alpha=0.5 )
    end
    topo = ax1.plot(xtopo.x, xtopo.h, color="black", lw=3)
    ax1.fill_between(xtopo.x, zeros(length(xtopo.x)), xtopo.h, color="white")
    hm_cb = fig2.colorbar(hm, ax=ax1)
    hm_cb.set_label(L"T~[\mathrm{K}]")
    PyPlot.gca().set_aspect("equal")

    fig = PyPlot.figure()
    ax1 = fig.add_subplot(111)
    if postfilter
        ax1.set_title(filtercond)
    else
        ax1.set_title("No filtering")
    end
    for k in 1:size(profilelocations, 1)
        ax1.plot(medianprof[:,k,2], dist_along_profile[:, k], label=string("x=", profilelocations[k]-profhalfwidth, " - ", profilelocations[k]+profhalfwidth, "m"))
        ax1.fill_betweenx(dist_along_profile[:, k], medianprof[:,k,1], medianprof[:,k,3], alpha=0.3)
    end
    #ax1.plot(medianprof[:,2,2], dist_along_profile[:, 2], label=string("x=", profilelocations[2]-profhalfwidth, " - ", profilelocations[2]+profhalfwidth, "m"))
    #ax1.fill_betweenx(dist_along_profile[:, 2], medianprof[:,2,1], medianprof[:,2,3], alpha=0.3)
    #ax1.plot(medianprof[:,3,2], dist_along_profile[:, 3], label=string("x=", profilelocations[3]-profhalfwidth, " - ", profilelocations[3]+profhalfwidth, "m"))
    #ax1.fill_betweenx(dist_along_profile[:, 3], medianprof[:,3,1], medianprof[:,3,3], alpha=0.3)
    #calculate mean and IQR for sonics
    if quantity_symbol == :T
        #t2vent_excerpt[:, :vent_air_temp]
        ventmed = quantile(filter(!isnan, scaledTair.+273.15), [0.25,0.5,0.75])
        ax1.errorbar(ventmed[2], 1.0, xerr=[[abs(ventmed[1]-ventmed[2])], [abs(ventmed[3]-ventmed[2])]], fmt="o", label=L"T_{vent}~\mathrm{(scaled)}")
    end
    if loadsonics
        if "kaijo" in listsonics
            sonicmed = quantile(filter(!isnan, kaijodata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="kaijo", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="Kaijo (0.25-0.75)")
        end
        if "t1irg" in listsonics
            sonicmed = quantile(filter(!isnan, t1irgdata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t1irg", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="T1IRG (0.25-0.75)")
        end
        if "t2irg" in listsonics
            sonicmed = quantile(filter(!isnan, t2irgdata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t2irg", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="T2IRG (0.25-0.75)")
        end
        if "t2lcsat" in listsonics
            sonicmed = quantile(filter(!isnan, t2lcsatdata[:, sonic_symbol]), [0.25,0.5,0.75])
            ax1.errorbar(sonicmed[2], sonic_heights[findfirst(x->x=="t2lcsat", listsonics)], xerr=[[abs(sonicmed[1]-sonicmed[2])], [abs(sonicmed[3]-sonicmed[2])]], fmt="o", label="T2LCSAT (0.25-0.75)")
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
    ax1.set_xlabel(x_label_text)
    ax1.set_ylabel(L"h~\mathrm{[m]}")
    ax1.set_ylim(0,2.0)
    ax1.grid()
    ax1.legend()
    #######################################################

elseif plottype == "points"
    take1 = func.chooseval3D(p1, hcat(x, y, z))
    take2 = func.chooseval3D(p2, hcat(x, y, z))
    take3 = func.chooseval3D(p3, hcat(x, y, z))

    println(string("point1 (x, y, z) [m] = (", x[take1], ", ", y[take1], ", ", z[take1], ")"))
    println(string("point2 (x, y, z) [m] = (", x[take2], ", ", y[take2], ", ", z[take2], ")"))
    println(string("point3 (x, y, z) [m] = (", x[take3], ", ", y[take3], ", ", z[take3], ")"))

    take_point1 = vals[take1, skiptoidx:end]
    take_point2 = vals[take2, skiptoidx:end]
    take_point3 = vals[take3, skiptoidx:end]
    times_points = times[skiptoidx:end]

    if forcingdata_read_sucessfully
        #get height in forcing data corresponding to choice of z for profile
        if forcingheight == -1
            correspforcheightidx = findfirst(x->x==minimum(abs.(forcing_heights .- z[take1])), abs.(forcing_heights .- z[take1]))
            correspforcwind = forcing_u[correspforcheightidx, :, :]
            correspforcT = forcing_T[correspforcheightidx, :]
            correspforcheight = forcing_heights[correspforcheightidx]
        else
            correspforcheightidx = findfirst(x->x==minimum(abs.(forcing_heights .- forcingheight)), abs.(forcing_heights .- forcingheight))
            correspforcwind = forcing_u[correspforcheightidx, :, :]
            correspforcT = forcing_T[correspforcheightidx, :]
            correspforcheight = forcing_heights[correspforcheightidx]
        end

        #convert to DataFrame
        corrforce = DataFrame("time" => forcing_times, "U0" => correspforcwind[1,:], "U1" => correspforcwind[2,:], "U2" => correspforcwind[3,:], "T" => correspforcT)
    end

    #time series

    #plot wind time series with wind forcing time series
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("time [s]")
    ax.set_ylabel(x_label_text)
    if forcingdata_read_sucessfully .&& quantity_symbol in [:T, :U_0, :U_1, :U_2]
        ax.plot(forcing_times, corrforce[:,quantity_symbol], label=string("forcing at z=", round(correspforcheight, digits=2), "m"))
    end
    ax.plot(times_points, take_point1, label=string("x=", round(x[take1], digits=2), ", y=", round(y[take1], digits=2), ", z=", round(z[take1], digits=2)), alpha=0.5)
    ax.plot(times_points, take_point2, label=string("x=", round(x[take2], digits=2), ", y=", round(y[take2], digits=2), ", z=", round(z[take2], digits=2)), alpha=0.5)
    ax.plot(times_points, take_point3, label=string("x=", round(x[take3], digits=2), ", y=", round(y[take3], digits=2), ", z=", round(z[take3], digits=2)), alpha=0.5)
    ax.legend()
    ax.grid(true)

    #spectra

    #evaluate time steps
    timesteps = func.gettimesteps(times_points)

    #choose new time step
    timesteptotake = 0.064 #0.045 #[s]

    #plot histogram of time steps
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Time steps")
    ax.set_xlabel("time step [1e-3 s]")
    ax.set_ylabel("density")
    ax.hist(timesteps .* 1e3, bins=collect(0:0.05:100), density=true, label="distribution");
    ax.axvline(timesteptotake *1e3, 0, 1, color="red", label="time step taken")
    ax.grid(true)
    ax.legend()

    println(string("mean(timesteps) [s] = ", mean(timesteps), "; median(timesteps) [s] = ", median(timesteps)))

    #make new time series with equal distance between times
    finaltimes = func.createequaltimesteps(timesteptotake, times_points)
    equalprofile1 = func.linearmaptotimeseries(take_point1, times_points, finaltimes)
    equalprofile2 = func.linearmaptotimeseries(take_point2, times_points, finaltimes)
    equalprofile3 = func.linearmaptotimeseries(take_point3, times_points, finaltimes)

    if forcingdata_read_sucessfully
        forcetmp = copy(corrforce)
    end
    evaldftmp1 = copy(equalprofile1)
    evaldftmp2 = copy(equalprofile2)
    evaldftmp3 = copy(equalprofile3)

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

    #do it for forcing data separately
    if forcingdata_read_sucessfully .&& quantity_symbol in [:T, :U_0, :U_1, :U_2]
        ftdata = forcetmp[:, quantity_symbol] #evaldftmp.w .* evaldftmp.T
        dur = forcetmp.time[end] - forcetmp.time[1]
        rawvec = ftdata
        vec = ft.detrend(rawvec) 
        vect = ft.belltaper(vec)
        (FTraw, FTvec) = ft.fourierplus(vect)
        (Soff, freq) = ft.spectralenergydensity(FTvec, dur)
        freqtimesSoff = freq.*Soff
        (freq, freqtimesSoff) = ft.logavg(freq, freqtimesSoff, 0.05)
    end

    ftdata1 = evaldftmp1 #evaldftmp.w .* evaldftmp.T
    ftdata2 = evaldftmp2 #evaldftmp.w .* evaldftmp.T
    ftdata3 = evaldftmp3 #evaldftmp.w .* evaldftmp.T

    #duration in seconds
    dur1 = finaltimes[end] - finaltimes[1]
    dur2 = finaltimes[end] - finaltimes[1]
    dur3 = finaltimes[end] - finaltimes[1]

    println("Preprocessing...")
    #select vector from dataset
    rawvec1 = ftdata1
    rawvec2 = ftdata2
    rawvec3 = ftdata3
    #detrend the vector (subtract linear fit)
    vec1 = ft.detrend(rawvec1)
    vec2 = ft.detrend(rawvec2)
    vec3 = ft.detrend(rawvec3)
    #apply bell taper
    vect1 = ft.belltaper(vec1)
    vect2 = ft.belltaper(vec2)
    vect3 = ft.belltaper(vec3)

    println("Applying Fourier-Trafo...")
    (FTraw1, FTvec1) = ft.fourierplus(vect1)
    (Soff1, freq1) = ft.spectralenergydensity(FTvec1, dur1)
    (FTraw2, FTvec2) = ft.fourierplus(vect2)
    (Soff2, freq2) = ft.spectralenergydensity(FTvec2, dur2)
    (FTraw3, FTvec3) = ft.fourierplus(vect3)
    (Soff3, freq3) = ft.spectralenergydensity(FTvec3, dur3)

    freqtimesSoff1 = freq1.*Soff1
    (freq1, freqtimesSoff1) = ft.logavg(freq1, freqtimesSoff1, 0.05)
    freqtimesSoff2 = freq2.*Soff2
    (freq2, freqtimesSoff2) = ft.logavg(freq2, freqtimesSoff2, 0.05)
    freqtimesSoff3 = freq3.*Soff3
    (freq3, freqtimesSoff3) = ft.logavg(freq3, freqtimesSoff3, 0.05)

    PyPlot.pygui(true)
    fig2 = PyPlot.figure()
    ax = PyPlot.gca()
    ax.set_title(string("Fourier-Decomposition ", quantity_to_plot))
    ax.set_xlabel("f [Hz]")
    ax.set_ylabel("fS(f)")
    #=if forcingdata_read_sucessfully
        ax.plot(freq, func.movingaverage(freqtimesSoff, 8), label=string("forcing, z=", round(correspforcheight, digits=2), "m"))
    end=#
    tmp1 = func.movingaverage(freqtimesSoff1, 8)
    tmp2 = func.movingaverage(freqtimesSoff2, 8)
    tmp3 = func.movingaverage(freqtimesSoff3, 8)
    ax.plot(freq1, tmp1, label=string("x=", round(x[take1], digits=2), ", y=", round(y[take1], digits=2), ", z=", round(z[take1], digits=2), "m"))
    ax.plot(freq2, tmp2, label=string("x=", round(x[take2], digits=2), ", y=", round(y[take2], digits=2), ", z=", round(z[take2], digits=2), "m"))
    ax.plot(freq3, tmp3, label=string("x=", round(x[take3], digits=2), ", y=", round(y[take3], digits=2), ", z=", round(z[take3], digits=2), "m"))
    #do it for sonics separately
    if loadsonics
        if "kaijo" in listsonics
            datatmp = copy(kaijodata[:, [:time, sonic_symbol]])
            dur_sonic = Dates.value(Second(datatmp.time[end] - datatmp.time[1]))
            rawvec_sonic = datatmp[:, sonic_symbol]
            vec_sonic = ft.detrend(rawvec_sonic) 
            vect_sonic = ft.belltaper(vec_sonic)
            (FTraw_sonic, FTvec_sonic) = ft.fourierplus(vect_sonic)
            (Soff_sonic, freq_sonic) = ft.spectralenergydensity(FTvec_sonic, dur_sonic)
            freqtimesSoff_sonic = freq_sonic .* Soff_sonic
            (freq_sonic, freqtimesSoff_sonic) = ft.logavg(freq_sonic, freqtimesSoff_sonic, 0.05)
            ax.plot(freq_sonic, func.movingaverage(freqtimesSoff_sonic, 8), label=string("Kaijo @ ", sonic_heights[findfirst(x->x=="kaijo", listsonics)], "m"))
        end
        if "t1irg" in listsonics
            datatmp = copy(t1irgdata[:, [:time, sonic_symbol]])
            dur_sonic = Dates.value(Second(datatmp.time[end] - datatmp.time[1]))
            rawvec_sonic = datatmp[:, sonic_symbol]
            vec_sonic = ft.detrend(rawvec_sonic) 
            vect_sonic = ft.belltaper(vec_sonic)
            (FTraw_sonic, FTvec_sonic) = ft.fourierplus(vect_sonic)
            (Soff_sonic, freq_sonic) = ft.spectralenergydensity(FTvec_sonic, dur_sonic)
            freqtimesSoff_sonic = freq_sonic .* Soff_sonic
            (freq_sonic, freqtimesSoff_sonic) = ft.logavg(freq_sonic, freqtimesSoff_sonic, 0.05)
            ax.plot(freq_sonic, func.movingaverage(freqtimesSoff_sonic, 8), label=string("T1IRG @ ", sonic_heights[findfirst(x->x=="t1irg", listsonics)], "m"))
        end
        if "t2irg" in listsonics
            datatmp = copy(t2irgdata[:, [:time, sonic_symbol]])
            dur_sonic = Dates.value(Second(datatmp.time[end] - datatmp.time[1]))
            rawvec_sonic = datatmp[:, sonic_symbol]
            vec_sonic = ft.detrend(rawvec_sonic) 
            vect_sonic = ft.belltaper(vec_sonic)
            (FTraw_sonic, FTvec_sonic) = ft.fourierplus(vect_sonic)
            (Soff_sonic, freq_sonic) = ft.spectralenergydensity(FTvec_sonic, dur_sonic)
            freqtimesSoff_sonic = freq_sonic .* Soff_sonic
            (freq_sonic, freqtimesSoff_sonic) = ft.logavg(freq_sonic, freqtimesSoff_sonic, 0.05)
            ax.plot(freq_sonic, func.movingaverage(freqtimesSoff_sonic, 8), label=string("T2IRG @ ", sonic_heights[findfirst(x->x=="t2irg", listsonics)], "m"))
        end
        if "t2lcsat" in listsonics
            datatmp = copy(t2lcsatdata[:, [:time, sonic_symbol]])
            dur_sonic = Dates.value(Second(datatmp.time[end] - datatmp.time[1]))
            rawvec_sonic = datatmp[:, sonic_symbol]
            vec_sonic = ft.detrend(rawvec_sonic) 
            vect_sonic = ft.belltaper(vec_sonic)
            (FTraw_sonic, FTvec_sonic) = ft.fourierplus(vect_sonic)
            (Soff_sonic, freq_sonic) = ft.spectralenergydensity(FTvec_sonic, dur_sonic)
            freqtimesSoff_sonic = freq_sonic .* Soff_sonic
            (freq_sonic, freqtimesSoff_sonic) = ft.logavg(freq_sonic, freqtimesSoff_sonic, 0.05)
            ax.plot(freq_sonic, func.movingaverage(freqtimesSoff_sonic, 8), label=string("T2LCSAT @ ", sonic_heights[findfirst(x->x=="t2lcsat", listsonics)], "m"))
        end
        if "t2ucsat" in listsonics
            datatmp = copy(t2ucsatdata[:, [:time, sonic_symbol]])
            dur_sonic = Dates.value(Second(datatmp.time[end] - datatmp.time[1]))
            rawvec_sonic = datatmp[:, sonic_symbol]
            vec_sonic = ft.detrend(rawvec_sonic) 
            vect_sonic = ft.belltaper(vec_sonic)
            (FTraw_sonic, FTvec_sonic) = ft.fourierplus(vect_sonic)
            (Soff_sonic, freq_sonic) = ft.spectralenergydensity(FTvec_sonic, dur_sonic)
            freqtimesSoff_sonic = freq_sonic .* Soff_sonic
            (freq_sonic, freqtimesSoff_sonic) = ft.logavg(freq_sonic, freqtimesSoff_sonic, 0.05)
            ax.plot(freq_sonic, func.movingaverage(freqtimesSoff_sonic, 8), label=string("T2UCSAT @ ", sonic_heights[findfirst(x->x=="t2ucsat", listsonics)], "m"))
        end
        if "tjk" in listsonics
            datatmp = copy(tjkdata[:, [:time, sonic_symbol]])
            dur_sonic = Dates.value(Second(datatmp.time[end] - datatmp.time[1]))
            rawvec_sonic = datatmp[:, sonic_symbol]
            vec_sonic = ft.detrend(rawvec_sonic) 
            vect_sonic = ft.belltaper(vec_sonic)
            (FTraw_sonic, FTvec_sonic) = ft.fourierplus(vect_sonic)
            (Soff_sonic, freq_sonic) = ft.spectralenergydensity(FTvec_sonic, dur_sonic)
            freqtimesSoff_sonic = freq_sonic .* Soff_sonic
            (freq_sonic, freqtimesSoff_sonic) = ft.logavg(freq_sonic, freqtimesSoff_sonic, 0.05)
            ax.plot(freq_sonic, func.movingaverage(freqtimesSoff_sonic, 8), label=string("TJK @ ", sonic_heights[findfirst(x->x=="tjk", listsonics)], "m"))
        end
    end
    xtheory=collect(2:0.05:8)
    ytheory=exp.(-2*log.(xtheory)/3)./15
    ax.plot(xtheory,ytheory, color="black", label="e^{-2/3}")
    PyPlot.xscale("log")
    PyPlot.yscale("log")
    #ax.set_xlim([0.01, 10])
    ax.grid()
    ax.legend()
end