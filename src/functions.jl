######################################################
###               FUNCTION MODULE                  ###
###            author: Michi Haugeneder            ###
######################################################
#=
Module containing functions for creating and evaluating
openFOAM LES simulations
=#

module func

using DataFrames, NCDatasets, Dates, StatsBase, Statistics, LsqFit, CSV

export avgflux, chooseheight, chooselistval, chooseval,
chooseval3D, createbottom, createequaltimesteps,
createtimeidcs, createxyzidcs, csvtodataframe, 
diagnosesibldepth, drdf!, filterfrontpoints, findmultipleentries,
getanimidcs, getdomainlength, gettimesteps,
linearmaptotimeseries, logscaletemperature,
logscalewindspeed, movingaverage, readaerialprofilefromnetcdf, readbinarypoints, 
readforcingfromnetcdf, readpoints, readprofilefromnetcdf,
readT, readtotalturbasnetcdf, readU, reducedoubleentries, reducenonnans,
saveasnetcdf, savelinasnetcdf,
saveforcingasnetcdf, timedepvectoroutstring, turbflux, vectoroutstring

"""
    avgflux(data::DataFrame, peri::Period)

Average the turbulent fluxes in the DataFrame
"""
function avgflux(data::DataFrame, peri::Period)
    ele = round(Int, Millisecond(peri) / Millisecond(50))
    fluxavg = similar(data)
    fluxavg.time = data.time
    for iname in names(data)[2:end]
        fluxavg[:, iname] = movingaverage(data[:, iname], ele)
    end
    return fluxavg
end

"""
    chooseheight(desiredvalue::Number, heights::Vector)::Int64

Return index with closest match to desired height.
"""
function chooseheight(desiredvalue::Number, heights::Vector)::Int64
    minvector = abs.(heights .- desiredvalue)
    mini = minimum(minvector)
    miniidx = findfirst(x->x==mini, minvector)
    println(string("User input height: ", desiredvalue, "; height taken Z = ", heights[miniidx]))
    return miniidx
end

"""
    chooselistval(profilealong::Char, desiredminvalue::Number, desiredmaxvalue::Number, points::Vector)::Vector

Return row indices in pts array matching desired value within range.
"""
function chooselistval(profilealong::String, desiredminvalue::Number, desiredmaxvalue::Number, points::Array)::Vector
    if !(profilealong in ["x", "z"])
        @error("Wrong profilealong variable. Must be x or z")
        return [false]
    elseif profilealong == "x"
        idcs = desiredminvalue .<= points[:, 3]  .<= desiredmaxvalue
    elseif profilealong == "z"
        idcs = desiredminvalue .<= points[:, 1]  .<= desiredmaxvalue
    end
    return idcs
end

"""
    chooseval(desiredvalue::Number, vals::Vector)::Int64

Return index with closest match to desired height.
"""
function chooseval(desiredvalue::Number, vals::Vector)::Int64
    minvector = abs.(vals .- desiredvalue)
    mini = minimum(minvector)
    miniidx = findfirst(x->x==mini, minvector)
    #println(string("User input height: ", desiredvalue, "; height taken Z = ", vals[miniidx]))
    return miniidx
end

"""
    chooseval3D(desiredvalue::Vector, vals::Array)::Int64

Return indices with closest match to desired height.
"""
function chooseval3D(desiredvalue::Vector, vals::Array)::Int64
    dist_tmp = fill(NaN, size(vals, 1), size(vals, 2))
    for i in 1:size(vals, 2)
        dist_tmp[:, i] = abs.(vals[:, i] .- desiredvalue[i])
    end
    xyz_dist = fill(NaN, size(vals, 1))
    for i in 1:size(vals, 1)
        xyz_dist[i] = sum(dist_tmp[i, :])
    end
    mini = minimum(xyz_dist)
    return findfirst(x->x==mini, xyz_dist)
end

"""
    createbottom(xin, zin, interval=0.05)::DataFrame

Find the lowest grid points and create a floor line for the plots
"""
function createbottom(xin, zin, interval=0.05)::DataFrame
    minx = minimum(xin)
    maxx = maximum(xin)
    xvec = collect(minx:interval:maxx)
    xround = round.(xin, digits=2)
    xout = zeros(Float64, 0)
    hout = zeros(Float64, 0)

    for i in 1:length(xvec)
        if i==1
            idcs = findall(x->x<=xvec[i], xround)
        elseif i== length(xvec)
            idcs = findall(x->x>xvec[i], xround)
        else
            idcs = findall(x->xvec[i]-((xvec[i]-xvec[i-1])/2) < x <= ((xvec[i+1]-xvec[i])/2) + xvec[i], xround)
        end
        if length(idcs)>0
            push!(hout, minimum(zin[idcs]))
            push!(xout, xvec[i])
        end
    end
    out = DataFrame(x=xout, h=hout)
    return out
end

"""
    createequaltimesteps(new_timestep::Number, old_times::Vector)::Vector

Calculate vector with new time steps equally spaced with new_timestep spacing from old_times
"""
function createequaltimesteps(new_timestep::Number, old_times::Vector)::Vector
    mintimefactor = old_times[1] / new_timestep
    maxtimefactor = old_times[end] / new_timestep
    return collect(ceil(mintimefactor):1:floor(maxtimefactor)).*new_timestep
end

"""
    createtimeidcs(inittimefolders2::Vector)

Return indices for time folders to assign in array
"""
function createtimeidcs(inittimefolders2::Vector)
    timesparsed = parse.(Float64, inittimefolders2)
    timessorted = sort(timesparsed)
    time_idx = zeros(Int64, size(timesparsed, 1))

    for i in 1:length(inittimefolders2)
        time_idx[i] = findfirst(x->x==timesparsed[i], timessorted)
    end

    return time_idx
end

"""
    createxyzidcs(pts::Array)

Create idx vectors (x,y, and z) for points to assign into array
"""
function createxyzidcs(pts::Array)
    unique_x = sort!(unique(pts[:,1]))
    unique_y = sort!(unique(pts[:,2]))
    unique_z = sort!(unique(pts[:,3]))


    x_idx = zeros(Int64, size(pts, 1))
    y_idx = zeros(Int64, size(pts, 1))
    z_idx = zeros(Int64, size(pts, 1))

    #create indices x,y,z of every point in pts
    for i in 1:size(pts, 1)
        x_idx[i] = findfirst(x->x==pts[i,1], unique_x)
        y_idx[i] = findfirst(x->x==pts[i,2], unique_y)
        z_idx[i] = findfirst(x->x==pts[i,3], unique_z)
    end

    return x_idx, y_idx, z_idx, unique_x, unique_y, unique_z
end

"""
    csvtodataframe(source::String)

Read data from .csv-file and convert it to a DataFrame that is returned.
Usual way of importing data after .csv-file is being created once.
"""
function csvtodataframe(source::String)
    @info("Reading data into DataFrame", source)
    df = DataFrame(CSV.File(source))
    allowmissing!(df)
    for col in eachcol(df)
        replace!(col, NaN => missing)
    end
    return df
end

"""
    detectnanandgap(data::DataFrame, gapthresh::Period)

Detect also nan-filled gaps addiitonally to above function 'detectgaps'
"""
function detectnanandgap(data::DataFrame, gapthresh::Period)
    gaps = DataFrame(idx_before_gap=Int64[], idx_after_gap=Int64[])
    for i in 1:size(data, 1)-1
        if data.time[i+1] - data.time[i] > gapthresh
            push!(gaps, [i, i+1])
        end
    end

    idxthresh = gapthresh/Millisecond(50)
    nanidcs = findall(x->(ismissing(x) | isnan(x)), data.u)

    lengthcurrser = 1
    serstartidx = 1

    for i in 2:length(nanidcs)
        if nanidcs[i] == nanidcs[i-1] +1
            lengthcurrser += 1
        else
            if lengthcurrser >= idxthresh
                push!(gaps, [serstartidx-1, nanidcs[i]])
            end
            lengthcurrser = 1
            serstartidx = nanidcs[i]
        end
    end

    sort(gaps, :idx_before_gap)
end

"""
    diagnosesibldepth(verticalwTprofile::Vector)

Diagnose depth of SIBL using either w'T' or z/L. Return index
"""
function diagnosesibldepth(verticalwTprofile::Vector, diagquantity::String)
    if diagquantity == "wT"
        try
            sibldepth_movavg = findfirst(x->x>0, movingaverage(verticalwTprofile, 10)[3:end])+2
        catch e
            sibldepth_movavg = 0
        end
        try
            sibldepth_normal = findfirst(x->x>0, verticalwTprofile[3:end])+2
        catch e
            sibldepth_normal = 0
        end
        sibldepth = maximum([sibldepth_movavg, sibldepth_normal])
        if sibldepth <= 3
            sibldepth = 0
        end
    elseif diagquantity == "z/L"
        sibldepth_movavg = findfirst(x->x<0, movingaverage(verticalwTprofile, 10)[3:end])+2
        sibldepth_normal = findfirst(x->x<0, verticalwTprofile[3:end])+2
        sibldepth = maximum([sibldepth_movavg, sibldepth_normal])
        if sibldepth <= 3
            sibldepth = 0
        end
    end
    return sibldepth
end

"""
    drdf!(data::DataFrame; blockdur=Minute(30), periodwise=true, gapthresh=Minute(10))

Perform in-place double rotation on DataFrame
(using col. names 'u','v','w') to set mean(v)=mean(w)=0., blockdur=duration of blocks, periodwise refers to gaps.
"""
function drdf!(data::DataFrame; blockdur=Minute(30), periodwise=true, gapthresh=Minute(10))
    endidcs = zeros(Int64, 0)
    startidcs = zeros(Int64, 0)
    currnanidx = 1
    println("Double rotation for blocks of ", blockdur)
    blockduridx = round(Int, blockdur/Millisecond(50))
    startidx = 1
    if periodwise
        gaps = detectnanandgap(data, gapthresh)
        nanendidcs = gaps.idx_before_gap
        while startidx < size(data, 1)
            push!(startidcs, startidx)
            if currnanidx<=size(nanendidcs, 1) && nanendidcs[currnanidx] <= startidx+blockduridx-1
                push!(endidcs, nanendidcs[currnanidx])
                startidx = gaps.idx_after_gap[currnanidx]
                currnanidx += 1
            else
                push!(endidcs, startidx+blockduridx-1)
                startidx+=blockduridx-1
            end
        end
        endidcs[end] = size(data, 1)
        @info("Double rotation period-wise. Considering data gaps (e.g. due to reposition).")
        println(size(gaps, 1) + 1, " periods")
    else
        @info("Performing in-place Double Rotation without considering data gaps (DR over gaps as well).")
        while startidx < size(data, 1)
            push!(startidcs, startidx)
            push!(endidcs, startidx+blockduridx-1)
            startidx+=blockduridx-1
        end
        endidcs[end] = size(data, 1)

    end

    for per in axes(endidcs, 1)
        startidx = startidcs[per]
        endidx = endidcs[per]

        datatouse = data[startidx:endidx, :]

        #creating necessary temporary arrays
        data1 = Array{Union{Missing,Float64}}(missing, size(datatouse, 1), 2)
        data2 = Array{Union{Missing,Float64}}(missing, size(datatouse, 1), 3)

        #calculating averages
        meanu = mean(filter(!isnan, skipmissing(datatouse[:, :u])))
        meanv = mean(filter(!isnan, skipmissing(datatouse[:, :v])))
        meanw = mean(filter(!isnan, skipmissing(datatouse[:, :w])))

        #=
        #calculate mean windspeed and -direction
        mean_wndspd = sqrt(meanu^2 + meanv^2 + meanw^2)
        mean_dir = mod(rad2deg(atan(-meanv, meanu)), 360)
        @info("Info:", mean_wndspd, mean_dir)
        =#

        #calculating first rotation angle alpha[rad]
        alpha = atan(meanv, meanu)
        #@show rad2deg(alpha)
        #rotating the windfield to obtain meanv=0 (w-component stays the same)
        data1[:, 1] = datatouse[:, :u] .* cos(alpha) + datatouse[:, :v] .* sin(alpha)
        data1[:, 2] = -datatouse[:, :u] .* sin(alpha) + datatouse[:, :v] .* cos(alpha)
        #calculating second rotation angle beta[rad]
        beta = atan(meanw, mean(filter(!isnan, skipmissing(data1[:, 1]))))
        #@show rad2deg(beta)
        #rotating the windfield to obtain meanw=0 (v-component stays the same)
        data2[:, 1] = data1[:, 1] .* cos(beta) + datatouse[:, :w] .* sin(beta)
        data2[:, 2] = data1[:, 2]
        data2[:, 3] = -data1[:, 1] .* sin(beta) + datatouse[:, :w] .* cos(beta)

        #overwrite the input-data
        data[startidx:endidx, :u] = data2[:, 1]
        data[startidx:endidx, :v] = data2[:, 2]
        data[startidx:endidx, :w] = data2[:, 3]
    end
end

"""
    filterfrontpoints(points::Array)::Array

Take in all points and filter out just inlet (x==0)
"""
function filterfrontpoints(points::Array)::Array
    frontidcs = findall(x->x<1e-7, points[:,1])
    return points[frontidcs, :]
end

"""
    findmultipleentries(x_idx::Vector, y_idx::Vector, z_idx::Vector)

Find entries in arrays that get assigned to multiple times
"""
function findmultipleentries(x_idx::Vector, y_idx::Vector, z_idx::Vector)
    tmparray = fill(0, length(unique(x_idx)), length(unique(y_idx)), length(unique(z_idx)))
    multiplejs = zeros(Int64, 0)
    multiple_idcs = zeros(Int64, 0, 3)
    for j in 1:length(x_idx)
        tmparray[x_idx[j], y_idx[j], z_idx[j]] += 1
        if tmparray[x_idx[j], y_idx[j], z_idx[j]] > 1
            multiplejs = vcat(multiplejs, [j])
            multiple_idcs = vcat(multiple_idcs, [x_idx[j] y_idx[j] z_idx[j]])
        end
    end

    j_to_idx = 1
    new_multiple_idcs = zeros(Int64, 1, 3) 
    new_multiple_idcs[1, :] = multiple_idcs[1,:]
    for k in 2:size(multiple_idcs, 1)
        already_found = 0
        for l in 1:size(new_multiple_idcs, 1)
            if multiple_idcs[k, :] == new_multiple_idcs[l, :]
                already_found = l
                break
            end
        end
        if already_found != 0
            j_to_idx = vcat(j_to_idx, already_found)
        else
            new_multiple_idcs = cat(new_multiple_idcs, transpose(multiple_idcs[k,:]), dims = 1)
            j_to_idx = vcat(j_to_idx, size(new_multiple_idcs, 1))
        end
    end

    return multiplejs, j_to_idx, new_multiple_idcs
end

"""
    getanimidcs(tr::Vector, intv::Int)::Vector

Get indices for animation. intv is interval [ms] for animation given in plotting command
"""
function getanimidcs(tr::Vector, intv::Int)::Vector
    animidcs = zeros(Int, size(collect(1:intv:tr[end]*1000), 1))
    currminidx = 1
    animidcs[1] = currminidx
    for i in 2:size(animidcs, 1)
        if currminidx == length(tr)
            animidcs[i:end] .= currminidx
            break
        elseif abs(tr[currminidx] - (i-1)*intv/1000) > abs(tr[currminidx+1] - (i-1)*intv/1000)
            newmin = abs(tr[currminidx+1] - (i-1)*intv/1000)
            k=1
            if currminidx+k >= length(tr)
                newmin = abs(tr[end]-(i-1)*intv/1000)
            else
                while (abs(tr[currminidx+k+1] - (i-1)*intv/1000) < newmin)
                    newmin = abs(tr[currminidx+k+1]-(i-1)*intv/1000)
                    k += 1
                    if currminidx+k+1 > length(tr)
                        break;
                    end
                end
            end
            currminidx += k
        end
        animidcs[i] = currminidx
    end
    return animidcs
end

"""
    getdomainlength(points::Array)::Number

Get the length of the domain from the points-file.
"""
function getdomainlength(points::Array)::Number
    return maximum(points[:,1])
end

"""
    gettimesteps(profiledata::Vector)::Vector

Get the time steps of the sampled profile
"""
function gettimesteps(profiledata::Vector)::Vector
    timesteps = zeros(Float64, size(profiledata, 1)-1)
    for i in 1:length(timesteps)
        timesteps[i] = profiledata[i+1] - profiledata[i]
    end
    return timesteps
end

"""
    linearmaptotimeseries(datain::Vector, timesin::Vector, timesout::Vector)::Vector

Create new synthetic time series with equal time steps by linear interpolation.
"""
function linearmaptotimeseries(datain::Vector, timesin::Vector, timesout::Vector)::Vector
    dataout = zeros(eltype(datain), length(timesout))
    for jrow in 1:length(timesout)
        timediffvec = timesin .- timesout[jrow]
        isexact = findfirst(x->x==0, timediffvec)
        if !isnothing(isexact) #value is contained in initial data set
            dataout[jrow] = datain[isexact]
        else #linear interpolation
            prevtimeidx = findlast(x->x<0, timediffvec)
            prevtime = timesin[prevtimeidx]
            nexttimeidx = findfirst(x->x>0, timediffvec)
            nexttime = timesin[nexttimeidx]
            prevvalue = datain[prevtimeidx]
            nextvalue = datain[nexttimeidx]

            #linear interpolation y = m*x + t
            tmptime = timesout[jrow] .- prevtime
            tmpnexttime = nexttime .- prevtime

            m = (nextvalue - prevvalue) / tmpnexttime
            dataout[jrow] = m * tmptime + prevvalue
        end
    end
    return dataout
end

"""
    logscaletemperature(t_zref::Number, zref::Number, t_srf::Number, z0::Number, zend)

Scale air temperature from measurement T(zref) to T(zend) with T(0) given by t_srf.
"""
function logscaletemperature(t_zref::Number, zref::Number, t_srf::Number, z0::Number, zend)
    #equation: T(z) = a*log(z/z0)+t_srf
    a = (t_zref-t_srf)/log(zref/z0) #t_srf actually at z=z0, not z=0
    logfactor = log.(zend ./ z0)
    logfactor[logfactor .< 0] .= 0.0
    return a .* logfactor .+ t_srf
end

"""
    scalewindspeed(u_zref::Number, zref::Number, z0::Number, zend)

Scale wind speed from measurement u(zref) to u(zend)
"""
function logscalewindspeed(u_zref::Number, zref::Number, z0::Number, zend)
    u_fric_over_kappa = log(zref/z0)^(-1)*u_zref
    logfactor = log.(zend ./ z0)
    logfactor[logfactor .< 0] .= 0.0
    return u_fric_over_kappa.*logfactor
end

"""
    movingaverage(X::Vector, numofele::Integer)

Create moving average (michi 16.09.2021); omitting NaNs!!
"""
function movingaverage(X::Vector, numofele::Integer)
    BackDelta = div(numofele, 2)
    ForwardDelta = isodd(numofele) ? div(numofele, 2) : div(numofele, 2) - 1
    len = size(X, 1)
    #create vector with vec_isnan[i]=1 if X[i]=NaN, 0 otherwise
    vec_isnan = isnan.(X)
    if numofele >= len
        #println("#avg elements >= size(vector). Returning mean.")
        return ones(Float64, len) .* mean(filter(!isnan, X))
    else
        Y = ones(Float64, size(X, 1))
        firstnonnan = findfirst(x -> x == false, vec_isnan)
        len = findlast(x -> x == false, vec_isnan)
        if isnothing(firstnonnan)
            Y .= NaN
            return Y
        elseif firstnonnan > 1
            Y[1:firstnonnan-1] .= NaN
        end
        if len < length(X)
            Y[len:end] .= NaN
        end
        if (len - firstnonnan) +1 <= numofele
            return Y.* mean(filter(!isnan, X))
        end
        n = firstnonnan
        summed = sum(filter(!isnan, X[(0:ForwardDelta-1).+n]))
        curr_nans = sum(vec_isnan[(0:ForwardDelta-1).+n])
        curr_nonnans = ForwardDelta - curr_nans #count(x -> x == 0, vec_isnan[(0:ForwardDelta-1).+n])
        if n <= (BackDelta + firstnonnan)
            for n in (0:BackDelta) .+ firstnonnan
                #@info("Loop1")
                #if not NaN
                if !vec_isnan[n+ForwardDelta]
                    summed += X[n+ForwardDelta]
                    curr_nonnans += 1
                    Y[n] = summed / curr_nonnans
                    #if NaN
                else
                    curr_nans += 1
                    if n > 1
                        Y[n] = Y[n-1]
                    else
                        Y[n] = NaN
                    end
                end
            end
        end
        n = BackDelta - firstnonnan + 1
        for n in (BackDelta+firstnonnan+1):(len-ForwardDelta)
            #@info("Loop2")
            curr_nans += vec_isnan[n+ForwardDelta]
            curr_nans -= vec_isnan[n-BackDelta-1]
            curr_nonnans += !vec_isnan[n+ForwardDelta]
            curr_nonnans -= !vec_isnan[n-BackDelta-1]
            if !vec_isnan[n+ForwardDelta]
                summed += X[n+ForwardDelta]
            end
            if !vec_isnan[n-BackDelta-1]
                summed -= X[n-BackDelta-1]
            end
            if curr_nonnans > 0
                Y[n] = summed / curr_nonnans
            else
                Y[n] = NaN
            end
        end
        n = len - ForwardDelta
        for n in len-ForwardDelta+1:len
            #@info("Loop3")
            curr_nans -= vec_isnan[n-BackDelta-1]
            curr_nonnans -= !vec_isnan[n-BackDelta-1]
            if !vec_isnan[n-BackDelta-1]
                summed -= X[n-BackDelta-1]
            end
            if curr_nonnans > 0
                Y[n] = summed / curr_nonnans
            else
                Y[n] = NaN
            end
        end
    end
    return Y
end

"""
    readaerialprofilefromnetcdf(source::String)

Read x, y, z, time, T, U from netcdf (arrays)
"""
function readaerialprofilefromnetcdf(source::String)
    ds = Dataset(source, "r")
    vars = varbyattrib(ds)
    x = ds["x"][:]
    y = ds["y"][:]
    z = ds["z"][:]
    time = ds["time"][:]
    T = ds["T"][:]
    U = ds["U"][:]
    close(ds)
    return x, y, z, time, T, U
end

"""
    readbinarypoints(sourcefile::String, startnr::Int=833, hexperpts::Int=8)

Read binary points file from polyMesh and convert to array
"""
function readbinarypoints(sourcefile::String, startnr::Int=833, hexperpts::Int=8)
    @show("Convert binary constant/polyMesh/points to BC-points-file")
    a=read(sourcefile)

    #start of number in 833
    startnr = 833
    endnr = findfirst(x->x==0x0a, a[833:end])+831
    nrcells = parse(Int, String(a[startnr:endnr]))

    pts = ones(Float64, nrcells, 3)

    idx = endnr+3
    hexperpts = 8

    #convert binary floats to Float64 and fill array
    for i in 1:nrcells
        for j in 1:3
            pts[i,j] = reinterpret(Float64, a[(0:hexperpts-1).+idx])[1]
            idx += (hexperpts)
        end
    end

    return pts
end

"""
    readforcingfromnetcdf(source)

Read data from forcing .nc file
"""
function readforcingfromnetcdf(source)
    ds = Dataset(source, "r")
    heights_out = ds["heights"][:]
    times_out = ds["times"][:]
    u_out = ds["wind"][:,:,:]
    T_out = ds["tair"][:,:]
    close(ds)
    return heights_out, times_out, u_out, T_out
end

"""
    readpoints(source::String)

Read points from file
"""
function readpoints(source::String)
    a = readlines(source)[21:end]
    nrcells = parse(Int, a[1])
    pts = ones(Float64, nrcells, 3)

    a=a[3:end-4]

    for i in 1:length(a)
        tmp = split(a[i][2:end-1], " ")
        pts[i, 1] = parse(Float64, tmp[1])
        pts[i, 2] = parse(Float64, tmp[2])
        pts[i, 3] = parse(Float64, tmp[3])
    end

    return pts
end

"""
    readprofilefromnetcdf(source::String)::DataFrame

Read the concentated profiles from netCDF4 file
"""
function readprofilefromnetcdf(source::String)::DataFrame
    ds = Dataset(source, "r")
    df = DataFrame()
    cols = map(name, varbyattrib(ds))
    for col in cols
        df[:,col] = ds[col][:]
    end
    close(ds)
    return df
end

"""
    readT(source::String)

read temperature values for aerial cross section
"""
function readT(source::String)
    a = readlines(source)[21:end-4]
    nrcells = parse(Int, a[1])
    vals = ones(Float64, nrcells)
    a=a[3:end]

    for i in 1:length(a)
        vals[i] = parse(Float64, a[i])
    end

    return vals
end

"""
    readturbasnetcdf(source::String, perstart::DateTime, perend::DateTime)::DataFrame

Read turbulence data given start and endtime from NetCDF4-file.
"""
function readtotalturbasnetcdf(source::String)::DataFrame
    @info("Reading turbulence data from NetCDF-file")
    ds = Dataset(source, "r")
    #load
    df = DataFrame()
    cols = map(name, varbyattrib(ds))
    for col in cols
        df[:, col] = ds[col][:]
    end
    close(ds)
    return df
end

"""
    readU(source::String)

read wind speeds values for aerial cross section
"""
function readU(source::String)
    a = readlines(source)[21:end-4]
    nrcells = parse(Int, a[1])
    vals = ones(Float64, nrcells, 3)
    a=a[3:end]

    for i in 1:length(a)
        tmp = split(a[i][2:end-1], " ")
        vals[i, 1] = parse(Float64, tmp[1])
        vals[i, 2] = parse(Float64, tmp[2])
        vals[i, 3] = parse(Float64, tmp[3])
    end

    return vals
end

"""
    reducedoubleentries(profile::Array, dist_along::Vector, tolerance::Number=0.01)

Reduce profile lines with multiple entries per dist_along. Take mean instead.
"""
function reducedoubleentries(profile::Array, dist_along::Vector, tolerance::Number=0.01)
    newprof = fill(NaN, size(profile, 1), size(profile, 2))
    newdist = fill(NaN, length(dist_along))
    curridx = 1
    i = 1
    while i <= size(profile, 1)
        if !isnan(dist_along[i])
            idx = findall(x->dist_along[i]-tolerance/2 <= x <= dist_along[i]+tolerance/2, dist_along)
        else
            break
        end
        if size(idx, 1)>1
            for j in 1:size(profile, 2)
                newprof[curridx, j] = mean(profile[idx, j])
            end
            newdist[curridx] = mean(dist_along[idx])
        else
            newprof[curridx, :] = profile[i,:]
            newdist[curridx] = dist_along[i]
        end
        curridx += 1
        i += size(idx, 1)
    end
    return newprof, newdist
end

"""
    reducenonnans(x_idx::Vector, y_idx::Vector, z_idx::Vector)

Return boolean vectors for x-, y-, and z-dimension of array to reduce if too many NaNs are there.
"""
function reducenonnans(x_idx::Vector, y_idx::Vector, z_idx::Vector, nanthresh::Number=0.9)
    tmparray = fill(false, length(unique(x_idx)), length(unique(y_idx)), length(unique(z_idx)))
    for j in 1:length(x_idx)
        tmparray[x_idx[j], y_idx[j], z_idx[j]] = true
    end
    totalentriesperrow = size(tmparray, 2) * size(tmparray, 3)
    keeprows = fill(true, size(tmparray, 1))
    for i in 1:size(tmparray, 1)
        nr_no_entries = count(x->x==0, tmparray[i,:,:])
        if nr_no_entries/totalentriesperrow .> nanthresh
            keeprows[i] = false
        end
    end
    totalentriesperydim = size(tmparray, 1) * size(tmparray, 3)
    keepydim = fill(true, size(tmparray, 2))
    for i in 1:size(tmparray, 2)
        nr_no_entries = count(x->x==0, tmparray[:,i,:])
        if nr_no_entries/totalentriesperydim .> nanthresh
            keepydim[i] = false
        end
    end
    totalentriespercol = size(tmparray, 1) * size(tmparray, 2)
    keepcols = fill(true, size(tmparray, 3))
    for i in 1:size(tmparray, 3)
        nr_no_entries = count(x->x==0, tmparray[:,:,i])
        if nr_no_entries/totalentriespercol .> nanthresh
            keepcols[i] = false
        end
    end

    discard_xidcs = findall(x->x==0, keeprows)
    discard_yidcs = findall(x->x==0, keepydim)
    discard_zidcs = findall(x->x==0, keepcols)
    readline = fill(true, length(x_idx))
    for k in 1:length(x_idx)
        if !isnothing(findfirst(x->x==x_idx[k], discard_xidcs))
            readline[k] = false
        elseif !isnothing(findfirst(x->x==y_idx[k], discard_yidcs))
            readline[k] = false
        elseif !isnothing(findfirst(x->x==z_idx[k], discard_zidcs))
            readline[k] = false
        end
    end

    return keeprows, keepydim, keepcols, readline
end

"""
    saveasnetcdf(T::Array, U::Array, xvals::Vector, yvals::Vector, zvals::Vector, timefolders::Vector, target::String, deflatelvl::Int64=5)

Save DataFrame as netcdf-file
"""
function saveasnetcdf(T::Array, U::Array, xvals::Vector, yvals::Vector, zvals::Vector, timefolders::Vector, target::String, deflatelvl::Int64=5)
    timevals = sort(parse.(Float64, timefolders))
    ds = NCDataset(target, "c")
    defDim(ds, "x", length(xvals))
    defDim(ds, "y", length(yvals))
    defDim(ds, "z", length(zvals))
    defDim(ds, "time", length(timevals))
    defDim(ds, "wind component", 3)
    defVar(ds, "x", xvals, ("x",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "y", yvals, ("y",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "z", zvals, ("z",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "time", timevals, ("time",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "T", T, ("x", "y", "z", "time",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "U", U, ("x", "y", "z", "time","wind component",); shuffle=true, deflatelevel=deflatelvl)
    close(ds)
end

"""
    saveaslinnetcdf(T::Array, U::Array, xvals::Vector, yvals::Vector, zvals::Vector, timefolders::Vector, target::String, deflatelvl::Int64=5)

Save DataFrame as netcdf-file (data linearly as in points and T&U output)
"""
function savelinasnetcdf(T::Array, U::Array, xvals::Vector, yvals::Vector, zvals::Vector, timefolders::Vector, target::String, deflatelvl::Int64=5)
    timevals = sort(parse.(Float64, timefolders))
    ds = NCDataset(target, "c")
    defDim(ds, "length", length(xvals))
    defDim(ds, "time", length(timevals))
    defDim(ds, "wind component", 3)
    defVar(ds, "x", xvals, ("length",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "y", yvals, ("length",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "z", zvals, ("length",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "time", timevals, ("time",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "T", T, ("length", "time",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "U", U, ("length", "time","wind component",); shuffle=true, deflatelevel=deflatelvl)
    close(ds)
end

"""
    saveforcingasnetcdf(frontpts_in::Array, uout_in::Array, Tout_in::Array, totaltimes_in::Array, target::String, deflatelvl::Int64=5)

Save wind and temperature forcing as netcdf.
"""
function saveforcingasnetcdf(frontpts_in::Array, uout_in::Array, Tout_in::Array, totaltimes_in::Array, target::String, deflatelvl::Int64=5)
    ds = NCDataset(target, "c")
    ds.dim["height"] = size(frontpts_in, 1)
    ds.dim["windcomponent"] = size(uout_in, 2)
    ds.dim["time"] = size(uout_in, 3)
    defVar(ds, "heights", frontpts_in[:,3], ("height",); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "times", totaltimes_in, ("time",); shuffle=true, deflatelevel = deflatelvl)
    defVar(ds, "wind", uout_in, ("height", "windcomponent", "time"); shuffle=true, deflatelevel=deflatelvl)
    defVar(ds, "tair", Tout_in, ("height", "time"); shuffle=true, deflatelevel=deflatelvl)
    close(ds)
end


function sortprofileentries(profilealong::String, points::Array, ptstotake::Vector)::Vector
    if !(profilealong in ["x", "y", "z"])
        @error("Variable profilealong has wrong value. x, y, or z allowed!")
        return [1]
    elseif profilealong == "x"
        sortidx = sortperm(points[ptstotake, 1])
    elseif profilealong == "y"
        sortidx = sortperm(points[ptstotake, 2])
    elseif profilealong == "z"
        sortidx = sortperm(points[ptstotake, 3])
    end
    return sortidx
end

"""
    timedepvectoroutstring(p::Array, times::Vector, outdir::String, filename::String)

Create output for time-dependent vector field p (nxmxtimes)
"""
function timedepvectoroutstring(p::Array, times::Vector, outdir::String, filename::String)
    mkpath(outdir)
    for itime in 1:length(times)
        currdir = joinpath(outdir, string(times[itime]))
        mkpath(currdir)
        if ndims(p) == 2
            towrite = vectoroutstring(p[:,itime])
        elseif ndims(p) == 3
            towrite = vectoroutstring(p[:,:,itime])
        else
            @error("No method for writing array with dimensions != 2 or 3! Exiting w/o writing.")
            return;
        end
        write(joinpath(currdir, filename), towrite)
    end
end

"""
    turbflux(data::DataFrame, reyavgtime::Period, p_over_p0::Float64=1013 / 798)::DataFrame

Calculate turbulent flux and Obukhov length
"""
function turbflux(data::DataFrame, reyavgtime::Period, p_over_p0::Float64=1013 / 798)::DataFrame
    leng = size(data, 1)
    fluxout = DataFrame("time" => data.time, "wT" => fill(NaN, leng), "uv" => fill(NaN, leng), "wq" => fill(NaN, leng),
        "uw" => fill(NaN, leng), "vw" => fill(NaN, leng), "uT" => fill(NaN, leng),
        "u_star" => fill(NaN, leng), "uu" => fill(NaN, leng), "vv" => fill(NaN, leng),
        "ww" => fill(NaN, leng), "TT" => fill(NaN, leng), "tke" => fill(NaN, leng), "T_pot_20Hz" => fill(NaN, leng),
        "L_highfreq" => fill(NaN, leng))

    timestep = data.time[2] - data.time[1]
    #check
    checkidx = round(Int, (size(data, 1) - 1) * rand()) + 1
    timestepcheck = data.time[checkidx] - data.time[checkidx-1]
    if timestep != timestepcheck
        @warn("timestep check (sonic data) was not successfull!")
    end
    avgidcs = maximum([1 round(Int, Millisecond(reyavgtime) / Millisecond(timestep))])
    avg_u = movingaverage(data.u, avgidcs)
    avg_v = movingaverage(data.v, avgidcs)
    avg_w = movingaverage(data.w, avgidcs)
    avg_T = movingaverage(data.T, avgidcs)
    avg_T_pot = avg_T .* (p_over_p0^(2 / 7))
    dev_u = data.u .- avg_u
    dev_v = data.v .- avg_v
    dev_w = data.w .- avg_w
    dev_T = data.T .- avg_T
    fluxout.wT = dev_w .* dev_T
    fluxout.uv = dev_u .* dev_v
    fluxout.uw = dev_u .* dev_w
    fluxout.vw = dev_v .* dev_w
    fluxout.uT = dev_u .* dev_T
    fluxout.u_star = (fluxout.uw .^ 2 + fluxout.vw .^ 2) .^ (1 / 4)
    fluxout.uu = dev_u .^ 2
    fluxout.vv = dev_v .^ 2
    fluxout.ww = dev_w .^ 2
    fluxout.TT = dev_T .* dev_T
    fluxout.tke = (fluxout.uu .+ fluxout.vv .+ fluxout.ww) .* 0.5
    fluxout.T_pot_20Hz = data.T .* (p_over_p0^(2 / 7))
    fluxout.L_highfreq = -((fluxout.u_star .^ 3) .* avg_T_pot) ./ fluxout.wT ./ (0.4 * 9.81)#Obukhov-length
    if "h2o" in names(data)
        avg_q = movingaverage(data.h2o, avgidcs)
        dev_q = data.h2o .- avg_q
        fluxout.wq = dev_w .* dev_q
    end
    return fluxout
end

"""
    vectoroutstring(p::Array)

Create String from n x m-array containing vector data (points, u) for writing to file
"""
function vectoroutstring(p::Array)
    outstring = String[]
    if ndims(p) == 1
        push!(outstring, "(\n")
        for j in 1:size(p, 1)
                push!(outstring, string(string(p[j, 1]), " "))
            if j != size(p, 1)
                push!(outstring, "\n")
            else
                push!(outstring, "\n)")
            end
        end
    elseif ndims(p) == 2
        push!(outstring, "(\n( ")
        for j in 1:size(p, 1)
            for i in 1:size(p, 2)
                push!(outstring, string(string(p[j, i]), " "))
            end
            if j != size(p, 1)
                push!(outstring, ")\n( ")
            else
                push!(outstring, ")\n)")
            end
        end
    end
    return join(outstring)
end

end #module
