#=
xxxdocumentation
script for openFOAM-LES-runs: Concentating the profiles
written by the simulation to time folders -> one file
=#

print("JULIA: Concentating profiles\n")

using CSV, DataFrames, NCDatasets, Dates, ProgressMeter

#path of simulation
casedir = "/home/haugened/Documents/openfoam/real_run_state/cscs/second_run/"

#load functions
duerrpath = "/home/haugened/Documents/openfoam/duerr_les/"
srcpath = joinpath(duerrpath, "scripts", "src", "functions.jl")
if !(@isdefined func)
    include(srcpath)
    #import .func
    println("Including functions")
end

#path to inital profile time folders
initprofdir2 = joinpath(casedir, "postProcessing", "surfaces", "yNormal")

#names of single profile files
proffiles2 = ["T_yNormal.raw", "U_yNormal.raw"]

#logfile
logfiledir2 = joinpath(casedir, "log.of_concentate_surfaces_julia")

#regex matching the desired folders
folderregex = r"(\b\d+\b)|(\b\d+\.\d+\b)"

#where to write the final files
writeloc2 = initprofdir2

#file appendix to "yNormal"
file_appendix = "_time8th"

#concentate surfaces
inittimefolders2 = readdir(initprofdir2)
inittimefolders2 = inittimefolders2[occursin.(folderregex, inittimefolders2)]
if length(inittimefolders2) < 1
    exit()
end

#sort it
init2parsed = parse.(Float64, inittimefolders2)
init2sortperm = sortperm(init2parsed)
inittimefolders2 = inittimefolders2[init2sortperm]

#cut it, because too Long
#inittimefolders2 = inittimefolders2[1:54512]
inittimefolders2 = inittimefolders2[1:8:end]
time_idx = func.createtimeidcs(inittimefolders2)

#remove old log file
rm(logfiledir2, force=true)

logfile2 =  open(logfiledir2,"a")

write(logfile2, string("number of folders to work on: ", length(inittimefolders2), "\n"));
write(logfile2, string("number of profiles in each folder: ", length(proffiles2), "\n\n"));

write(logfile2, "Starting concentation...\n\n");
close(logfile2)

pts1 = func.readpoints(joinpath(initprofdir2, "points1"))
smallerzeros1 = findall(x->x<0, pts1)
pts1[smallerzeros1] .= 0.0
pts2 = func.readpoints(joinpath(initprofdir2, "points2"))
smallerzeros2 = findall(x->x<0, pts2)
pts2[smallerzeros2] .= 0.0

#define when to switch between pts1 and pts2
firstnewpts = findfirst(x->x>47, parse.(Float64, inittimefolders2))

#pts_round = round.(pts, digits=2)

#(x_idx, y_idx, z_idx, unique_x, unique_y, unique_z) = func.createxyzidcs(pts_round)

#select x filter
filterx1 = fill(true, size(pts1, 1)) #7 .<=  pts1[:,1] .&& mod.(round.(pts1[:,1], digits=2) .* 10, 1) .== 0 #fill(true, size(pts1, 1))
filtery1 = fill(true, size(pts1, 1))
filterz1 = fill(true, size(pts1, 1)) #.&& mod.(collect(1:length(unique_z)), 2) .== 0
filter1_total = filterx1 .&& filtery1 .&& filterz1

filterx2 = fill(true, size(pts2, 1)) #7 .<=  pts1[:,1] .&& mod.(round.(pts1[:,1], digits=2) .* 10, 1) .== 0 #fill(true, size(pts2, 1))
filtery2 = fill(true, size(pts2, 1))
filterz2 = fill(true, size(pts2, 1)) #.&& mod.(collect(1:length(unique_z)), 2) .== 0
filter2_total = filterx2 .&& filtery2 .&& filterz2

#sort points
pts1df = DataFrame(x=pts1[filter1_total,1], y=pts1[filter1_total,2], z=pts1[filter1_total,3])
pts1sort = sortperm(pts1df, [:x, :z])
pts2df = DataFrame(x=pts2[filter2_total,1], y=pts2[filter2_total,2], z=pts2[filter2_total,3])
pts2sort = sortperm(pts2df, [:x, :z])

filter_total_count = maximum([count(filter1_total), count(filter2_total)])

Ttmp = fill(NaN, size(pts1, 1))
Utmp = fill(NaN, size(pts1, 1), 3)
Tarray = fill(NaN, filter_total_count, length(unique(time_idx)))
Uarray = fill(NaN, filter_total_count, length(unique(time_idx)), 3)

@showprogress for i in 1:length(inittimefolders2)
    if i < firstnewpts
        srt = pts1sort
        fltr = filter1_total
    else
        srt = pts2sort
        fltr = filter2_total
    end
    Ttmp[:] = func.readT(joinpath(initprofdir2, inittimefolders2[i], "T"))
    Utmp[:,:] = func.readU(joinpath(initprofdir2, inittimefolders2[i], "U"))
    Tarray[:, i] = Ttmp[fltr][srt]
    Uarray[:,i,:] = Utmp[fltr, :][srt, :]

    if i%100 == 0
        logfileloop =  open(logfiledir2,"a")
        write(logfileloop, string(now(), ": ", "timestep ", i, "/", length(inittimefolders2), "\n"));
        close(logfileloop)
    end
end

#save arrays to netcdf

@info("writing file...")

#func.saveasnetcdf(Tarray, Uarray, unique_x, unique_y, unique_z, inittimefolders2, joinpath(writeloc2, string(proffiles2[1][3:end-3], "nc")))
func.savelinasnetcdf(Tarray, Uarray, pts1df[pts1sort, 1], pts1df[pts1sort, 2], pts1df[pts1sort, 3], inittimefolders2, joinpath(writeloc2, string(proffiles2[1][3:end-4], file_appendix, ".nc")))

@info("file written.")

write(logfile2, string("\n Done.\n"));

close(logfile2)

#rm.(joinpath.(initprofdir2, inittimefolders2), recursive=true)
#rm(joinpath(initprofdir2, "points"))