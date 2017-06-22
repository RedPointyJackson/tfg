#!/usr/bin/env julia

const Tlist = [0.9, 0.8, 0.7]
const tw1list = map(Int64,[1e4,1e6,1e8])

const reps=8

function rectify_times!(array)
    L = length(array)
    for i in 2:L
        if array[i] == array[i-1]
            for j in i:L
                array[j] += 1
            end
        end
    end
end

function generate_times()
    tws = logspace(log10(1),log10(3e8+10),100)
    tws = map(x->floor(Int64,x), tws)
    tws = vcat(tws, tw1list) |> sort
    rectify_times!(tws)

    tws = unique(tws) |> sort
    ts = logspace(4,8,10)
    ts = map(x->floor(Int64,x), ts) |> sort

    allt = [tw + t for tw in tws, t in [0; ts]] |> vec |> unique |> sort

    writedlm("tiempos_log.dat", allt)
    writedlm("tiempos_lineal.dat",[])
end

if isdir("jobs")
    error("Cowardly refusing to use existing directory 'jobs'")
else
    mkdir("jobs")
end

function generate_jobs(;T=0.9)
    suffix = "quench_T_$(T)"
    mkdir("jobs/$suffix")
    for i in 1:reps
        jobid  = join(rand('a':'z',5))
        open("jobs/$(suffix)/$(jobid).job", "w") do f
            println(f, "/home/sergio/output/$(suffix)_replica_$(i)    Nombre del directorio")
            println(f, "400000001                                Timewindow (maximo numero pasos MC simulacion) ")
            println(f, "0                                        FlagJ                                          ")
            println(f, "0                                        FlagS                                          ")
            println(f, "666                                      Flag_h (now, flag_dummy: not used)             ")
            println(f, "0                                        Flag_continua                                  ")
            println(f, "$T                                      T1                                             ")
            println(f, "$T                                      T2                                             ")
            println(f, "1                                     tw1                                            ")
            println(f, "2                                    tw2                                            ")
            println(f, "0                                        Hfield                                         ")
            println(f, "1                                        TimeHon                                        ")
        end
    end
end


generate_times()


info("Generating jobs")
f = 0
for T in Tlist
    generate_jobs(T=T)
    f+=1
    print("$f*$reps files generated  \r")
end
println()
