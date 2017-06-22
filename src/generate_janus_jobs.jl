#!/usr/bin/env julia

const T1list = [0.9, 0.8, 0.7]
const T2list = [0.8, 0.7, 0.6, 0.5]
const tw1list = map(Int64,[1e4,1e5,1e6,1e7,1e8])

const reps=1

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

    for tw1 in tw1list
        newlogs = logspace(log10(tw1), log10(3e8), 50)
        tws = vcat(tws, newlogs)
    end

    endlogs = logspace(log10(2e8), log10(3e8), 50)
    tws = vcat(tws, endlogs)

    tws = map(x->floor(Int64,x), tws) |> sort
    rectify_times!(tws)

    tws = vcat(tw1list, tws) |> unique |> sort
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

function generate_conts(;T1=0.9, T2=0.6, tw1=10_000)
    suffix = "T1_$(T1)_T2_$(T2)_tw1_$(tw1)"
    mkdir("jobs/$suffix")
    for i in 1:reps
        jobid  = join(rand('a':'z',10))
        open("jobs/$(suffix)/$(jobid).job", "w") do f
            println(f, "/home/sergio/output/$(suffix)_$(jobid)_replica_$(i)    Nombre del directorio")
            println(f, "400000001                                Timewindow (maximo numero pasos MC simulacion) ")
            println(f, "0                                        FlagJ                                          ")
            println(f, "0                                        FlagS                                          ")
            println(f, "666                                      Flag_h (now, flag_dummy: not used)             ")
            println(f, "0                                        Flag_continua                                  ")
            println(f, "$T1                                      T1                                             ")
            println(f, "$T2                                      T2                                             ")
            println(f, "$tw1                                     tw1                                            ")
            println(f, "200000000                                tw2                                            ")
            println(f, "0                                        Hfield                                         ")
            println(f, "1                                        TimeHon                                        ")
        end
    end
end

function generate_quenchs(;T=0.9)
    suffix = "quench_T_$(T)"
    mkdir("jobs/$suffix")
    for i in 1:reps
        jobid  = join(rand('a':'z',10))
        open("jobs/$(suffix)/$(jobid).job", "w") do f
            println(f, "/home/sergio/output/$(suffix)_$(jobid)_replica_$(i)    Nombre del directorio")
            println(f, "400000001                                Timewindow (maximo numero pasos MC simulacion) ")
            println(f, "0                                        FlagJ                                          ")
            println(f, "0                                        FlagS                                          ")
            println(f, "666                                      Flag_h (now, flag_dummy: not used)             ")
            println(f, "0                                        Flag_continua                                  ")
            println(f, "$T                                      T1                                             ")
            println(f, "$T                                      T2                                             ")
            println(f, "1                                       tw1                                            ")
            println(f, "2                                        tw2                                            ")
            println(f, "0                                        Hfield                                         ")
            println(f, "1                                        TimeHon                                        ")
        end
    end
end


generate_times()


info("Generating jobs")
f = 0
for T1 in T1list
    generate_quenchs(T=T1)
    for T2 in T2list
        T2 ≥ T1 && continue
        for tw1 in tw1list
            generate_conts(T1=T1, T2=T2, tw1=tw1)
            f+=1
        end
    end
    println()
end
println("$f*$reps files generated")
