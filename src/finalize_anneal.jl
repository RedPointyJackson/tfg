#!/usr/bin/env julia

# For the files given in ARGS, create all the needed finalization jobs
# for the third phase. That is, continue them at T₁ until tw=3e8.

const param = "-m 33 -t 1e4 -T 1e8 -n 10 -L"

function get_fields(file)
    regex = r"continue_at_T2_([0-9.]+)_tw1_([0-9]+)_for_T1_([0-9.]+)_L_([0-9]+)_id_([a-z]+-R[AB]).net"
    fields = match(regex, basename(file)).captures
    length(fields) == 0 && error("Couldn't get metadata from $file")
    T2 = parse(Float64, fields[1])
    tw1 = parse(Int64, fields[2])
    T1 = parse(Float64, fields[3])
    L = parse(Int64, fields[4])
    id = fields[5]
    return T1, T2, tw1, L, id
end

# Return the string containing the finalization jobs for .net file
function get_the_job(file)
    T1, T2, tw1, L, id = get_fields(file)
    outname = "finish_T2_$(T2)_tw1_$(tw1)_from_T1_$(T1)_L_$(L)_id_$(id).net"
    needediters = 3e8 - (tw1+1e8)
    cmd = "./exe/annealer -c $(tw1+1e8) $param -p -l $L -i $needediters $T1"
    job = "cat $file | $cmd \"$(outname)\""
    cleanup = "lzma $(outname)"
    return job * "\n" * cleanup * "\n" * "echo 'All done! ($outname)'"
end

jobs = String[]
for file in ARGS
    push!(jobs,get_the_job(file))
end
info("$(length(jobs)) jobs processed.")

n = 1
for job in jobs
    header=
    """
    #!/bin/sh
    #\$ -o \$HOME/glass_finish_$(n).out -j y
    #\$ -N glass_$(n)

    cd tfg/

    """
    open("finish_$n.sh", "w") do io
        println(io,header)
        println(io,job)
    end
    n+=1
end
