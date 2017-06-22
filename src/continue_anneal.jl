#!/usr/bin/env julia

# For the files given in ARGS, create all the needed continuation jobs
# for the second phase. That includes all the accessible T2 and all
# the tw₁ starts.

const T2iters = 1e8
const availableT2 = [0.9, 0.8, 0.7, 0.6, 0.5]
const tw1list = map(Int64,[1e4,1e5,1e6,1e7,1e8])

const param = "-m 33 -t 1e4 -T 1e8 -n 10 -L"

function get_fields(file)
    regex = r"quench_T_([0-9.]+)_L_([0-9]+)_([a-z]+)_R([AB]).net"
    fields = match(regex, basename(file)).captures
    length(fields) == 0 && error("Couldn't get metadata from $file")
    T = parse(Float64, fields[1])
    L = parse(Int64, fields[2])
    id = fields[3] * "-R" * fields[4]
    return T, L, id
end

# Return strings containing the continuation jobs for .net file
function get_job_list(file)
    T1, L, id = get_fields(file)
    jobs = String[]
    T2list = [t for t ∈ availableT2 if t<T1]
    for T2 in T2list
        for tw1 in tw1list
            outname = "continue_at_T2_$(T2)_tw1_$(tw1)_for_T1_$(T1)_L_$(L)_id_$(id).net"
            cmd = "./exe/annealer -c $tw1  $param -p -l $L -i $T2iters $T2"
            job = "cat $file | $cmd \"$(outname)\""
            cleanup = "lzma $(outname)"
            push!(jobs, job * "\n" * cleanup * "\n" * "echo 'All done! ($outname)'")
        end
    end
    return jobs
end

jobs = String[]
for file in ARGS
    for j in get_job_list(file)
        push!(jobs,j)
    end
end
info("$(length(jobs)) jobs processed.")

n = 1
for job in jobs
    header=
    """
    #!/bin/sh
    #\$ -o \$HOME/glass_continue_$(n).out -j y
    #\$ -N glass_$(n)

    cd tfg/

    """
    open("continueT2_$n.sh", "w") do io
        println(io,header)
        println(io,job)
    end
    n+=1
end
