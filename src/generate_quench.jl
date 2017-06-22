#!/usr/bin/env julia

const ITERS = 3e8
const Tlist = [0.9, 0.8, 0.7]#, 0.6, 0.5]
const Llist = [20]#[4, 6, 8, 12, 16, 20]

const param = "-m 100 -t 1e4 -T 1e8 -n 10 -L"

if isdir("jobs")
    error("Cowardly refusing to use existing directory 'jobs'")
else
    mkdir("jobs")
end

function write_job_tuples(;T=0.9, iters=3e8, L=8)
    globalid    = join(rand('a':'z',5))
    open("jobs/$(globalid).job", "w") do f
        println(f, "#!/bin/sh")
        println(f, "#\$ -o \$HOME/$(globalid).out -j y")
        println(f, "#\$ -N $(globalid)")
        println(f, "#\$ -pe mp16 16")
        println(f)
        println(f, "cd tfg/")
        println(f)
        for i in 1:16/2
            id    = join(rand('a':'z',5))
            meta  = "quench_T_$(T)_L_$(L)_$(id)"
            jseed = rand(UInt32)
            println(f, "./exe/annealer -j $jseed $param -p -l $L -i $iters $T \"$(meta)_RA.net\"")
            println(f, "./exe/annealer -j $jseed $param -p -l $L -i $iters $T \"$(meta)_RB.net\"")
        end
    end
end

f = 0
for T in Tlist, L in Llist
    write_job_tuples(iters=ITERS, T=T, L=L)
    f+=1
    print("$f files generated  \r")
end
println()
