#!/usr/bin/julia

function time_parallelism(N; L=8, iters=1e4)
    cmd = "./exe/annealer -m 2 -n 2 -i $iters -l $L 1 > /dev/null"
    # Create shell script
    open("shell.sh", "w") do io
        for i in 1:N
            println(io, cmd * " &")
        end
        println(io, "wait")
    end
    # Run it and time it
    val, time = @timed run(`sh shell.sh`)
    rm("shell.sh")
    return time
end


println("N,L,t")
for L in [8,12,20,30]
    for N in 1:20
        println(N,"\,",L,"\,",time_parallelism(N,L=L))
    end
end
