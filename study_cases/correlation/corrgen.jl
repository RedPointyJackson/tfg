#!/usr/bin/julia --color=yes

for n in 1:30
    header=
    """
    #!/bin/sh
    #\$ -o \$HOME/corrthing_$(n).out -j y
    #\$ -N corrthing_$(n)

    cd tfg/

    ./exe/annealer -p -t 1e6 -T 1e8 -n 100 -i 1e9 -m 100 0.8 corrthing_$(n).net

    """
    open("job_$n","w") do io
        println(io,header)
    end
end
