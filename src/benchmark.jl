#!/usr/bin/julia --color=yes
using ProgressMeter

    noflags = filter(x->'-' ∉x, ARGS)
if length(noflags) == 1
    info("Will use $(noflags[1]) as L.")
    const L        = parse(Int64,noflags[1])
else
    const L = 8
    info("Defaulting to L=$L")
end

function run_annealer(iters; T=1)
    cmd = `./exe/annealer -l $L -i $iters -m 2 -n 2 1 /tmp/debuganneal`
    val, time = @timed run(cmd)
    return time
end

if "-f" ∉ ARGS
    tt = linspace(1e2,1e5,10) |> floor
    times = @showprogress "Benchmarking..." [run_annealer(t) for t in tt]
else
    info("Skipping measurement.")
    tt = linspace(1e2,1e5,10) |> floor
    times = 0*tt
end

function fit_to_line(x, y, yerr=ones(x))
    L = length(x)

    br(x,σ) = sum(  x./(σ.^2) / length(x) )
    brone = br(ones(L), yerr)
    brx   = br(x      , yerr)
    brxx  = br(x.*x   , yerr)
    bry   = br(y      , yerr)
    brxy  = br(x.*y   , yerr)

    m_mean = ( brone*brxy - brx*bry )/
            ( brone*brxx - brx*brx )
    n_mean = mean(y) - m_mean * mean(x)

    covmatrix = zeros(2,2)

    D = brxx*brone - brx*brx

    covmatrix[2,2] = +1/(L*D) * brxx
    covmatrix[1,2] = -1/(L*D) * brx
    covmatrix[2,1] = -1/(L*D) * brx
    covmatrix[1,1] = +1/(L*D) * brone

    σ = diag(covmatrix) |> sqrt

    return m_mean, σ[1], n_mean, σ[2]
end


μm, σm, μn, σn = fit_to_line(tt,times)

@printf("t ∼ %.2lf ± %.2lf × mcsteps × L³  ns\n",
        1e9*μm/L^3, 1e9*σm/L^3)

println()

# Based in empirical data:
μm_memento30 = 126.99*L^3 / 1e9
μm_memento100 = 281.37*L^3 / 1e9

println("Approximate iterations (L=$L):")

print_with_color(:bold,"                                           Memento  \n")
print_with_color(:bold,"        Time         This PC        30 jobs        100 jobs\n")
@printf("    In 1 minute      %.1e        %.1e        %.1e\n" , 60/μm        , 60/μm_memento30        , 60/μm_memento100        )
@printf("    In 1 hour        %.1e        %.1e        %.1e\n" , 3600/μm      , 3600/μm_memento30      , 3600/μm_memento100      )
@printf("    In 1 day         %.1e        %.1e        %.1e\n" , 24*3600/μm   , 24*3600/μm_memento30   , 24*3600/μm_memento100   )
@printf("    In 3 days        %.1e        %.1e        %.1e\n" , 3*24*3600/μm , 3*24*3600/μm_memento30 , 3*24*3600/μm_memento100 )
@printf("    In 1 week        %.1e        %.1e        %.1e\n" , 7*24*3600/μm , 7*24*3600/μm_memento30 , 7*24*3600/μm_memento100 )
