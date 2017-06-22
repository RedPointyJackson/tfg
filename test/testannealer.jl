using Base.Test

function test_output(L,iters,m)

    const mint = 2
    const maxt = 20
    run(`./exe/annealer -l $L -t $mint -T $maxt -i $iters -m $m 42 debugfile`);
    io = open("debugfile")

    # Firstly, there should be some fields with information.
    header = map(Char, read(io,Int8,8)) |> join
    @test header == "ANNEALER"

    is_big_endian = read(io,Int64)
    @test is_big_endian ∈ [0,1]

    rL = read(io,Int64)
    @test rL == L

    ntws = read(io,Int64)
    tws = read(io,Int64,ntws)
    @test minimum(tws) == 1
    @test maximum(tws) == iters
    @test all(1≤tw≤iters for tw in tws)

    nts = read(io,Int64)
    ts = read(io,Int64,nts)
    @test minimum(ts) == mint
    @test maximum(ts) == maxt
    @test all(mint≤t≤maxt for t in ts)

    total_meas = read(io, Int64)
    @test total_meas ≤ ntws*nts

    T = read(io, Float64)
    @test T ≈ 42

    Jup = read(io,UInt64,L*L*L)
    Jright = read(io,UInt64,L*L*L)
    Jfront = read(io,UInt64,L*L*L)

    for i in 1:L*L*L # Check J∈0,1
        @test all(J∈('0','1') for J in oct(Jup[i]))
        @test all(J∈('0','1') for J in oct(Jright[i]))
        @test all(J∈('0','1') for J in oct(Jfront[i]))
    end

    mcs = Int64[]
    for m in 1:total_meas
        mc = read(io,Int64)
        push!(mcs,mc)
        spins = read(io,UInt64,L*L*L)
        for i in 1:L*L*L # Check σ∈0,1
            @test all(σ∈('0','1') for σ in oct(spins[i]))
        end
    end

    @test minimum(mcs) == 1
    @test maximum(mcs) == maximum(tws)+maximum(ts)
    @test iters ∈ mcs


    footer = map(Char, read(io,Int8,3)) |> join
    @test footer == "END"

    @test eof(io)

    close(io)
    rm("debugfile")
end


@testset "Output format" begin
    const iters = 10
    for m in 2:iters
        for L in [2,8]
            test_output(L, iters, m)
        end
    end
end


@testset "tw and t specification" begin

    function getfields(cmd)
        run(cmd);
        io = open("debugfile")
        header = map(Char, read(io,Int8,8)) |> join
        is_big_endian = read(io,Int64)
        L = read(io,Int64)
        ntws = read(io,Int64)
        tws = read(io,Int64,ntws)
        nts = read(io,Int64)
        ts = read(io,Int64,nts)
        total_meas = read(io, Int64)
        T = read(io, Float64)
        Jup = read(io,UInt64,L*L*L)
        Jright = read(io,UInt64,L*L*L)
        Jfront = read(io,UInt64,L*L*L)
        mcs = Int64[]
        for m in 1:total_meas
            mc = read(io,Int64)
            push!(mcs,mc)
            spins = read(io,UInt64,L*L*L)
        end
        return mcs,tws,ts
    end

    function glasslogspace(from,to,N)
        v = logspace(log10(from),log10(to),N)
        v = floor(Int64,v)
        v[1] = from
        v[end] = to
        return unique(v)
    end

    glasslinspace(from,to,N) = floor(Int64,linspace(from,to,N)) |> unique

    sortunique(x) = x |> unique |> sort

    for m in [10,50,150], M in m+[10,50,150]
        # Normal scale for t and tws
        mcs,tws,ts = getfields(`./exe/annealer -m 84 -i 423 -t $m -T $M -n 15 42 debugfile`)
        @test tws == glasslogspace(1,423,84)
        @test ts  == glasslinspace(m,M,15)
        tw_plus_t = [tw+t for tw in tws, t in vcat(0,ts)]
        @test sortunique(mcs) == sortunique(tw_plus_t)

        # Specify log scale for t
        mcs,tws,ts = getfields(`./exe/annealer -L -m 84 -i 423 -t $m -T $M -n 15 42 debugfile`)
        @test tws == glasslogspace(1,423,84)
        @test ts  == glasslogspace(m,M,15)
        tw_plus_t = [tw+t for tw in tws, t in vcat(0,ts)]
        @test sortunique(mcs) == sortunique(tw_plus_t)
        # Only a t
        mcs,tws,ts = getfields(`./exe/annealer -m 84 -i 423 -t $m -n 1 42 debugfile`)
        @test tws == glasslogspace(1,423,84)
        @test ts  == [m]
        tw_plus_t = [tw+t for tw in tws, t in vcat(0,ts)]
        @test sortunique(mcs) == sortunique(tw_plus_t)
        # No t
        mcs,tws,ts = getfields(`./exe/annealer -m 84 -i 423 -n 0 42 debugfile`)
        @test tws == glasslogspace(1,423,84)
        @test ts  == [0]
        @test mcs == tws
    end
end

@testset "Output format" begin
    const iters = 10
    for m in 2:iters
        for L in [2,8]
            test_output(L, iters, m)
        end
    end
end


@testset "Invalid input" begin
    # Normaly, return 0
    cmd = `./exe/annealer -i 10 -m 2 42 debugfile`
    io, proc = open(pipeline(cmd))
    wait(proc)
    @test proc.exitcode == 0
    # You can't do more measurements than iters
    cmd = `./exe/annealer -i 20 -m 21 42 debugfile`
    io, proc = open(pipeline(cmd,stderr=DevNull))
    wait(proc)
    @test proc.exitcode != 0
    # But you can measure every iter
    cmd = `./exe/annealer -i 20 -m 20 42 debugfile`
    io, proc = open(pipeline(cmd,stderr=DevNull))
    wait(proc)
    @test proc.exitcode == 0
    # You need to provide at least a temperature
    cmd = `./exe/annealer -i 20 -m 2 debugfile`
    io, proc = open(pipeline(cmd,stderr=DevNull))
    wait(proc)
    @test proc.exitcode != 0
    # You can't provide more than a temperature
    cmd = `./exe/annealer -i 20 -m 2 4 4 debugfile`
    io, proc = open(pipeline(cmd,stderr=DevNull))
    wait(proc)
    @test proc.exitcode != 0
    # It must be a positive T
    cmd = `./exe/annealer -i 20 -m 2 -2 debugfile`
    io, proc = open(pipeline(cmd,stderr=DevNull))
    wait(proc)
    @test proc.exitcode != 0
    # Don't even try to cheat argp
    cmd = `./exe/annealer -i 20 -m 2 -- -2 debugfile`
    io, proc = open(pipeline(cmd,stderr=DevNull))
    wait(proc)
    @test proc.exitcode != 0
    # measurements ≥ 1
    cmd = `./exe/annealer -i 20 -m 1 -- -2 debugfile`
    io, proc = open(pipeline(cmd,stderr=DevNull))
    wait(proc)
    @test proc.exitcode != 0
    # iters ≥ 1
    cmd = `./exe/annealer -i 1 -m 0 -- -2 debugfile`
    io, proc = open(pipeline(cmd,stderr=DevNull))
    wait(proc)
    @test proc.exitcode != 0

    rm("debugfile")
end

@testset "ferromagnet" begin
    run(`./exe/annealer -f -i 20 -m 10 42 debugfile`);
    io = open("debugfile")

    header = map(Char, read(io,Int8,8)) |> join
    is_big_endian = read(io,Int64)
    L = read(io,Int64)
    ntws = read(io,Int64)
    tws = read(io,Int64,ntws)
    nts = read(io,Int64)
    ts = read(io,Int64,nts)
    total_meas = read(io, Int64)
    T = read(io, Float64)

    Jup = read(io,UInt64,L*L*L)
    Jright = read(io,UInt64,L*L*L)
    Jfront = read(io,UInt64,L*L*L)

    # 0 001 ⋯ 001 001 001
    const ONES_BY_TRIPLETS = 0x1249249249249249

    for i in 1:L*L*L # Check J∈0,1
        @test Jup[i]    == ONES_BY_TRIPLETS
        @test Jright[i] == ONES_BY_TRIPLETS
        @test Jfront[i] == ONES_BY_TRIPLETS
    end

    close(io)
    rm("debugfile")
end


@testset "J seeding" begin

    function getJ(seed)
        run(`./exe/annealer -j $seed -i 20 -m 10 42 debugfile`);
        io = open("debugfile")

        header = map(Char, read(io,Int8,8)) |> join
        is_big_endian = read(io,Int64)
        L = read(io,Int64)
        ntws = read(io,Int64)
        tws = read(io,Int64,ntws)
        nts = read(io,Int64)
        ts = read(io,Int64,nts)
        total_meas = read(io, Int64)
        T = read(io, Float64)

        JupA = read(io,UInt64,L*L*L)
        JrightA = read(io,UInt64,L*L*L)
        JfrontA = read(io,UInt64,L*L*L)

        close(io)
        rm("debugfile")

        return JupA, JrightA, JfrontA
    end

    seed = rand(UInt32)

    @test getJ(seed) == getJ(seed) == getJ(seed) == getJ(seed)
    @test getJ(seed) != getJ(seed*2) != getJ(seed*4) != getJ(seed*5)

end
