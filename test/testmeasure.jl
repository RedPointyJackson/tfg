using Base.Test
using DataFrames

@testset "All spins down" begin
    const L = 4
    const meas = 10
    const theT = 1

    ferrofileA = tempname()
    ferrofileB = tempname()
    run(`./exe/annealer -j $(rand(UInt32)) -f -l $L -i 1e3 -m $meas $theT $ferrofileA`)
    run(`./exe/annealer -j $(rand(UInt32)) -f -l $L -i 1e3 -m $meas $theT $ferrofileB`)
    const ONES_BY_TRIPLETS = UInt64(0x1249249249249249)
    ONES = [ONES_BY_TRIPLETS for i in 1:L*L*L] :: Vector{UInt64}

    # Overwrite all the spins with "↑"
    fakefileA = tempname()
    fakefileB = tempname()

    # Copy the header
    ioferroA = open(ferrofileA,"r")
    iofakeA  = open(fakefileA,"w")

    ioferroB = open(ferrofileB,"r")
    iofakeB  = open(fakefileB,"w")

    header_A        = read(ioferroA,Int8,8)
    is_big_endian_A = read(ioferroA,Int64)
    rL_A            = read(ioferroA,Int64)
    ntws_A          = read(ioferroA,Int64)
    tws_A           = read(ioferroA,Int64,ntws_A)
    nts_A           = read(ioferroA,Int64)
    ts_A            = read(ioferroA,Int64,nts_A)
    total_meas_A    = read(ioferroA,Int64)
    T_A             = read(ioferroA,Float64)

    header_B        = read(ioferroB,Int8,8)
    is_big_endian_B = read(ioferroB,Int64)
    rL_B            = read(ioferroB,Int64)
    ntws_B          = read(ioferroB,Int64)
    tws_B           = read(ioferroB,Int64,ntws_B)
    nts_B           = read(ioferroB,Int64)
    ts_B            = read(ioferroB,Int64,nts_B)
    total_meas_B    = read(ioferroB,Int64)
    T_B             = read(ioferroB,Float64)

    write(iofakeA, header_A)
    write(iofakeA, is_big_endian_A)
    write(iofakeA, rL_A)
    write(iofakeA, ntws_A)
    write(iofakeA, tws_A)
    write(iofakeA, nts_A)
    write(iofakeA, ts_A)
    write(iofakeA, total_meas_A)
    write(iofakeA, T_A)

    write(iofakeB, header_B)
    write(iofakeB, is_big_endian_B)
    write(iofakeB, rL_B)
    write(iofakeB, ntws_B)
    write(iofakeB, tws_B)
    write(iofakeB, nts_B)
    write(iofakeB, ts_B)
    write(iofakeB, total_meas_B)
    write(iofakeB, T_B)

    for i in 1:3 # J's
        Js_A = read(ioferroA, UInt64, L*L*L)
        Js_B = read(ioferroB, UInt64, L*L*L)
        write(iofakeA, Js_A)
        write(iofakeB, Js_B)
    end

    # Now, for every measure set, copy it with spins down except for
    # the first iteration.
    flagfirst = true
    for i in 1:total_meas_A
        mcA = read(ioferroA,Int64)
        write(iofakeA, mcA)

        mcB = read(ioferroB,Int64)
        write(iofakeB, mcB)

        spinsA = read(ioferroA,UInt64,L*L*L)
        spinsB = read(ioferroB,UInt64,L*L*L)
        if flagfirst
            write(iofakeA, ONES)
            write(iofakeB, ONES)
            flagfirst = false
        else
            write(iofakeA, zeros(UInt64,L*L*L))
            write(iofakeB, zeros(UInt64,L*L*L))
        end
    end

    footerA = read(ioferroA,Int8,3)
    write(iofakeA,footerA)

    footerB = read(ioferroB,Int8,3)
    write(iofakeB,footerB)

    close(ioferroA); close(iofakeA)
    close(ioferroB); close(iofakeB)

    # Measure things
    measfile = tempname()

    run(pipeline(`./exe/measure -cCemq $fakefileA $fakefileB`,stdout=measfile))
    data = readtable(measfile)

    # The first iteration is completelly uncorrelated with everything
    # else except itself.
    # The rest are completelly correlated.
    Cdata = data[data[:observable] .== "Correlation",:]
    for i in 1:size(Cdata,1)
        mc = Cdata[:mc][i]
        t  = Cdata[:parameter][i]
        C  = Cdata[:value][i]
        if mc == 1 && t≠0
            @test C ≈ -1
        else
            @test C ≈ 1
        end
    end

    # The spatial correlation is complete.
    sCdata = data[data[:observable] .== "Spatial_correlation",:]
    @test all(sCdata[:value] .== 1)


    # The energy is always -1.
    Edata = data[data[:observable] .== "Energy",:]
    @test all(Edata[:value] .== -1)

    # The magnetization is always -1 (except for the first iteration).
    Mdata = data[data[:observable] .== "Magnetization",:]
    for i in 1:size(Mdata,1)
        if Mdata[:mc][i] == 1
            @test Mdata[:value][i] == 1
        else
            @test Mdata[:value][i] == -1
        end
    end

    # Overlap is absolute.
    Qdata = data[data[:observable] .== "Overlap",:]
    @test all(Qdata[:value] .== 1)

    # The only ts are the intended ones.
    U_ts  = Cdata[:parameter] |> unique
    @test U_ts == vcat(ts_A)
    @test U_ts == vcat(ts_B)

    # mc = tw + t
    tw = data[:mc]
    tw_plus_t = [tw+t for tw in tws_A, t in vcat(0,ts_A)]
    tw_plus_t = [tw+t for tw in tws_B, t in vcat(0,ts_B)]
    @test tw |> unique |> sort == tw_plus_t |> unique |> sort
end
