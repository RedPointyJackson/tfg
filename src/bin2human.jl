#!/usr/bin/julia --color=yes

# Usage: call it with .net files as ARGS.

function parse_file(file)

    io = open(file,"r")
    println(file,":")

    header = read(io,Int8,8)
    @assert join(map(Char,header),"") == "ANNEALER"
    IS_BIG_ENDIAN = ENDIAN_BOM == 0x01020304

    is_big_endian = read(io,Int64)
    @assert is_big_endian == IS_BIG_ENDIAN

    L = read(io,Int64)
    V = L^3

    nmeas = read(io,Int64)
    tws = read(io,Int64,nmeas)

    nts = read(io,Int64)
    ts = read(io,Int64,nts)

    total_meas = read(io, Int64)

    T = read(io, Float64)

    Jup    = read(io,UInt64,V)
    Jright = read(io,UInt64,V)
    Jfront = read(io,UInt64,V)

    println("    L             : ", L)
    print("    tw's          : ")
    @printf("∈ [%.4g, %.4g] in %d steps\n",
            minimum(tws), maximum(tws), length(tws))
    print("    t's           : ")
    @printf("∈ [%.4g, %.4g] in %d steps\n",
            minimum(ts), maximum(ts), length(ts))
    println("    Endianness    : ", is_big_endian==1? "Big Endian":"Little Endian")
    println("    T             : ", T)
    println("    J hash        : ",
            hex(hash(Jup) + hash(Jright) + hash(Jfront)))


    mcs = Int64[]

    for i in 1:total_meas
        try
            mc    = read(io,Int64)
            spins = read(io,UInt64,V)
            push!(mcs,mc)
        catch e
            if e == EOFError()
                warn("EOF was reached.")
                print()
                print_with_color(:red, "Latest mc's: ")
                print_with_color(:red, "$(mcs[end-0]), ")
                print_with_color(:red, "$(mcs[end-1]), ")
                print_with_color(:red, "$(mcs[end-2]), ")
                print_with_color(:red, "$(mcs[end-3])\n")
                break
            else
                throw(e)
            end
        end
    end

    mcdigest =
    @sprintf("mc steps ∈ (%.4g,%.4g), %.4g values of %.4g\n",
             minimum(mcs), maximum(mcs), length(mcs), total_meas)

    print_with_color(:blue, mcdigest)

    close(io)
end

for file in ARGS
    parse_file(file)
end
