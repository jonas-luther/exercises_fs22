function shift(N, n, s::Int64)
    mask_lo =  2^N - 1
    mask_hi = mask_lo << N

    s = s << n
    println(bitstring(s))
    s |= (s & mask_hi) >> (N - 1)
    println(bitstring(s)) 
    s &= mask_lo
    println(bitstring(s))
    
end  