
module Ex5_1
function main()
    vup = 1 / sqrt(2) * [0 1]
    vdown = 1 / sqrt(2) * [1 0]
    Aup = [0 1; 0 0]
    Adown = [0 0; 1 0]
    uup = [1, 0]
    udown = [0, 1]

    vup * Adown * Aup * Adown * Aup * udown
end
end

module Ex5_2
function main()
    
end
end