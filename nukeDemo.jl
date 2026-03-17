using DelimitedFiles, Cbc, JuMP, SparseArrays

H = vec(readdlm("data/channel_heights.txt"))

K = [
300 140 40
]


function constructA(H,K)
    n = length(H)
    A = zeros(Float64, n, n)

    k0 = K[1]
    k1 = K[2]
    k2 = K[3]

    for i in 1:n
        A[i,i] = k0
        if i-1 >= 1
            A[i,i-1] = k1
        end
        if i+1 <= n
            A[i,i+1] = k1
        end
        if i-2 >= 1
            A[i,i-2] = k2
        end
        if i+2 <= n
            A[i,i+2] = k2
        end
    end

    return A
end

function solveIP(H, K)
    h = length(H)
    myModel = Model(Cbc.Optimizer)
    # If your want ot use GLPK instead use:
    #myModel = Model(GLPK.Optimizer)

    A = constructA(H,K)

    @variable(myModel, x[1:h], Bin )
    @variable(myModel, R[1:h] >= 0 )

    @objective(myModel, Min, sum(x[j] for j=1:h) )

    @constraint(myModel, [j=1:h],R[j] >= H[j] + 10 )
    @constraint(myModel, [i=1:h],R[i] == sum(A[i,j]*x[j] for j=1:h) )

    optimize!(myModel)

    if termination_status(myModel) == MOI.OPTIMAL
        xval = JuMP.value.(x)
        Rval = JuMP.value.(R)
        bomb_idx = findall(v -> v > 0.5, xval)
        bomb_dist = 250 .* (bomb_idx .- 1)

        writedlm("data/bomb_positions.txt", bomb_dist)
        writedlm("data/bomb_indices.txt", bomb_idx)
        writedlm("data/R_values.txt", Rval)
    else
        println("Optimize was not succesful. Return code: ", termination_status(myModel))
    end
end

solveIP(H,K)
