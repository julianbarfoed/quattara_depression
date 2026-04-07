using DelimitedFiles, Cbc, JuMP, SparseArrays

H = vec(readdlm("./data/channel_heights.txt"))

K = [300, 140, 40]
CHD = 10

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

function solveSmoothIP(H, K, CHD)
    n = length(H)
    myModel = Model(Cbc.Optimizer)
    set_optimizer_attribute(myModel, "seconds", 300.0)

    A = constructA(H, K)

    @variable(myModel, x[1:n], Bin)          
    @variable(myModel, R[1:n] >= 0)          
    @variable(myModel, z[1:n] >= 0)  # absolute value calculation

    # maximum smoothness
    @objective(myModel, Min, sum(z[i] for i in 1:n))

    # Depth definitiion
    @constraint(myModel, [i=1:n], R[i] >= H[i] + CHD)

    # Depth accumulation definition
    @constraint(myModel, [i=1:n], R[i] == sum(A[i,j] * x[j] for j in 1:n))

    # absolute value constraint z[i] = | Ri-Hi-CHD |
    @constraint(myModel, [i=1:n], z[i] >= R[i] - (H[i] + CHD))
    
    # add constraint that no adjacent bomb
    @constraint(myModel, [i=1:n-1], x[i] + x[i+1] <= 1)
    
    optimize!(myModel)

    if termination_status(myModel) == MOI.OPTIMAL || termination_status(myModel) == MOI.TIME_LIMIT
        xval = JuMP.value.(x)
        Rval = JuMP.value.(R)
        bomb_idx = findall(v -> v > 0.5, xval)
        bomb_dist = 250 .* (bomb_idx .- 1)

        writedlm("no_adjacent_smooth_bomb_positions.txt", bomb_dist)
        writedlm("no_adjacent_smooth_R_values.txt", Rval)
    
        println("Total nukes used: $(length(bomb_idx)) Nukes")
    else
        println("Optimize was not succesful. Return code: ", termination_status(myModel))
    end
end

solveSmoothIP(H, K, CHD)