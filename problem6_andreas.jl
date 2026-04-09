<<<<<<< HEAD
using DelimitedFiles, Cbc, JuMP, SparseArrays, LinearAlgebra

H = vec(readdlm("./data/channel_heights.txt"))


CHD = 10

function constructA(H)
    n = length(H)
    A = zeros(Float64, n, n,4)
    K1 = [300, 140, 40]
    K2 = [500, 230, 60]
    K3 = [1000, 400, 70]

    for i in 1:n
        A[i,i,2] = K1[1]
        A[i,i,3] = K2[1]
        A[i,i,4] = K3[1]
        if i-1 >= 1
            A[i,i-1,2] = K1[2]
            A[i,i-1,3] = K2[2]
            A[i,i-1,4] = K3[2]
        end
        if i+1 <= n
            A[i,i+1,2] = K1[2]
            A[i,i+1,3] = K2[2]
            A[i,i+1,4] = K3[2]
        end
        if i-2 >= 1
            A[i,i-2,2] = K1[3]
            A[i,i-2,3] = K2[3]
            A[i,i-2,4] = K3[3]
        end
        if i+2 <= n
            A[i,i+2,2] = K1[3]
            A[i,i+2,3] = K2[3]
            A[i,i+2,4] = K3[3]
        end
    end

    return A
end

function solveSmoothIP(H, CHD)
    n = length(H)
    myModel = Model(Cbc.Optimizer)
    set_optimizer_attribute(myModel, "seconds", 300.0)

    A = constructA(H)

    @variable(myModel, x[1:n,1:5], Bin)
    @variable(myModel, R[1:n] >= 0)
=======
using DelimitedFiles, Cbc, JuMP, SparseArrays

H = vec(readdlm("./data/channel_heights.txt"))
H_part = H[1:250]

# Rows are settings, colu mns are [k0 k1 k2]
K = [
     300 140 40
     500 230 60
    1000 400 70
]
CHD = 10

function constructA3D(H, K)
    n = length(H)
    s = size(K, 1)
    A3 = zeros(Float64, n, n, s)

    for setting in 1:s
        k0 = K[setting, 1]
        k1 = K[setting, 2]
        k2 = K[setting, 3]

        for i in 1:n
            A3[i, i, setting] = k0
            if i - 1 >= 1
                A3[i, i - 1, setting] = k1
            end
            if i + 1 <= n
                A3[i, i + 1, setting] = k1
            end
            if i - 2 >= 1
                A3[i, i - 2, setting] = k2
            end
            if i + 2 <= n
                A3[i, i + 2, setting] = k2
            end
        end
    end

    return A3
end

function solveDialYieldIP(H, K, CHD)
    n = length(H)
    s = size(K, 1)
    myModel = Model(Cbc.Optimizer)
    set_optimizer_attribute(myModel, "seconds", 300.0)

    A3 = constructA3D(H, K)
    @assert size(A3) == (n, n, s)

    # y[i,setting] = 1 if setting is used at location i
    @variable(myModel, y[1:n, 1:s], Bin)          
    @variable(myModel, R[1:n] >= 0)          
>>>>>>> origin/main
    @variable(myModel, z[1:n] >= 0)  # absolute value calculation

    # maximum smoothness
    @objective(myModel, Min, sum(z[i] for i in 1:n))

<<<<<<< HEAD
    # Depth definitiion
    @constraint(myModel, [i=1:n], R[i] >= H[i] + CHD)

    @constraint(myModel, [i=1:n], sum(x[i,k] for k in 1:4) <= 1)

    # Depth accumulation definition
    @constraint(myModel, [i=1:n], R[i] == sum(sum(A[i,j,k] * x[j,k] for j in 1:n) for k in 1:4))
=======
    # One setting per location
    @constraint(myModel, [i=1:n], sum(y[i,setting] for setting in 1:s) <= 1)

    # No adjacent bomb
    @constraint(myModel, [i=1:n-1],
        sum(y[i,setting] for setting in 1:s) + sum(y[i+1,setting] for setting in 1:s) <= 1
    )

    # Depth definitiion
    @constraint(myModel, [i=1:n], R[i] >= H[i] + CHD)

    # Depth accumulation definition
    @constraint(myModel, [i=1:n],
        R[i] == sum(A3[i,j,setting] * y[j,setting] for j in 1:n, setting in 1:s)
    )
>>>>>>> origin/main

    # absolute value constraint z[i] = | Ri-Hi-CHD |
    @constraint(myModel, [i=1:n], z[i] >= R[i] - (H[i] + CHD))

<<<<<<< HEAD
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
=======
    optimize!(myModel)

    if termination_status(myModel) == MOI.OPTIMAL || termination_status(myModel) == MOI.TIME_LIMIT
        yval = JuMP.value.(y)
        Rval = JuMP.value.(R)
        chosen = findall(v -> v > 0.5, yval)
        bomb_idx = [idx[1] for idx in chosen]
        bomb_setting = [idx[2] for idx in chosen]
        bomb_dist = 250 .* (bomb_idx .- 1)

        writedlm("data/dial_yield_positions.txt", bomb_dist)
        writedlm("data/dial_yield_settings.txt", bomb_setting)
        writedlm("data/dial_yield_R_values.txt", Rval)
>>>>>>> origin/main

        println("Total nukes used: $(length(bomb_idx)) Nukes")
    else
        println("Optimize was not succesful. Return code: ", termination_status(myModel))
    end
end

<<<<<<< HEAD
solveSmoothIP(H, CHD)
=======
solveDialYieldIP(H_part, K, CHD)
>>>>>>> origin/main
