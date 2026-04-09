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
    @variable(myModel, z[1:n] >= 0)  # absolute value calculation

    # maximum smoothness
    @objective(myModel, Min, sum(z[i] for i in 1:n))

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

    # absolute value constraint z[i] = | Ri-Hi-CHD |
    @constraint(myModel, [i=1:n], z[i] >= R[i] - (H[i] + CHD))

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

        println("Total nukes used: $(length(bomb_idx)) Nukes")
    else
        println("Optimize was not succesful. Return code: ", termination_status(myModel))
    end
end

solveDialYieldIP(H_part, K, CHD)
