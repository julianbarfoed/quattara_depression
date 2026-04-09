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
    @variable(myModel, z[1:n] >= 0)  # absolute value calculation

    # maximum smoothness
    @objective(myModel, Min, sum(z[i] for i in 1:n))

    # Depth definitiion
    @constraint(myModel, [i=1:n], R[i] >= H[i] + CHD)

    @constraint(myModel, [i=1:n], sum(x[i,k] for k in 1:4) <= 1)

    # Depth accumulation definition
    @constraint(myModel, [i=1:n], R[i] == sum(sum(A[i,j,k] * x[j,k] for j in 1:n) for k in 1:4))

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

solveSmoothIP(H, CHD)