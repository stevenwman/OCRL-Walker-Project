include("hybrid_traj_opt_5link/DynamicSystems/biped5link.jl")
include("hybrid_traj_opt_5link/DynamicSystems/constraints.jl")

# open file called "results.jld2"

# read the file and store the data in a variable called data
data1 = JLD2.load("results.jld2")
data10 = JLD2.load("results_10xR.jld2")
data100 = JLD2.load("results_100xR.jld2")
data1000 = JLD2.load("results_1000xR.jld2")

data_vec = [data1, data10, data100, data1000]
Rx = [1, 10, 100, 1000]

plot()

for i in 1:4
    data = data_vec[i]
    Xm = data["Xm"]
    Um = data["Um"]

    Um = Um[1:4,:]

    totalInput = norm(vcat(Um...))

    legendText = string(Rx[i]) * "xR, Input = " * string(round(totalInput, digits=1))

    # t_vec = data["params"].t_vec

    plot!(Xm[1,:],Xm[2,:], label = legendText, xlabel = "Torso x-position (m)", ylabel = "Torso y-position (m)", legend=:outertopright)
    # display(plot(t_vec[1:end-1], Um',
    #              xlabel = "time (s)", ylabel = "U",
    #              label = ["u1" "u2" "u3" "u4"], title = "Controls"))
end

xs = -1:0.01:4
ys = height_stairs.(xs)
display(plot!(xs, ys, linestyle=:solid, label="Ground", color="black"))

png("torso_plot.png")


model = (m1 = .5,  m2 = .5,  m3 = 1,  m4 = .5,  m5 = .5,  m6 = .01,
            l12 = 1, l23 = 1, l34 = 1, l45 = 1, l36 = 1, g = 9.81)
            
X = data1000["X"]

animate_walker(X, model)



# print the keys of the data
# println(keys(data))

# Print value at key "X"
# println(data["X"])

# Xm = data["Xm"]
# Um = data["Um"]
# t_vec = data["params"].t_vec

# display(plot(Xm[1,:],Xm[2,:], label = "body"))
# display(plot(t_vec[1:end-1], Um',
#              xlabel = "time (s)", ylabel = "U",
#              label = ["u1" "u2" "u3" "u4"], title = "Controls"))