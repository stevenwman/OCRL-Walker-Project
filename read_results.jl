
# open file called "results.jld2"

# read the file and store the data in a variable called data
data = JLD2.load("results.jld2")

# print the keys of the data
# println(keys(data))

# Print value at key "X"
# println(data["X"])

Xm = data["Xm"]
Um = data["Um"]
t_vec = data["params"].t_vec

display(plot(Xm[1,:],Xm[2,:], label = "body"))
display(plot(t_vec[1:end-1], Um',
             xlabel = "time (s)", ylabel = "U",
             label = ["u1" "u2" "u3" "u4"], title = "Controls"))