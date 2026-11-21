
using AirfoilSimulation

p = read_setting("setting.dat")
s1 = allocate_sim(p)
calbeta!(s1)

# Check fs0!
fs0!(s1)
println("Step 0 u max: ", maximum(abs.(s1.U)))
println("Step 0 ut max: ", maximum(abs.(s1.UT)))
