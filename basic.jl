using CSV, Tables
using LinearAlgebra, SparseArrays 

# read data
w1 = CSV.read("./data/w1.csv", Tables.matrix)
w2 = CSV.read("./data/w2.csv", Tables.matrix)
x = CSV.read("./data/x.csv", Tables.matrix)

beta = CSV.read("./data/beta.csv", Tables.matrix)[:,1]
nu = CSV.read("./data/nu.csv", Tables.matrix)[:,1]
e = CSV.read("./data/e.csv", Tables.matrix)[:,1]
y = CSV.read("./data/y.csv", Tables.matrix)[:,1]

n = length(y)
nx = size(x,2)

struct Parameters
    beta::Vector{Float64}
    rho1::Float64
    rho2::Float64
    sigmasq::Float64
    dev::Float64
end