using CSV, Tables
using LinearAlgebra, SparseArrays 
using Optim, ForwardDiff, LineSearches

include("./lnam.jl")

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
rho = [0.2,0.1]
sigmasq = 0.01

nll = make_nll(n,nx,w1,w2,y,x)
n2ll_rho = make_n2ll_rho(n,nx,w1,w2,y,x,beta,sigmasq)

nll([beta...,rho...,sigmasq]) ≈ n2ll_rho(rho)/2