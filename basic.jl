using CSV, Tables
using LinearAlgebra, SparseArrays 
using Optim, ForwardDiff

# notes:
# 1. it seems like In the R code we can compute Xb once for muhat and n2ll.rho

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

mutable struct Parameters
    beta::Vector{Float64}
    rho1::Float64
    rho2::Float64
    sigmasq::Float64
    dev::Float64
end

Parameters(nx) = Parameters(zeros(nx), 0, 0, 1, Inf)

parm = Parameters(nx)
W1 = deepcopy(w1)
W2 = deepcopy(w2)

# first iteration of estimate

parm_new = deepcopy(parm)

# aggregate weight matrices
W1a = I(n) - W1*parm_new.rho1
W2a = I(n) - W2*parm_new.rho2

# assume covariates given, estimate beta | rho
tXtW2aW2a = transpose(x) * transpose(W2a) * W2a
parm_new.beta .= (tXtW2aW2a * x) \ (tXtW2aW2a * W1a * y)

#Estimate sigma | beta, rho
parm_new.sigmasq = sigmasqhat(muhat(y,x,W1a,W2a,parm_new.beta))

# (if not final iteration) estimate rho | beta, sigma
n2ll_rho = make_n2ll_rho(n,parm_new.beta,parm_new.sigmasq,W1,W2,y,x)
# n2ll_rho([parm_new.rho1;parm_new.rho2])

rho = [parm_new.rho1,parm_new.rho2]
temp = optimize(n2ll_rho, rho, BFGS(), autodiff=:forward)


# fn to make n2ll.rho
function make_n2ll_rho(n,beta,sigmasq,W1,W2,y,x)
    function n2ll_rho(rho)
        # rho1 (AR)
        W1a = I(n) - W1*rho[1]
        W1ay = W1a * y
        adetW1a = abs(det(W1a))

        # rho2 (MA)
        W2a = I(n) - W2*rho[2]
        tpW2a = transpose(W2a) * W2a
        adetW2a = abs(det(W2a))

        # nx
        Xb = x * beta

        # negloglike
        return n*(log(2*π)+log(sigmasq)) + transpose(W1ay-Xb)*tpW2a*(W1ay-Xb)/sigmasq - 2*(log(adetW1a)+log(adetW2a))
    end
end

# Estimate predicted means, conditional on other effects
function muhat(y,X,W1a,W2a,betahat)
    Xb = X * betahat
    return W2a * (W1a * y - Xb)
end

#Estimate innovation variance, conditional on other effects
sigmasqhat = function(muhat)
    transpose(muhat) * muhat / prod(size(muhat))
end