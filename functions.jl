using ITensors
import ITensors: op


function expτ(s1::Index, s2::Index,w,τ)
        h = w*op("Sz", s1) * op("Sz", s2)
  return exp(τ * h)
end

function expτh(s::Index,h,τ)
        hh = h*op("Sz", s)
  return exp(τ * hh)
end

function corr(rho_in::MPO,i::Int,j::Int)
    N = length(rho_in)
    s = [siteinds(rho_in)[i][2] for i in 1:N]
    if i != j
        O = op("Sz",s,i)*op("Sz",s,j)
    else
        O = 0.25*op("I",s,i)
    end
    res = tr(apply([O],rho_in))
    return res
end


function mag(rho_in::MPO,i::Int)
    N = length(rho_in)
    s = [siteinds(rho_in)[i][2] for i in 1:N]
    O = op("Sz",s,i)
    res = tr(apply([O],rho_in))
    return res
end

function corr_mat(rho)
    N = length(rho)
    c_mat = zeros(N,N)    
    for i in 1:N
        for j in 1:N
            c_mat[i,j] = corr(rho,i,j)
        end
    end
    
    return c_mat
end

function mag_mat(rho)
    N = length(rho)
    mag_v = [mag(rho,i) for i in 1:N]
    return mag_v
end
    
function overlap(state,rho)
    N = length(rho)
    s = [siteinds(rho)[i][2] for i in 1:N]
    psi_d = MPS(s,state)
    rho_d = outer(psi_d',psi_d)
    return inner(rho,rho_d)
end

function thermal_state(N,w,θ)
    cutoff = 1E-10
    δτ = 0.01 
    beta_max = 1

    s = siteinds("S=1/2", N)

    gates = ITensor[]

    for i in 1:N
        for j in i+1:N
            e = expτ(s[i],s[j],w[i,j],-δτ)
            push!(gates,e)
        end
    end

    for i in 1:N
        eh = expτh(s[i],θ[i],-δτ)
        push!(gates,eh)
    end

    rho = MPO(s, "Id") ./ √2

    terms = OpSum()

    for i in 1:N
        for j in i+1:N
            terms += w[i,j],"Sz", i, "Sz", j
        end
    end

    for i in 1:N
        terms += θ[i],"Sz",i
    end

    H = MPO(terms, s)
    
    for β in 0:δτ:beta_max
        rho = apply(gates,rho;cutoff)
        rho = rho / tr(rho)
    end
    
    return rho
end
