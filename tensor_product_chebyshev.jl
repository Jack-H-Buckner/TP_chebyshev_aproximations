module TP_chebyshev


###################################################################
###                                                             ###
### define the basic machinery needed for chebyshev polynomials ### 
###                                                             ###
###################################################################



# the nth degree chebyshev polynomial
# evaluated at x
function T(x,n)
    return cos(n*acos(x))
end 

    
# the product of chebyshev polynomials
function T_alpha_i(alpha_i,x)
    return prod(T.(x,alpha_i)) 
end
    
# the product of chebyshev polynomials
# summed for each value of alpha
function T_alpha(x,alpha, coefs)
    T_i = broadcast(a -> T_alpha_i(a,x), alpha)
    return sum(T_i .* coefs)
end
    
    
    
# collect terms for tensor products
function collect_alpha(m,d)
    if d == 1
        alpha = 0:m
    elseif d == 2
        alpha = Iterators.product(0:m,0:m)  
    elseif d == 3
        alpha = Iterators.product(0:m,0:m,0:m)
    elseif d == 4
        alpha = Iterators.product(0:m,0:m,0:m,0:m)
    elseif d == 5
        alpha = Iterators.product(0:m,0:m,0:m,0:m,0:m)
    elseif d == 6
        alpha = Iterators.product(0:m,0:m,0:m,0:m,0:m,0:m)
    elseif d == 7
        alpha = Iterators.product(0:m,0:m,0:m,0:m,0:m,0:m,0:m)
    elseif d == 8
        alpha = Iterators.product(0:m,0:m,0:m,0:m,0:m,0:m,0:m,0:m)
    end
    alpha = collect(alpha)[sum.(alpha) .<= m]
    return alpha
end
    

# collect z grid 
# creates and array of that stores
# the nodes for each dimension 
function collect_z(m,d)
    f = n -> -cos((2*n-1)*pi/(2*m))
    z = f.(1:m)
    if d == 1
        z_array = z
    elseif d == 2
        z_array = Iterators.product(z,z)  
    elseif d == 3
        z_array = Iterators.product(z,z,z) 
    elseif d == 4
        z_array = Iterators.product(z,z,z,z)
    elseif d == 5
        z_array = Iterators.product(z,z,z,z,z) 
    elseif d == 6
        z_array = Iterators.product(z,z,z,z,z,z)
    elseif d == 7
        z_array = Iterators.product(z,z,z,z,z,z,z)
    elseif d == 8
        z_array = Iterators.product(z,z,z,z,z,z,z) 
    end
    return collect(z_array)
end

function compute_coefs(d,m,values,alpha)
    # get size of 
    @assert length(values) == m^d
    @assert all(size(values) .== m)
    # get constants for each term
    d_bar = broadcast(alpha_i -> sum(alpha_i .> 0), alpha)
    c = 2 .^ d_bar ./(m^d)
    # compute sum for each term
    z = collect_z(m,d) 
    # this is complicated. it initally broadcasts T_alpha_i over the grid z and takes the dot 
    # product with values array. It then broad casts this function over each set of valeus in alpha
    vT_sums = broadcast(alpha_i -> sum(broadcast(x -> T_alpha_i(alpha_i,x), z).*values), alpha)   
    return c.*vT_sums
end 

    
    
################################################################
###                                                          ###
### define a data structure and assocaited methods to store, ###
### initlaize, updated and evaulate the interpolations.      ###
###                                                          ###
################################################################
    

    
# define a data structure to save 
# informaiton for interpolation
mutable struct interpolation
    d::Int64 # dimensions
    a::AbstractVector{Float64} # lower bounds
    b::AbstractVector{Float64} # upper bounds
    m::Int # nodes per dimension
    nodes::AbstractMatrix{Float64} # m by d matrix with nodes in each dimension
    values # value associated with each node (m^d entries)
    coeficents::AbstractVector{Float64} # coeficents for computing the polynomial
    alpha::AbstractVector{Any}  # vector of tuples with m integer arguments 
end  

# define a function to initialize the 
# interpolation data structure with zeros
function init_interpolation(a,b,m)
    @assert length(a) == length(b)
    d = length(a)
    # calcualte nodes
    f = n -> -cos((2*n-1)*pi/(2*m))
    z = f.(1:m)
    nodes = (z.+1).*transpose((b.-a)./2 .- a)
    # initialize values as zero
    values = zeros(ntuple(x -> m, d))
    coefs = zeros(binomial(m+d,m))
    alpha = collect_alpha(m,d)
    return interpolation(d,a,b,m,nodes,values,coefs,alpha)
end 


# define a function to update the values and 
# coeficents for the interpolation
function update_interpolation(current, new_values)
    # check size of new values
    @assert length(new_values) == current.m^current.d
    @assert all(size(new_values) .== current.m)
    # update values
    current.values = new_values
    # compute coeficents
    m = current.m
    d = current.d
    alpha = current.alpha
    coefs = compute_coefs(d,m,new_values,alpha)
    current.coeficents = coefs
    return current
end 

    
# evaluluates the interpolation at all points inclueded x
function evaluate_interpolation(x,interpolation)
    a = interpolation.a
    b = interpolation.b   
    z = broadcast(x -> (x.-a).*2 ./(b.-a) .- 1,x)
    v = broadcast(x -> T_alpha(x,interpolation.alpha, interpolation.coeficents),z)
    return v
end 


end #modle