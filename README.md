# tensor_product_chebyshev
This repository provides a data structure and set of methods for aproximating multidimensional functions with chebyshev poynomials. These function aproximation methods should be very useful for dynamic programming and other applications that solve for arbitrary continuous functions. 

The library is based around a single data structure, a mutable struct call `interpolation`. It stores data that describe the number of dimensions of the function the domain of the function and the number of nodes used for the aproximation. The value of the function at the nodes and a vector of polynomial coeficents used to represent the function are also stored and can be used to evaluate the function at points that lie off of the grid. 

The `interpolation` data structure has three assocaited methods `init_interpolation`, `update_interpolaiton` and `evaluate_interpolation`.  `init_interpolation` is used to generate an instance of `interpolation`,  `update_interpolaiton` allows the user to change the value of the function at the grid points and `evaluate_interpolation` allows to user to retriev the function value any where in its domain. 

In addition to these primary methods there are several helper functions used to define the primary methods that can be accessed by the user. These functions are used to define the chebyshev basis functions and to calcualte the optimal locaitons for the nodes. 

## `interpolation(d,a,b,m,nodes,values,alpha,coefs)` mutable struct

d - an integer specifying the number of dimensions of the function domain

a - a vector of floating point numbers. Gives the lower bounds in each dimension 

b - a vector of floating point numbers such that a[i] < b[i] . Gives the upper bounds in each dimension 

m - an integer specifying the number of nodes in each dimension

nodes - a d by m matrix of floating point number specifying the locaiton of the nodes in each dimension

values - a m by m by ... m array of floating point numbers specifying the value of the function at each node

alpha - a m + d choose d vector of tuples of length d. used for the tensor product

coeficnets - a vector of length m + d choose d that parameterizes the chebyshev polynomial 


## `init_interpolation(a,b,m)` function

### arguments

a - vecor of float the give lower bound

b - vector of floats that give upper bounds

m - integer the give number of nodes 

### value

instance of `interpolation` with values set to zeros

## `update_interpolation(itp,values)` function

### arguments

itp - an instrnace of `interpolation`

values - an arry or new values for each node in itp

### value

an instance of `interpolation` with updated values and coeficents


## `evaluate_interpolation(x, itp)` function

### arguments

itp - an instrnace of `interpolation`

x - a vector of vectors or tuples 

### value

a vector of floating point numbers giving the value of the aproximated function for each value given in x.


