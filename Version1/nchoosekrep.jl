# Defines n choose k with repetition
# Generates all of the permutations.  For my purposes: n hospitals can be allocated among
# No Intensive or Intermediate, Intermediate and Intensive in nchoosekrep(n,k) different ways
# which are returned by the function. Then nchoosekrep(5,3) will return (5,0,0), (4,1,0), ...
# but it won't return them in any particular order.




#=
Right now the function below is basically non-functional above n = 13, which takes something like 12
GB of memory to compute.  Should rewrite to make recursive somehow, like the Matlab version, then
it will hopefully take up less memory.
=#
#=
function nchoosekrep(n::Int64,k::Int64)
	# This function takes two inputs (n, k) and returns the number of different combinations of 
	# values which can be drawn on [0:n].  So (5, 3) should return (5, 0, 0), (1, 2, 2),...
	# For the purpose of my paper I only need k = 3 case
	if n == 0
		return zeros(Int, 1, k)
	else 
		if k == 0
			return TypeError
			return println("You have to choose a positive number")	
		elseif k == 1
			return transpose( collect(0:n))'
		else
			a = collect(partitions(vcat(zeros(Int64, k), ones(Int64,n)), k))
			b = unique([[sum(y) for y in x] for x in a])
			H = zeros(Int, 1,k)
			for el in b
				container = Array(Int, maximum(size(unique(permutations(el)))) , k)
				for x in 1:length(unique(permutations(el)))
					container[x, :] = unique(permutations(el))[x]
					container
				end
				H = vcat(H, container)
			end
			H = unique(H,1)
			H = H[2:end, :]
			return H
		end
	end 
end

=#


