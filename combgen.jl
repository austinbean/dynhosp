#= Combinations generator - takes arrays (n, 3) and generates
   all of the combinations of rows taken separately from 
   each of the arrays.  So (r11, r21, r31), (r11, r21, r32)...
   where rij has i - array and j row of array i. 

   There is surely a smarter, recursive way to do this for any 
   arbitrary collection of arrays
   =#

function combgen3(A1::Array, A2::Array, A3::Array)
	lenA1 = size(A1)[1]
	lenA2 = size(A2)[1]
	lenA3 = size(A3)[1]
	# all of my applications will have three types, but just in case
	width = size(A1)[2] + size(A2)[2] + size(A3)[2]
	container = zeros(Int, 1, width)
	for i = 1:lenA1
		for j = 1:lenA2
			for k = 1:lenA3
				container = [container; [A1[i,:] A2[j,:] A3[k,:]]]
			end
		end
	end
	container = container[2:end, :]
	return container
end

function combgen2(A1::Array, A2::Array)
	lenA1 = size(A1)[1]
	lenA2 = size(A2)[1]
	width = size(A1)[2] + size(A2)[2]
	container = zeros(Int, 1, width)
	for i = 1:lenA1
		for j = 1:lenA2
			container = [container; [A1[i,:] A2[j,:]]]
		end
	end
	container = container[2:end, :]
	return container
end






