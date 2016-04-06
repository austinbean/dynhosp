#= Does not really do n choose k with repetions, but
does enumerate the ways N hospitals can be allocated across
3 categories.=#

function nchoosekrep(v::Int64)
	y = zeros(Int64, 1, 3)
		for i = 0:v
			for j = 0:v
				for k = 0:v
					if i + j + k <= v # this allows exit (less than v)
						y = [y ; [ i j k]]
					end
				end
			end
		end
	return y[2:end, :]
end

function nckr(v::Int64)
	y = zeros(Int64, 1, 3)
		for i = 0:v
			for j = 0:v
				for k = 0:v
					if i + j + k == v # this forbids exit
						y = [y ; [ i j k]]
					end
				end
			end
		end
	return y[2:end, :]
end


#=
What do I want?
- Track exit: when (ii), (ij), (ik) < v[], let me know and by how much
- Permit Entry: finite number of patterns

=#
function nckrexen(v::Array{Int64,2})
	#= The form of the data below is 2 sets of 3 integers followed by 3 sets of 4
	   integers.  Then an additional set of three integers.
		 - The first three are the number of Lev1, Lev2, Lev3 in the current
		 period.
		 - The second are the number of Lev1, Lev2, Lev3 in the next period.
		 - The three sets of four integers to follow are choices of those lev1, lev2,
		 or lev3 among: (ii) (ij) (ik) (iEX) transitions.
		 - The last set of three covers entry possibilities, allowing up to 2.
		 =#
	rows1 = nexit(v[1])
	rows2 = nexit(v[2])
	rows3 = nexit(v[3])
	y = zeros(Int64, 1, 22) # change to 22 to allow a column for "no entry"
	for r1 in 1:size(rows1)[1]
		for r2 in 1:size(rows2)[1]
			for r3 in 1:size(rows3)[1]
				if ( rows1[r1, 2] + rows1[r1, 3] + rows2[r2, 1] + rows2[r2, 3] + rows3[r3, 1] + rows3[r3, 2]  <= 2) # here we permit only two to change
					for x in 1:size(nchoosekrep(2),1)
						entrants = slice(nchoosekrep(2),x,:)
						if sum(entrants) == 0
							entrants = vcat(1, entrants)
						else
							entrants = vcat(0, entrants)
						end
						a = [(rows1[r1,1] + rows2[r2,1] + rows3[r3,1]) (rows1[r1,2] + rows2[r2,2] + rows3[r3,2]) (rows1[r1,3] + rows2[r2,3] + rows3[r3,3]) ]
						y = vcat(y, hcat(v, [a[1] + entrants[2] a[2] + entrants[3] a[3] + entrants[4]], rows1[r1,:], rows2[r2,:], rows3[r3,:], entrants' ) )
					end
				end
			end
		end
	end
	y = y[2:end, :]
	return y
end


# (v - i - j - k <= v - 2   )

function nexit(v::Int64)
	y = zeros(Int64, 1, 4)
		for i = 0:v
			for j = 0:v
				for k = 0:v
					if (i + j + k <= v) & (max(v- 2,0) <= i+j+k) # this allows exit (less than v) and up to two exits
						y = [y ; [ i j k (v - i - j - k)]]
					end
				end
			end
		end
	return y[2:end, :]
end



#=
# This one attempts to do it recusively, so far without success
# another idea - append everything and then reshape it all?

function nckr(v::Array{Int64,1}, n::Int64)
	if n == 1
		y = v
	else
		println(n)
		m = length(v)
		y = zeros(Int64, m)
		if m == 1
			y = v
		else
			for i = 1:m
				y_recr = nckr(v[i:length(v)], n - 1)
				s_repl = fill(v[i], size(y_recr, 1))
				println(size(y), size(s_repl), size(y_recr))
				y = hcat(y , hcat(s_repl, y_recr) )
			end
		end
	end
	return y
end
=#
