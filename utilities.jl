

"""
`MktSize(d::DynState)`
How many patients of each type are there in each firm's market?
"""
function MktSize(d::DynState)
  emp::Array{Int64,1} = Array{Int64,1}()
  for el in d.all 
    mp::Int64 = 0
    pp::Int64 = 0
    for z in el.mk.m  
      pp += sum(z.pcounts)
      if sum(z.pcounts) == 0
        push!(emp,z.zp)
      end 
      mp += sum(z.mcounts)
      if sum(z.mcounts) == 0
        push!(emp, z.zp)
      end 
    end 
    println(el.fid, " private: ", pp, " medicaid: ", mp)
  end 
  println(unique(emp))
end 




"""
`MarketPrint(mkt::Market)`
Prints the elements of the market record: name, neighbors, choice probabilities.
For debugging purposes to make sure things look right.
"""
function MarketPrint(mkt::Market)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
  println(mkt.fipscode)
  for el in mkt.config
    println("*******************")
    println(el.name)
    println(el.neigh)
    println(el.hood)
    println(el.chprobability)
  end
end





"""
`NeighborsPrint(mkt::Market)`
Prints the name and neighbors of every hospital in the market `mkt`
"""
function NeighborsPrint(mkt::Market)
  for el in mkt.config
    println(el.name, " ", el.neigh)
  end
end




"""
`FacPrint(hosp::hospital)`
Prints out a hospital facility record: name, fid, fips, level, neighbors.
"""
function FacPrint(hosp::hospital)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
  println(hosp.name)
  println(hosp.fid)
  println("Fips: ", hosp.fipscode)
  println("Level: ", hosp.level)
  println("Neighbors: ", hosp.hood)
  println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
end




"""
`PrintZip(zi::zip)`
Prints the fid and the name of the facilities attached to the zips.
"""
function PrintZip(zi::zip)
  for el in keys(zi.facilities)
    println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
    println(el, "  ", zi.facilities[el].name)
  end
    println("⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒⭒")
    println("385 ", zi.ppatients.count385, " ", zi.mpatients.count385)
    println("386 ", zi.ppatients.count386, " ", zi.mpatients.count386)
    println("387 ", zi.ppatients.count387, " ", zi.mpatients.count387)
    println("388 ", zi.ppatients.count388, " ", zi.mpatients.count388)
    println("389 ", zi.ppatients.count389, " ", zi.mpatients.count389)
    println("390 ", zi.ppatients.count390, " ", zi.mpatients.count390)
    println("391 ", zi.ppatients.count391, " ", zi.mpatients.count391)
end



"""
`ZeroFind(mat::Array{Float64,2})`
looks for hospitals which were never demanded - there are some which have rows of zeros
will want to pop those out of the results matrix.
"""
function ZeroFind(mat::Array{Float64,2})
  zers = zeros(size(mat,1)) # vector of zeros for output
  for i = 1:size(mat,1)
    zercnt = 0
    for j = 1:size(mat, 2) # don't want to count all of them - fix this.
      if mat[i,j] == 0.0
        zercnt += 1
      end
    end
    zers[i] = zercnt
  end
  return zers
end


"""
`CheckPats(pats::patientcollection)`
This will compute the total sum of all patients in all zips in the patientcollection.
This is just to check that all are being added as expected.
"""
function CheckPats(pats::patientcollection)
  count::Int64 = 0
  for k in keys(pats.zips)
    count += sum(pats.zips[k].mpatients)
    count += sum(pats.zips[k].ppatients)
  end
  return count
end
