#Other paper computations: 

#=

Miscellaneous model computations for the paper:

=#



"""
`WTPchange(d::DynState)`
Take each firm.  Update the level.  Update the WTP.  Compute the change.

dyn = CounterObjects(50);


"""
function WTPchange(d::DynState)
    outp::Array{Float64,2} = zeros(286, 3)
    for ix in 1:size(dyn.all,1)
        outp[ix,1] = dyn.all[ix].fid
        if dyn.all[ix].level == 1
            println(dyn.all[ix].fid)
            outp[ix,2] = FindWTP(dyn.all[ix])
            println(FindWTP(dyn.all[ix]))
            dyn.all[ix].level = 3
            dyn.all[ix].tbu = true 
            # Ok - this doesn't update yet.  
            UpdateDUtil(dyn.all[ix]) 
            outp[ix,3] = FindWTP(dyn.all[ix])
            println(FindWTP(dyn.all[ix]))
            dyn.all[ix].level = 1
            dyn.all[ix].tbu = true 
            UpdateDUtil(dyn.all[ix])
            println("********")
        elseif dyn.all[ix].level == 2
            outp[ix,2] = FindWTP(dyn.all[ix])
            dyn.all[ix].level = 3 
            dyn.all[ix].tbu = true 
            UpdateDUtil(dyn.all[ix])
            outp[ix,3] = FindWTP(dyn.all[ix])
            dyn.all[ix].level = 2
            dyn.all[ix].tbu = true 
            UpdateDUtil(dyn.all[ix])
        elseif dyn.all[ix].level == 3
            # do nothing.
        end 
    end 
    return outp
end 

"""
`DUtilCheck`
See whether the utilities are updated correctly for "shortrecs"
"""
function DUtilCheck(h::simh)
    println("Facilities")
    for el in h.mk.m 
        println(el.zp, " ", size(unique(el.putils[1,:]))[1], "  ", size(unique(el.mutils[1,:]))[1])
    end 


    println("Own Utility Update Check: ")
    for el in h.mk.m 
        println(el.zp, "  ", el.putils[2, findin(el.putils[1,:], h.fid)])
    end 
    h.level = 3
    h.tbu = true 
    UpdateDUtil(h)
    println("********")
    for el in h.mk.m 
        println(el.zp, "  ", el.putils[2, findin(el.putils[1,:], h.fid)])
    end
    h.level = 1 
    h.tbu = true 
    UpdateDUtil(h) 
end 

#=
This computation relates to determining the market sizes of individual hositals by checking the total number of births in 
all zip codes attached to that hospital every year.  It is merged with some data from Stata.  Specifically, see TX Patient Uncertainty.do

=#
mx = 0;
for el in dyn.all 
    if size(el.mk.m,1) > mx 
        mx = size(el.mk.m,1)
    end 
end 
println("Max zips = ", mx)

zps = zeros(Int64, mx+1, size(dyn.all,1) ); # NB - each COLUMN is an fid 

for el in 1:size(dyn.all,1) # COLUMNS
    zps[1, el] = dyn.all[el].fid  
    for z in 1:size(dyn.all[el].mk.m,1) # ROWS 
        zps[z+1,el] = dyn.all[el].mk.m[z].zp # zip 
    end 
end 
zps = transpose(zps)
writecsv( "/Users/austinbean/Desktop/dynhosp/hospzips.csv", zps)




#=

typetest 


=#

struct WTPee
    w::Array{Float64, 1}
end

ab1 = WTPee(Array{Float64,1}())

ab2 = WTPee(Array{Float64,1}(50))

@time for el =1:50
    push!(ab1.w, rand())
end 

@time for el =1:50 
    ab2.w[el] = rand()
end 

for j = 1:10
@time for i = 1:50
 ab2.w[i] = i 
end 
for k = 1:50
    ab2.w[k] = 0.0
end 
end 


#=

=#