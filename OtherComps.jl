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