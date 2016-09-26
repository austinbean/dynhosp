# ParallelFunctions.jl
# Stolen from this Stack Overflow Comment: http://stackoverflow.com/questions/27677399/julia-how-to-copy-data-to-another-processor-in-julia

function sendto(p::Int; args...)
    for (nm, val) in args
        @spawnat(p, eval(Main, Expr(:(=), nm, val)))
    end
end


function sendto(ps::Vector{Int}; args...)
    for p in ps
        sendto(p; args...)
    end
end

getfrom(p::Int, nm::Symbol; mod=Main) = fetch(@spawnat(p, getfield(mod, nm)))


function passobj(src::Int, target::Vector{Int}, nm::Symbol;
                 from_mod=Main, to_mod=Main)
    r = RemoteRef(src)
    @spawnat(src, put!(r, getfield(from_mod, nm)))
    for to in target
        @spawnat(to, eval(to_mod, Expr(:(=), nm, fetch(r))))
    end
    nothing
end


function passobj(src::Int, target::Int, nm::Symbol; from_mod=Main, to_mod=Main)
    passobj(src, [target], nm; from_mod=from_mod, to_mod=to_mod)
end


function passobj(src::Int, target, nms::Vector{Symbol};
                 from_mod=Main, to_mod=Main)
    for nm in nms
        passobj(src, target, nm; from_mod=from_mod, to_mod=to_mod)
    end
end
