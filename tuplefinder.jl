#=
This function takes an array and a set of three columns as arguments and returns
all of the unique 3 tuples which can be made out of elements appearing in the
columns.  (Not rearranging the elements.)

tuplefinder( X, 1, 2, 3) applied to
[a b c]
[d e f]
[a b c]
will return (a,b,c) and (d,e,f) but NOT (a,e,f).
=#


function tuplefinder(combos::Array, col1::Int64, col2::Int64, col3::Int64)
  emplist = [("a", "b", "c")]
  for i = 1:size(combos)[1]
    emplist = vcat(emplist, (combos[i,col1], combos[i,col2], combos[i,col3]))
  end
  emplist = emplist[2:end, :]
  return unique(emplist)
end
