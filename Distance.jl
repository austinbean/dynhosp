# Function defines the great circle distance between two points given
# (latitude, longitude) (w,x) and (y,z)


function distance(w::Float64, x::Float64,y::Float64,z::Float64)
  rad = 3961
  conv = pi/180
  w = w*conv
  x = x*conv
  y = y*conv
  z = z*conv
  d = 2*rad*asin(sqrt((sin((w-y)/2))^2 + cos(w)*cos(y)*(sin((x-z)/2))^2))
  return d
end
