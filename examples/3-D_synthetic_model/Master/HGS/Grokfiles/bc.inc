clear chosen zones
clear chosen nodes
clear chosen elements
clear chosen faces
clear chosen segments


Choose nodes am
../Mesh/Selections/Upstream_BC
1,34


create face set
Upstream


boundary condition
  type
    fluid transfer
  name
    Upstream_Specified_head
  face set
    Upstream
  time value table
    0.0  99.47
  end


fluid transfer coefficients
 500
 1e3
end

clear chosen zones
clear chosen nodes
clear chosen elements
clear chosen faces
clear chosen segments


Choose nodes am
../Mesh/Selections/Downstream_BC
1,34


create face set
Downstream


boundary condition
  type
    fluid transfer
  name
    Downstream_Specified_head
  face set
    Downstream
  time value table
    0.0  93.22
  end


fluid transfer coefficients
 500
 1e3
end

