! ------- Well1 ---------

Choose node
175,250,85

create node set
Well1_node

boundary condition
    type
     flux nodal
    name
     Well1
    node set
     Well1_node
    time value table
     0.0 -0
    end 
     tecplot output
end

clear chosen zones
clear chosen nodes
clear chosen elements
clear chosen segments
clear chosen faces
!end

! ------- Well2 ---------

Choose node
175,150,85

create node set
Well2_node

boundary condition
    type
     flux nodal
    name
     Well2
    node set
     Well2_node
    time value table
     0.0 -0
    end 
     tecplot output
end

clear chosen zones
clear chosen nodes
clear chosen elements
clear chosen segments
clear chosen faces
!end
