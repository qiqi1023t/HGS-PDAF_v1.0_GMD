- in grok file:

         GRID DEFINITION:

         read 3d grid, ascii
         end grid definition

         GRID DEFINITION II(OPTIONAL):

         control volume


         AT THE END OF THE GROK FILE
         THE K VALUES FOR THE ELEMENTS
         SHOULD BE READ IN:

         Read elemental k from file
         elementalk.txt
         
- in `<name>.coordinates`  ! Coordinates of the nodes
         
         <number of nodes>  ! In the first line the number of nodes is written.
         x_1 y_1 z_1        ! In the follwoing lines all the coordinates of the nodes are written.
         ...                ! First, the x-coordinates are written, then the y and z coordinates
         x_nrows y_nrows z_nrows
         0 0 2 0            ! 1st: max number of grid lines in x-direction, 2nd: max number of grid lines in y direction. (1st and 2nd are neglected for unstructured grids, which is the case here, therefore, both numbers can also be set to 0) 3rd: min number in z-direction (min=2 or higher), 4th: number of species for transport (sic!)
         65                 ! Number of elements or triangels in 2-D
         .false.            ! .true == tetrahedral mesh
         .false.            ! .true == Galerkin method is used for tetramesh
         .false.            ! .true == do_write_face_seg (?), HGS-Manual says: "should always be .false"
         9                  ! max. number of nodes connected to each node for a 2-D triangular mesh, "9" is a number given by Aquanty. However, "7" should be fine in 2-D and "8" in 3-D
         
- in `<name>.elements`

         <n_nodes_per_element>b ! 1st line number of nodes per element
         <n_elements>           ! 2nd line number of elements
         node_id_1 node_id_2 node_id_3 node_id_4 node_id_5 node_id_6 ! All node ids of the element. id_1, 2 and 3 correspond to the lower base ares (z-direction) and id_4, 5, 6 to the higher base area
         ...
         element_id_1 ! After all the node ids, in the first column all elemental ids follow directly, i.e. v = 1, ..., n , with n = number of elements
         ...
         
- in `elementalk.txt`

         id_1 1 1.0e-04 1.0e-04 1.0e-04 ! First column: elemental ids (integer), than K in x-direction, K in y-direction, K in z direction for all elements, i.e. number of rows = number of elements
         ...
         id_n n ...


However, all this specifications will give an error, when grok is running. To overcome the error, in the newly generated array_sizes.default file the number of "flow: material zones  " has to equal at least 
the number of elements. In the example case:

flow: material zones                                        
          65


Now everything should run fine!

