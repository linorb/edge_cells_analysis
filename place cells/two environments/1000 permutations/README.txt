The structure of <mouse name>_place_cells.pkl:
min_number_of_events: the threshold of number of complex events for calculating 
    the place cells analysis. <int>
number_of_permutations: the number of permutations for the significance test of 
    the spatial information <int>
velocity_threshold: the absoulte velocity threshold for calculating 
    the place cells analysis. <int>
sessions_place_cells: a list of the place cells from all sessions:
    the first branch is the sessions
    the second branch is the session type [0]: environment a
                                          [1]: environment b
    the thirsd bracnch is the place cells function output: 
    [0]: place cells global indices (according to the cell registration "cellRegisteredFixed")
    [1]: place cells significance
    [2]: place cells information content
    for example: 
    sessions_place_cells[2][0][0] are the place cells of day 3 for environment a
    sessions_place_cells[2][1][0] are the place cells of day 3 for environment b
    sessions_place_cells[4][1][2] are the place cells' information of day 5 (4+1) for environment b
                 
    