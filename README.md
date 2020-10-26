# Shearsort with MPI
Program sorts cuboids by their volume from a given list of cuboids with id, length, width and height using shearsort with odd-even sort\
and writes output to console and to a text file.\
Parallel programing with MPI creating a new communicator using a cartesian topology - \
each process gets data of a different cuboid, calculates its volume and compares with neighbor.

# Shearsort  
For nxn matrix sorting requires sqrt(n)*(ceil(log(n)+1)) steps.\
Snake sort - even rows sorted in increasing order, odd rows sorted in decreasing order. Columns sorted in increasing order always.\
Rows and columns sorted with odd-even sort.

# Odd-Even sort 
For nxn matrix requires n steps for each row/column (n^2 for worst case) and operates in two alternating phases - \
Even phases - even row/column exchange value with right/bottom neighbor (unless last row/column, doesn't have a right/bottom neighbor).\
Odd phases - even row/column exchange value with left/upper neighbor (unless first row/column, doesn't have a left/upper neighbor).
