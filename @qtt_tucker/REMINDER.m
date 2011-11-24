%Two simple functions are required
%The storage for tt_tensor is the following:
%Dimension d
%size vector n
%rank vector r
%cores, stored contigiously in one "long" array core
%with enumeration a_{k-1} i_k a_k
%The main function is to convert "ordinary" TT-cell representation 
%to such kind of representation and back
%from_cell(cell) -- realized as constructor
%to_cell(tt_tensor)