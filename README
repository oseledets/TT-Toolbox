TT-Toolbox (TT=Tensor Train) Version 2.2.2

TT(Tensor Train) format is an efficient way for low-parametric
representation of high-dimensional tensors. The TT-Toolbox
is a MATLAB implementation of basic operations with 
tensors in TT-format. It includes: 

    * tt_tensor and tt_matrix classes for storing vectors and operators
    * Basic linear algebra subroutines (addition, matrix-by-vector product, 
    elementwise multiplication and many others) using standard MATLAB syntax, 
    linear complexity in the dimension, reshape function
    * Fast rounding procedure with a prescribed accuracy
    * Advanced approximation and solution techniques:
        * Approximate solution of linear systems and eigenvalue problems 
        * Cross methods to approximate ``black-box'' tensors
        * Wavelet tensor train decomposition
    * Construction of basic operators and functions (Laplace operator, function of a TT-tensor)
    * Computation of maximal and minimal elements of a tensor
    * and several others
	
New in Version 2.2.2
    * one-element MATLAB-style subsref and subsassgn now possible for tt_matrix:
	current_value=ttm(123,4567) and ttm(123,4567)=new_value
    * submatrix extraction is now possible: 
	ttm2=ttm([1 2 7],[1 2 1;1 1 2]) 
	For additional help see code of @tt_matrix/subsref.m and core/tt_submatrix.m
	
New in Version 2.2.1
    * Blockwise assembly of tt_matrices added:
	ttM=tt_blockwise({ttA ttB ttC;ttD ttE ttF})
    * subtensor extraction is now possible: 
	ttv2=ttv(1,2,:,2,2,:,:,:,1) or ttv2=ttv({1,2,':',2,2,':',':',':',1})
	For additional help see code of @tt_tensor/subsref.m and core/tt_subtensor.m

New in Version 2.2
    * Better documentation
    * Mixed QTT-Tucker format (qtt_tucker class)
    * reshape function for a TT-tensor/TT-matrix
    * dmrg_cross method for black-box tensor approximation
    * Convolution in QTT-format

 TT-Toolbox 2.2. can be downloaded from http://spring.inm.ras.ru/osel
