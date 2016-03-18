function z=mtkron(varargin)
    %"Tensor" Kronecker product (in "TT" order) of multiple matrices
    %   Z=MTKRON(A,B,C,...) Takes Kronecker product of the arguments
    %   See also TT_TENSOR/TKRON, TT_MATRIX/TKRON 
    %
    %   
    %   it also takes a cell array of tt_tensors as a single input
    %
    % TT-Toolbox 2.2.1, 2009-2016
    %
    %This is TT Toolbox, written by Ivan Oseledets et al.
    %Institute of Numerical Mathematics, Moscow, Russia
    %webpage: http://spring.inm.ras.ru/osel
    %
    %For all questions, bugs and suggestions please mail
    %ivan.oseledets@gmail.com
    %---------------------------
    %Last tweaks by Alexey Boyko of Oseledets' group
    %alexey.boyko@skolkovotech.ru

    z=varargin{1};
    if numel(z)~=1
        % in case it a cell array of tt_tensors as a single input
        varargin=z;
        z=varargin{1};
    end
    for i=2:numel(varargin)
      z=tkron(z,varargin{i});
    end
end
