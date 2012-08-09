classdef tt_function
    %TT_FUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fun;
        coeff;
    end
    
    methods
        %-----------------------------------------------------------------%
        function obj = tt_function(fun, coeff)
            obj.fun = fun;
            obj.coeff = coeff;
        end
        
        %-----------------------------------------------------------------%
        function tt = to_tt_tensor(obj, x)
        %TO_TT_TENSOR Converts tt_function to tt_tensor
        %   x - d x 1 cell array of coordinates for evaluation
        
        d = obj.coeff.d;
        r = obj.coeff.r;
        n = obj.coeff.n;
        ps = obj.coeff.ps;
        cr = obj.coeff.core;
        n_tt = cellfun(@numel, x);
        ps_tt = cumsum([1 ; n_tt.*r(1:d).*r(2:d+1)]);
        cr_tt = zeros(ps_tt(d+1)-1, 1);
        
        % compute basis functions
        basis = tt_function.get_basis(x, obj.fun, n);
        
        % evaluate cores
        for i = 1: d
            c = reshape(cr(ps(i):ps(i+1)-1), [r(i) n(i) r(i+1)]);
            t = ten_conv(c, 2, basis{i});
            cr_tt(ps_tt(i):ps_tt(i+1)-1) = reshape(t, [r(i)*n_tt(i)*r(i+1) 1]);
        end
        
        tt = tt_tensor;
        tt.d = d;
        tt.r = r;
        tt.n = n_tt;
        tt.ps = ps_tt;
        tt.core = cr_tt;
        
        end
        
        %-----------------------------------------------------------------%
        function val = eval(obj, x)
        %EVAL Evaluates tt_function at a given d-dimensional point
        %   Inefficient code, simple solution
        
        d = numel(x);
        t = cell(d,1);
        for i = 1: d
            t{i} = x(i);
        end
        tt = obj.to_tt_tensor(t);
        val = tt(ones(1, d));
        
        end
        
    end
    
    methods (Static)
        %-----------------------------------------------------------------%
        function basis = get_basis(x, fun, n)
        %GET_BASIS Constructs a cell array of values of basis functions.

        d = numel(n);
        basis = cell(d, 1);
        
        for i = 1: d
            basis{i} = zeros(n(i), numel(x{i}));
            for j = 1: n(i)
                for k = 1: numel(x{i})
                    basis{i}(j,k) = fun(i, j, x{i}(k));
                end
            end
        end

        end
        
    end
    
end

