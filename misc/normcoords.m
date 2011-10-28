function [Z,x]=normcoords(u,varargin)
% function [Z,x]=normcoords(u,varargin)
% Returns an orthogonal matrix Z of normal coordinates axes in columns, and
% (optionally) the vector of extreme indices x.
% u must be the TT 1.0 vector, varargins are:
%   'x' - specifies the vector of extremal indices,
%   'h' - mesh step size

x = [];
h=1;

for i=1:2:length(varargin)-1
    if (~isempty(varargin{i+1}))
        switch lower(varargin{i})
            case 'x'
                x=varargin{i+1};
            case 'h'
                h=varargin{i+1};
            otherwise
                error('Unrecognized option: %s\n',varargin{i});
        end;
    end;
end;

if (isempty(x))
    [maxel,x]=tt_max(u);
end;

Hess=tt_hess(u, x, h);
maxel = mean(abs(diag(Hess)));
% Hess(abs(Hess)./maxel<1e-6)=0;

[Z,L]=eig(Hess);

end