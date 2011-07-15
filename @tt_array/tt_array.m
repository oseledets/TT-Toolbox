function t = tt_array(varargin)
% Constructs tt_array from different objects

if (nargin == 0)
    t1=tt_tensor;
    t.tt=t1;
    t.ns = zeros(0,0);
%     t.n=0; %n-dimensions
%     t.m=0; %m-dimensions
%     t.k=0; %k-dimensions
    t = class(t, 'tt_array');
    return;
end

% Copy CONSTRUCTOR
if (nargin == 1) && isa(varargin{1}, 'tt_array')
    t=tt_array;
    t.tt = varargin{1}.tt;
    t.ns = varargin{1}.ns;
%     t.m = varargin{1}.m;
%     t.k = varargin{1}.k;
    return;
end

% From tt_tensor
if ( nargin == 2 && isa(varargin{1},'tt_tensor') && isa(varargin{2},'double'))
    tt=varargin{1};
    ns=varargin{2};
%     m=varargin{3};
%     k=varargin{4};
    d=tt.d;
    if ( size(ns,1) == 1 )
      ns=ones(d,1)*ns;
    end
%     if ( numel(m) == 1 )
%       m=m*ones(d,1);
%     end
%     if ( numel(k) == 1 )
%       k=k*ones(d,1);
%     end    
    t=tt_array;
    t.ns=ns;
%     t.m=m;
%     t.k=k;
    t.tt=tt;
return
end

end