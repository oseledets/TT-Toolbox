function t = tt_matrix(varargin)
%TT_TENSOR Tensor stored as a TT-operator

if (nargin == 0)
    t1=tt_tensor;
    t.tt=t1;
    t.n=0; %n-dimensions
    t.m=0; %m-dimensions
    t = class(t, 'tt_matrix');
    return;
end

% Copy CONSTRUCTOR
if (nargin == 1) && isa(varargin{1}, 'tt_matrix')
    t=tt_matrix;
    t.tt = varargin{1}.tt;
    t.n = varargin{1}.n;
    t.m = varargin{1}.m;
    return;
end
if ( nargin == 3 && isa(varargin{1},'tt_tensor') && isa(varargin{2},'double') && isa(varargin{3},'double'))
    tt=varargin{1};
    n=varargin{2};
    m=varargin{3};
    d=tt.d;
    if ( numel(n) == 1 )
      n=n*ones(d,1);
    end
    if ( numel(m) == 1 )
      m=m*ones(d,1);
    end
    t=tt_matrix;
    t.n=n;
    t.m=m;
    t.tt=tt;
return
end
%From old format
if ( nargin == 1 ) && isa(varargin{1},'cell') 
    t=tt_matrix;
   tt=varargin{1};
   tt1=tt_mat_to_vec(tt); 
   t.tt=tt_tensor(tt1);
   d=t.tt.d;
   n=zeros(d,1);
   m=zeros(d,1);
   for i=1:d
      n(i)=size(tt{i},1);
      m(i)=size(tt{i},2);
   end
   t.n=n;
   t.m=m;
end
if ( nargin == 1 && isa(varargin{1},'double') || nargin ==2 && isa(varargin{1},'double') && isa(varargin{2},'double'))
    t=tt_matrix;
    b=varargin{1};
    if ( nargin == 2 && isa(varargin{2},'double')) 
      eps=varargin{2};
    else
      eps=1e-14;
      fprintf('Accuracy not specified, using %3.2e \n',eps);  
    end
    nm=size(b);
    d=numel(nm);
    d=d/2;
    n=nm(1:d);
    m=nm(d+1:2*d);
    prm=1:2*d; prm=reshape(prm,[d,2]); prm=prm';
    prm=reshape(prm,[1,2*d]); %Transposed permutation
    b=permute(b,prm); 
    if (numel(n.*m) == 1 )
      b=b(:);
    else
    b=reshape(b,n.*m);
    end
    tt=tt_tensor(b,eps);
    t.tt=tt;
    t.n=n';
    t.m=m';
end
if ( nargin == 3 && isa(varargin{1},'tt_tensor') && isa(varargin{2},'double') && isa(varargin{3},'double') )
    t=tt_matrix;
    tt=varargin{1};
    n=varargin{2};
    m=varargin{3};
    t.tt=tt;
    t.n=n;
    t.m=m;
end

% From a simple tt_matrix struct without class definition
if (nargin == 1) && isa(varargin{1}, 'struct')
    
    t.tt=tt_tensor(varargin{1}.tt);
    t.n=varargin{1}.n; %n-dimensions
    t.m=varargin{1}.m; %m-dimensions
    t = class(t, 'tt_matrix');
 
    return;
end;

return;
