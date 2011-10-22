function [ev,vec]=norm2(tt,varargin)
%[ev,vec]=norm2(tt,options)
%Computes and approximation to the 2-norm of a tt-matrix
nswp=40;
tol=0.1; %Relative norm change
eps=1e-2; %Computations accuracy
rmax=100; %Maximal rank
sz=size(tt); sz=sz(:,2);
r0=2;
x0=[];
matvec='full';
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'tol'
            tol=lower(varargin{i+1});
        case 'x0'
            x0=varargin{i+1};
        case 'matvec'
            matvec=varargin{i+1};
        case 'eps'
            eps=varargin{i+1};
        case 'rmax'
            rmax=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end
if ( isempty(x0) )
   x0=tt_rand(sz,ndims(tt),r0);
end
er=2*tol;
ev_old=inf;
x=x0;
while ( er > tol) 
   x=x/norm(x);
   if ( strcmp(matvec,'full') )
     y=tt*x; 
     y=round(y,eps); 
     ev=dot(y,x);
     er=abs(ev-ev_old)/abs(ev);
     ev_old=ev;
     x=y;
   else
     error('Krylov-type method is not supported yet! \n');
   end
end
return
end