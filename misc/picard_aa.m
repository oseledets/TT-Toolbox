function [sol]=picard_aa(x0,f,eps,varargin)
%[sol]=picard_aa(x0,fun,eps,options)
%Accelerates the convergence of the fixed point iteration
%x := f(x)
%by Anderson acceleration
%x0 --- starting vector, f --- function to evaluate
%eps -- tolerance
%Optional parameters (supplied as keys, i.e. 'niters', 100)
%niters --- maximal number of iterations, default: 100
%restart --- dimension of the AA space, default: 100
%warmup --- number of fixed point iterations before AA starts, default: 1
%verb --- print convergence information, default: true
%fixed_iter --- if set to true, only the fixed point iteration is run, default  --- false 
%process_fun --- if found, apply this function to compute some information about the solution
%(i.e. draw it, or compute some functional from it, etc)
warmup=1;
niters=100;
m=100;
verb=true;
draw=false;
fixed_iter=false;
condtol=1e12;
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'maxiters'
            niters=varargin{i+1};
        case 'warmup'
            warmup=varargin{i+1};
        case 'restart'
            m=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'fixed_iter'
            fixed_iter=varargin{i+1};
        case 'process_fun'
           draw=true;
           process_fun=varargin{i+1}; 
        case 'condtol'
           condtol=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

for q=1:warmup
  x0 = f(x0);
end
%Simple fixed iterations
x=f(x0); 
er=2*eps;
k=2;
if ( fixed_iter ) %Just do a simple fixed point iteration
while ( k < niters && er  > eps )
  k = k + 1;
  fx=f(x);
  er=norm(fx-x);
  x=fx;
  if ( verb )  
     fprintf('it=%d er=%3.2e  \n',k,er);
  end
end
sol=x;
return
end
%res=[f(x0),f(x)]; %The array for the residuals
%xc=[x0,x];
xc=[x0,x];
res=[f(x0)-x0,f(x)-x];
while ( k < niters && er > eps )
    er=norm(res(:,end));
    if ( verb ) 
              x1=xc(:,end); x0=xc(:,end-1);
       er_loc=norm(x1-x0)/norm(x1);
       fprintf('it=%d  er=%3.2e er_loc=%3.2e \n',k,er,er_loc);
    end
    if ( size(res,2) > m || cond(res) > condtol) %Too much memory
      res(:,1)=[];
      xc(:,1)=[];
    end
    %res=res(:,end-m0+1:end);
    %xc=xc(:,end-m0+1:end);
    m0=size(res,2);
    
    ed=eye(m0); zd=eye(m0);
    ed(m0,:)=[]; zd(1,:)=[];
    gd=zd-ed;
    Fk=res*gd'; 
    rhs=res(:,end);
    gm=Fk \ rhs;
    dx=xc*gd';
    %New x
    xnew=rhs + xc(:,end)-(dx+Fk)*gm;
    
    if ( draw ) 
       process_fun(xnew);
    end
    
    xc=[xc,xnew];
    res=[res,f(xnew)-xnew];
    k=k+1;
end
sol=xc(:,end);
return
end
