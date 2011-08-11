function [y]=tt_rc(d,n,elem_fun,eps,varargin)
%[Y]=TT_RC(D,N,ARR,ELEM_FUN,EPS,[OPTIONS])
%The TT-Renormalization-Cross algorithm
%input is: the dimension of the array d, the 
%size vector n, the function that computes the prescribed element
%of an array, accuracy parameter eps
%and also additional arguments specified as 'rmax',10,'nswp',10
%'verb',true and so on
% 'x0', x0
%The elem function has syntax elem_fun(ind); all other additional
%information
%has to be pass as "anonymous" arguments


%Default parameters
if ( numel(n) == 1 )
  n=n*ones(1,d);
end
nswp=10;
rmax=1000;
rin=10; 
x0=[];
verb=true;
vec=false; %If there is a vectorization in options, i.e. if elem_fun
%supports vectorized computations of many elements at once
for i=1:2:length(varargin)-1
    switch lower(varargin{i})
        case 'nswp'
            nswp=varargin{i+1};
        case 'rmax'
            rmax=lower(varargin{i+1});
        case 'x0'
            x0=varargin{i+1};
        case 'verb'
            verb=varargin{i+1};
        case 'vec'
            vec=varargin{i+1};
        otherwise
            error('Unrecognized option: %s\n',varargin{i});
    end
end

if ( isempty(x0) )
   x0=tt_rand(n,d,rin); 
end

%The initial setup is done. 
%What we do first is we generate an initial index sets that we will expand
%the will be generated (deliberatedly) from the initial approximation. 
%The algorithm is essentially right-to-left qr + maxvol. That would
%generate right index sets. Left index set is then initialized as empty.

y=x0; %This is the approximate solution

ry=y.r;
n=y.n;
psy=y.ps;
cry=y.core;

yold=y;
rm=1;

pos1=psy(d+1)-1;
j_right=cell(d,1);
r_right=cell(d,1);
%Right-to-left qr & maxvol
for i=d:-1:2
    %fprintf('i=%d psy(%d)=%d psy(%d)=%d ry(%d)=%d n(%d)=%d ry(%d)=%d \n',i,i,psy(i),i+1,psy(i+1),i,ry(i),i,n(i),i+1,ry(i+1));
    cr=cry(psy(i):psy(i+1)-1);
    cr=reshape(cr,[ry(i)*n(i),ry(i+1)]);
    cr=cr*rm; 
    %cry1=cry;
    %cry1(psy(i):psy(i+1)-1)=cr(:);
    %y.core=cry1;
    %norm(nm*y-f)
    %keyboard

    ry(i+1)=size(cr,2); 
    cr=reshape(cr,[ry(i),n(i)*ry(i+1)]);
    cr=cr.'; %It seems to be needed
    [cr,rm]=qr(cr,cr); 
    indu=maxvol2(u); 
    r1=cr(indu,:);
    cr=cr/r1;
    rm=r1*rm; rm=rm.';
    rnew=size(cr,2);
    cr=cr.'; %n(i)*ry(i+1)xry(i)->ry(i)xn(i)xry(i+1) 
    %ry(i+1)=rnew;
    
    cry(pos1-ry(i+1)*n(i)*rnew+1:pos1)=cr(:); 
    pos1=pos1-ry(i+1)*n(i)*rnew; %Mda
    
    %Convert the index into a pair 
    %ind_loc=r_right+(j_right-1)*n(i)
    %j_right = 1 + floor(ind_loc/n(i));
    j_right{i}=floor(indu/n(i))+1;
    r_right{i}=indu - (j_right{i}-1)*n(i);
    
end
cr=cry(psy(1):psy(2)-1);
cr=reshape(cr,[ry(1)*n(1),ry(2)]);
cr=cr*rm; 
ry(2)=size(cr,2); 
cry(pos1-ry(1)*n(1)*ry(2)+1:pos1)=cr(:);
pos1=pos1-ry(1)*n(1)*ry(2)+1; %Fix it to the start
%And now truncate cry
cry=cry(pos1:end);

%Test that here that all is okey. (y has not changed, supposedly)
psy=cumsum([1;(n).*ry(1:y.d).*ry(2:y.d+1)]);
y.core=cry;
y.r=ry;
y.ps=psy;

%Now set up the left index set by a one sweep of the ``left-cross''.
%In the case of the starting rank equal to 1, it is OK --- the matrix
%will be initialized well. However, this is not a very good idea when
%the initial approximation has ``large ranks''. In this case,
%it is more convenient to provide initial index sets along with 
%the approximation, which should be ``interpolative''
%Also we will have to add the  ``restart logic''(LATER)

%This is just left-to-right qr & maxvol
%Not true --- here we already compute elements
%i.e. the step is as follows
%Compute indices of elements to be computed
%Compute these elements. QR & maxvol new core

cry_new=[];
i_left=cell(d,1);
r_left=cell(d,1);
ind_left=cell(d+1,1);
ind_right=cell(d+1,1); %The multiindices computed
%No need to recompute them all the time (?) 
for i=1:d
    %Compute the core related to (i,i+1) supercore.
    %For this we will need: 
    %left index set (i or i-1,
    %say --- i, that reduces modes (1,...,i-1) to ry(i)
    %right index set that reduces modes (i+2) ... d 
    %to ry(i+2) 
    %The final index set is obtained by taking all possible
    %combinations (i.e., 3-cycle) 
    %The first mode is deliberatedly considered, i.e. in fact 
    %we have d+2 indices. The computation now is performed only
    %in the middle of the pie
    l_cur=ind_left{i};
    r_cur=ind_right{i+1}; 
    index_set=zeros(d,ry(i),n(i),ry(i+1));
    for s1=1:ry(i)
      for s2=1:ry(i+2)
          for i1=1:n(i)
                iloc=[l_cur(:,s1); i1; r_cur(:,s2)];
               index_set(:,s1,i1,s2)=iloc(2:d+1); %Truncate border indices
          end
      end
    end
    %Compute the core
    index_set=reshape(index_set,numel(index_set)/d,d);
    if (vec )
      cr=elem_fun(index_set);
    else
      m=size(index_set,1);
      cr=zeros(m,1);
      for q=1:m
         cr(q)=elem_fun(index_set(q,:));
      end
    end
    %The new core is then QR-ed & maxvol-ed, also new multiindex set is
    %computed 
   cr=reshape(cr,[ry(i)*n(i),ry(i+1)]);
   [cr,rm]=qr(cr,0);
   indL=maxvol2(cr);
   
   r1=cr(indL,:);
   cr=cr/r1;
   cry_new=[cry_new,cr(:)]; 
   i_left{i}=floor(indL/n(i))+1;
   r_left{i}=indL - (i_left{i}-1)*n(i);

   %After that, new index 
end
psy=cumsum([1;(n).*ry(1:y.d).*ry(2:y.d+1)]);

%After that, the main TT-RC cross iteration starts (which also 
%allows to increase/decrease the rank). 









return