%load ../../../../codes/gamess/pes_radical.txt;
%load pes_radical.txt; %Load QFF file in the format that I prepared by hand
                      %Actually, it can be done by a simple script further
                      %on the GAMESS output
%a=pes_radical; clear pes_radical;
load formad_full.txt;
a=formad_full; clear formad_full;
%The format is: MODE; GRID; Value
pos=1;
nmodes=12; 
ngrid=16;
nmodes=nmodes-6; f=nmodes;
ap=cell(nmodes,1);
for i=1:nmodes
   ap{i}=zeros(ngrid,1);
end
while ( a(pos,2) == -1 )
   mode=a(pos,1);
   pos=pos+1;
   grid=a(pos,1);
   pos=pos+1;
   val=a(pos,1);
   pos=pos+1;
   ap{mode-6}(grid)=val;
end


%Read the coupling
vcoup=cell(nmodes-6,nmodes-6);
for i=1:nmodes
  for j=1:nmodes
    vcoup{i,j}=zeros(ngrid,ngrid);
  end
end
while ( pos <= size(a,1))
  
  mode1=a(pos,1); mode2=a(pos,2); mode1=mode1-6; mode2=mode2-6;
  pos=pos+1;
  grid1=a(pos,1); grid2=a(pos,2);
  pos=pos+1;
  val=a(pos,1);
  pos=pos+1;
  vcoup{mode2,mode1}(grid2,grid1)=val;
end

%This was recovered from the file; Later on will do 
%it automatically

% MODE # 1   QRANGE=   32.668884     DEL-Q=    4.355851   
% MODE # 2   QRANGE=   33.055664     DEL-Q=    4.407422
% MODE # 3   QRANGE=   42.896804     DEL-Q=    5.719574
% MODE # 4   QRANGE=   45.595380     DEL-Q=    6.079384
% MODE # 5   QRANGE=   50.402975     DEL-Q=    6.720397
% MODE # 6   QRANGE=   52.113100     DEL-Q=    6.948413

qrange=[32.668884,33.055664,42.896804,45.595380,50.402975,52.113100];
qrange=qrange(end:-1:1);
hqrange=2*qrange/(ngrid-1);
grx=cell(nmodes,1);
for i=1:nmodes
   grx{i}=-qrange(i)+(i-1)*hqrange(i);
end


%Generate the V1 and V2 potential energies
e=tt_tensor(ones(ngrid,1));
V1=[];
for i=1:f
  w=tt_tensor(ap{i});
  for j=1:i-1
    w=kron(e,w);
  end
  for j=i+1:f
    w=kron(w,e);
  end
  V1=V1+w;
  V1=round(V1,1e-8);
end

ng=16*ones(1,f); %That is true somehow
%Create v2 coupling positions
V2=[];
eps=1e-8;
for i=1:f
  for j=i+1:f
      w=tt_tensor(vcoup{i,j},eps); 
      %Now put auxiliary ranks between i...j 
      %i1a1 a1i2a2 [a2b2 b2c2] 
      d=w.d; n=w.n; r=w.r; ps=w.ps; cr=w.core;
      u1=chunk(w,1,1); 
      u2=chunk(w,2,2);
      rq=rank(u1,2); 
      for k=i+1:j-1
         ep=tt_tensor; ep.d=1; ep.n=ng(k); ep.r=[rq;rq];
         cr=zeros(rq,ng(k),rq);
         for p=1:ng(k)
             cr(:,p,:)=eye(rq);
         end
         cr=cr(:); ep.core=cr; ep.ps=[1,rq*ng(k)*rq];
         u1=kron(u1,ep);
      end
      w=kron(u1,u2);
      for k=1:i-1
        w=kron(e,w);
      end
      for k=j+1:f
        w=kron(w,e);
      end
         V2=V2+w; V2=round(V2,eps);
  end
end

%Vfull=V2+(f-2)*V1; Vfull=round(Vfull,1e-8); %if V=V1+V2+V3; V(q1,q2)=V1+V2 V(q1,q3)=V1+V3 V(q2,q3)=V2+V3 2*V1+2*V2+2*V3 
 
%Vfull=V1;
Vfull=V2+V1; 
%full(Vfull(:,8,8,8,8,8))
%full(V1(:,8,8,8,8,8))
Vfull=round(Vfull,1e-8);
%return
%% Grid details
%Generate Laplacian
lp0=tt_qlaplace_dd(4); lp0=tt_matrix(lp0); lp0=full(lp0); lp0=tt_matrix(lp0);
e0=eye(16); e0=tt_matrix(e0);
lp=[];
for i=1:f
  w=lp0;
  for j=1:i-1
    w=kron(e0,w);
  end
  for j=i+1:f
    w=kron(w,e0);
  end
     lp=lp+w/hqrange(i)^2;
  lp=round(lp,1e-12);
end

%lp=kron(lp0,kron(e0,e0))/hx(1)^2+kron(e0,kron(lp0,e0))/hx(2)^2+...
%    kron(e0,kron(e0,lp0))/hx(3)^2;
%lp=round(lp,1e-7); lp=lp/2; 


mat=lp*0.5+diag(Vfull); mat=round(mat,eps);

%Start the solver

v=dmrg_eigb(mat,1,1e-6,[],[],10);

%Eigenvalue
%ev(q)=dot(mat*v,v);
%q=q+1;


%Convert to cm
%amu= 1829; cm=2.194747000000000e+05;
%hm= sqrt(ap{1}(1))+sqrt(ap{2}(1))+sqrt(ap{3}(1));
%hm=hm/sqrt(2);
%sqrt(ev/(2*amu))*cm
