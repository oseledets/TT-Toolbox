function [tt]=qr(tt,op)
%[TT]=QR(TT,OP)
%Implements the TT-QR algorithm for the orthogonalization of the
%TT-representation. OP can be 'LR' (for the left-right orthogonalization)
%and 'RL' for the right-left orthogonalization
d=tt.d;
r=tt.r;
cr=tt.core;
pos=tt.ps;
n=tt.n;
pos1=1;
core0=cr(1:r(1)*n(1)*r(2));
if ( strcmp(op,'lr') || strcmp(op,'LR')) 
  %Orthogonalization from left-to-tight
  for i=1:d-1
   core0=reshape(core0,[r(i)*n(i),r(i+1)]);
   [core0,ru]=qr(core0,0); %nrm(i+1)=norm(ru,'fro');
   %ru=ru;%./nrm(i+1);
   core1=cr(pos(i+1):pos(i+2)-1);
   core1=reshape(core1,[r(i+1),n(i+1)*r(i+2)]);
   core1=ru*core1;
   r(i+1)=size(core0,2);
   cr(pos1:pos1-1+r(i)*n(i)*r(i+1))=core0(:);
   cr(pos1+r(i)*n(i)*r(i+1):pos1+r(i)*n(i)*r(i+1)+r(i+1)*n(i+1)*r(i+2)-1)=core1(:);
   core0=core1;
   pos1=pos1+r(i)*n(i)*r(i+1);
  end
 pos1=pos1+r(d)*n(d)*r(d+1)-1;
 cr=cr(1:pos1); %Truncate storage if required
 tt.core=cr; tt.r=r; 
 tt.ps=cumsum([1;tt.n.*tt.r(1:d).*tt.r(2:d+1)]);
 return
else
  error('Option is not supported (yet)');    
end
