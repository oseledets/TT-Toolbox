% function [qtt]=tt_to_qtt3(tt, type, eps)
% Converts TT tensor to QTT format by compressing TT cores with
% approximation eps.
% type:
%   0 - vector 
%   1 - matrix


function [qtt]=tt_to_qtt3(tt, type, eps, max_r)

d_phys = size(tt,1);

if ((nargin<4)||(isempty(max_r)))
    max_r = [];
end;

n = zeros(d_phys,1);
new_d=zeros(d_phys,1);
positions=zeros(d_phys+1,1);
d=0;
for q=1:d_phys
    n(q)=size(tt{q}, 1);
    new_d(q)=round(log2(n(q)));
    d=d+new_d(q);
    positions(q+1)=d;    
end;
% if (type==0)
%     phys_ranks = tt_ranks(tt);
% else
%     phys_ranks = tt_ranks(tt_mat_to_vec(tt));
% end;

if (type==1)
    tt = tt_mat_to_vec(tt);    
end;
% QR to core 1
rv = 1; rnew=1;
for i=d_phys:-1:2
    core = tt{i};
    ncur = size(core,1); r1 = size(core,2); r2 = size(core,3);
    core = reshape(core, ncur*r1, r2);
    core = core*rv;
    core = reshape(permute(reshape(core, ncur, r1, rnew), [1 3 2]), ncur*rnew, r1);
    r2 = rnew;
    [q,rv]=qr(core,0); % size ncur*r2*rnew, rnew*r1
    rnew=min(r1,ncur*r2);
    tt{i}=permute(reshape(q, ncur, r2, rnew), [1 3 2]);
    rv = permute(rv, [2 1]);
end;
tt{1}=tt{1}*rv;
if (type==1)
    tt = tt_vec_to_mat(tt, n, n);
end;

if (type==0)
    tt{1}=reshape(tt{1}, size(tt{1},1),1,size(tt{1},2));
else
    tt{1}=reshape(tt{1}, size(tt{1},1), size(tt{1},2),1,size(tt{1},3));
end;

qtt = cell(d,1);

if (type==0)
    for q=1:d_phys
        tt{q}=permute(tt{q}, [2 1 3]);
        r1 = size(tt{q}, 1);
        r2 = size(tt{q}, 3);
        qsize = 2*ones(1,new_d(q));
        tt{q}=reshape(tt{q}, [r1 qsize r2]);
        
        cur_qtt = full_to_tt(tt{q}, eps/sqrt(d_phys-1), max_r);        
        if (size(cur_qtt,1)<(new_d(q)+2)) 
            dumme_kuh=cell(new_d(q)+2,1);
            for s=1:new_d(q)+1
                dumme_kuh{s}=cur_qtt{s};
            end;
            dumme_kuh{new_d(q)+2}=1;
            cur_qtt=dumme_kuh; 
        end;        
        
        r3 = size(cur_qtt{2}, 2);
        r4 = size(cur_qtt{2}, 3);
        cur_qtt{2} = permute(cur_qtt{2}, [2 1 3]);        
        cur_qtt{2} = cur_qtt{1}*reshape(cur_qtt{2}, r3, 2*r4);
        cur_qtt{2} = permute(reshape(cur_qtt{2}, r1, 2, r4), [2 1 3]);
        
%         cur_qtt{new_d(q)+2} has the sizes r2*r2' and is not orthogonal
%         convolve it with q+1 core of tt
        r2new = size(cur_qtt{new_d(q)+2}, 2);
        if (q<d_phys)
            core = tt{q+1};
            ncur = size(core,1); r3 = size(core,3);
            core = reshape(permute(core, [2 1 3]), r2, ncur*r3);
            core = (cur_qtt{new_d(q)+2}.')*core; % now size r2new*ncur*r3
            tt{q+1} = permute(reshape(core, r2new, ncur, r3), [2 1 3]);
        else
            % The last core: we can just convolve last fictitious r*r' core to
            % the whole train
            r3 = size(cur_qtt{new_d(q)+1}, 2);
            r4 = size(cur_qtt{new_d(q)+1}, 3);
            cur_qtt{new_d(q)+1} = reshape(cur_qtt{new_d(q)+1}, 2*r3, r4);
            cur_qtt{new_d(q)+1} = cur_qtt{new_d(q)+1}*conj(cur_qtt{new_d(q)+2}');
            cur_qtt{new_d(q)+1} = reshape(cur_qtt{new_d(q)+1}, 2, r3, r2);
        end;
        
        for p=(positions(q)+1):positions(q+1)
            qtt{p}=cur_qtt{p-positions(q)+1};
        end;
    end;
    qtt{1}=reshape(qtt{1}, size(qtt{1},1),size(qtt{1},3));
else
    for q=1:d_phys
        tt{q}=permute(tt{q}, [3 1 2 4]);
        r1 = size(tt{q}, 1);
        r2 = size(tt{q}, 4);
        qsize = 2*ones(1,2*new_d(q));
        tt{q}=reshape(tt{q}, [r1 qsize r2]);
        prm=zeros(1,new_d(q)*2+2);
        prm(1)=1;
        for i=1:new_d(q)
            prm(2*i-1+1)=i+1;
            prm(2*i+1)=new_d(q)+i+1;
        end;        
        prm(2*new_d(q)+2)=2*new_d(q)+2;
        tt{q}=permute(tt{q},prm);       
        sz1=4*ones(1,new_d(q)+2);
        sz1(1)=r1;
        sz1(new_d(q)+2)=r2;
        tt{q}=reshape(tt{q},sz1);
        
        cur_qtt=full_to_tt(tt{q},eps/sqrt(d_phys-1), max_r);
        if (size(cur_qtt,1)<(new_d(q)+2)) 
            dumme_kuh=cell(new_d(q)+2,1);
            for s=1:new_d(q)+1
                dumme_kuh{s}=cur_qtt{s};
            end;
            dumme_kuh{new_d(q)+2}=1;
            cur_qtt=dumme_kuh; 
        end;

        r3 = size(cur_qtt{2}, 2);
        r4 = size(cur_qtt{2}, 3);
        cur_qtt{2} = permute(cur_qtt{2}, [2 1 3]);        
        cur_qtt{2} = cur_qtt{1}*reshape(cur_qtt{2}, r3, 4*r4);
        cur_qtt{2} = permute(reshape(cur_qtt{2}, r1, 4, r4), [2 1 3]);

%         cur_qtt{new_d(q)+2} has the sizes r2*r2' and is not orthogonal
%         convolve it with q+1 core of tt
        r2new = size(cur_qtt{new_d(q)+2}, 2);
        if (q<d_phys)
            core = tt{q+1};
            ncur = size(core,1); mcur=size(core,2); r3 = size(core,4);
            core = reshape(permute(core, [3 1 2 4]), r2, ncur*mcur*r3);
            core = (cur_qtt{new_d(q)+2}.')*core; % now size r2new*ncur*r3
            tt{q+1} = permute(reshape(core, r2new, ncur, mcur, r3), [2 3 1 4]);
        else
            % The last core: we can just convolve last fictitious r*r' core to
            % the whole train
            r3 = size(cur_qtt{new_d(q)+1}, 2);
            r4 = size(cur_qtt{new_d(q)+1}, 3);
            cur_qtt{new_d(q)+1} = reshape(cur_qtt{new_d(q)+1}, 4*r3, r4);
            cur_qtt{new_d(q)+1} = cur_qtt{new_d(q)+1}*conj(cur_qtt{new_d(q)+2}');
            cur_qtt{new_d(q)+1} = reshape(cur_qtt{new_d(q)+1}, 4, r3, r2);
        end;
        
        for p=(positions(q)+1):positions(q+1)
            r1 = size(cur_qtt{p-positions(q)+1}, 2);
            r2 = size(cur_qtt{p-positions(q)+1}, 3);
            qtt{p}=reshape(cur_qtt{p-positions(q)+1}, 2, 2, r1, r2);
        end;
    end;
    qtt{1}=reshape(qtt{1}, size(qtt{1},1),size(qtt{1},2),size(qtt{1},4));    
end;

end