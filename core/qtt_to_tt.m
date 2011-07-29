% function [tt]=qtt_to_tt(qtt, new_d, type)
% Produce TT vector in physical dimensions from given QTT vector by
% a convolution of each new_d(q) cores in QTT to a full vector, which represents
% corresponding TT factor of q-th core.
% type:
%     0 - vector,
%     1 - matrix


function [tt]=qtt_to_tt(qtt, new_d, type)

d_phys = max(size(new_d));
% if (d_phys==1) d_phys=size(new_d,2); end;

d = size(qtt,1);
tt = cell(d_phys, 1);
if (type==0)
    qtt_ranks = tt_ranks(qtt);
else
    qtt_ranks = tt_ranks(tt_mat_to_vec(qtt));
end;

phys_ranks = zeros(d_phys-1,1);
positions = zeros(d_phys-1,1);
cur_pos=0;
for q=1:d_phys-1
    cur_pos=cur_pos+new_d(q);
    phys_ranks(q)=qtt_ranks(cur_pos);
    positions(q)=cur_pos;
end;

if (type==0)
    if (new_d(1)>1)
        tt{1}=zeros(2^(new_d(1)), phys_ranks(1));
        cur_qtt = cell(new_d(1), 1);
        for p=1:new_d(1)-1
            cur_qtt{p}=qtt{p};
        end;
        for k=1:phys_ranks(1)
            cur_qtt{new_d(1)}=qtt{new_d(1)}(:,:,k);
            tt{1}(:,k)=qtt_to_full(cur_qtt);
        end;
    else
        tt{1}=qtt{1};
    end;
    
    for q=2:d_phys-1
        if (new_d(q)>1)
            tt{q}=zeros(2^(new_d(q)), phys_ranks(q-1), phys_ranks(q));
            cur_qtt = cell(new_d(q), 1);
            for k1=1:phys_ranks(q-1)
                for k2=1:phys_ranks(q)
                    cur_qtt{1}=permute(qtt{positions(q-1)+1}(:,k1,:), [1 3 2]);
                    for p=2:new_d(q)-1
                        cur_qtt{p}=qtt{positions(q-1)+p};
                    end;
                    cur_qtt{new_d(q)}=qtt{positions(q)}(:,:,k2);
                    tt{q}(:,k1,k2)=qtt_to_full(cur_qtt);
                end;
            end;
        else
            tt{q}=qtt{positions(q)};
        end;
    end;
    
    if (new_d(d_phys)>1)
        tt{d_phys}=zeros(2^(new_d(d_phys)), phys_ranks(d_phys-1));
        cur_qtt = cell(new_d(d_phys), 1);
        for p=2:new_d(d_phys)
            cur_qtt{p}=qtt{positions(d_phys-1)+p};
        end;
        for k=1:phys_ranks(d_phys-1)
            cur_qtt{1}=permute(qtt{positions(d_phys-1)+1}(:,k,:), [1 3 2]);
            tt{d_phys}(:,k)=qtt_to_full(cur_qtt);
        end;
    else
        tt{d_phys}=qtt{positions(d_phys-1)+1};
    end;
    
else
    n1 = tt_size(qtt);
    n2 = tt_size(tt_transp(qtt));
    
    if (new_d(1)>1)
        tt{1}=zeros(prod(n1(1:new_d(1))), prod(n2(1:new_d(1))), phys_ranks(1));
        cur_qtt = cell(new_d(1), 1);
        for p=1:new_d(1)-1
            cur_qtt{p}=qtt{p};
        end;
        for k=1:phys_ranks(1)
            cur_qtt{new_d(1)}=qtt{new_d(1)}(:,:,:,k);
            tt{1}(:,:,k)=nd_to_full(cur_qtt);
        end;
    else
        tt{1}=qtt{1};
    end;
    
    for q=2:d_phys-1
        if (new_d(q)>1)
            tt{q}=zeros(prod(n1(sum(new_d(1:q-1))+1:sum(new_d(1:q)))),prod(n2(sum(new_d(1:q-1))+1:sum(new_d(1:q)))), phys_ranks(q-1), phys_ranks(q));
%             tt{q}=zeros(2^(new_d(q)),2^(new_d(q)), phys_ranks(q-1), phys_ranks(q));
            cur_qtt = cell(new_d(q), 1);
            for k1=1:phys_ranks(q-1)
                for k2=1:phys_ranks(q)
                    cur_qtt{1}=permute(qtt{positions(q-1)+1}(:,:,k1,:), [1 2 4 3]);
                    for p=2:new_d(q)-1
                        cur_qtt{p}=qtt{positions(q-1)+p};
                    end;
                    cur_qtt{new_d(q)}=qtt{positions(q)}(:,:,:,k2);
                    tt{q}(:,:,k1,k2)=nd_to_full(cur_qtt);
                end;
            end;
        else
            tt{q}=qtt{positions(q)};
        end;
    end;
    
    if (new_d(d_phys)>1)
        tt{d_phys}=zeros(prod(n1(end-new_d(d_phys)+1:end)), prod(n2(end-new_d(d_phys)+1:end)), phys_ranks(d_phys-1));
        cur_qtt = cell(new_d(d_phys), 1);
        for p=2:new_d(d_phys)
            cur_qtt{p}=qtt{positions(d_phys-1)+p};
        end;
        for k=1:phys_ranks(d_phys-1)
            cur_qtt{1}=permute(qtt{positions(d_phys-1)+1}(:,:,k,:), [1 2 4 3]);
            tt{d_phys}(:,:,k)=nd_to_full(cur_qtt);
        end;
    else
        tt{d_phys}=qtt{positions(d_phys-1)+1};
    end;
end;

end