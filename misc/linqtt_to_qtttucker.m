function [qtu] = linqtt_to_qtttucker(tt, dims, eps)
% function [qtu] = linqtt_to_qtttucker(tt, dims, eps)
% Computes the qtt_tucker representation from the linear (Q)TT with the
% accuracy eps.
% dims specifies, how many dimensions include in each tucker factor

dphys = numel(dims);
dtt = tt.d;
qtu = qtt_tucker;
qtu.dphys = dphys;
tucks = cell(dphys, 1);
sz = cell(dphys, 1);
qtucore = tt_tensor;
qtucore.d = dphys;
qtucore.n = zeros(dphys, 1);
qtucore.r = [1; zeros(dphys-1,1); 1];
qtucore.ps = ones(dphys+1,1);
qtucore.core = 0;

[tt,nrm]=qr(tt, 'rl');
tt{1} = nrm*tt{1};

for i=1:dphys
    cf = chunk(tt, 1, dims(i));
    if (i<dphys)
        tt = chunk(tt, dims(i)+1, dtt);
    end;
    dtt = dtt-dims(i);
    
    curn = cf.n;
    sz{i} = curn;
    rc1 = cf.r(1);
    rc2 = cf.r(dims(i)+1);
    cf = tt_reshape(cf, [rc1; curn; rc2], eps/sqrt(dphys*2));
    cf = move_tt_block(cf, 1, dims(i)+1, eps/sqrt(dphys*2));
    cf = tt_reshape(cf, [curn; rc1*rc2]);
    [cf, nrm]=qr(cf, 'lr');
    rtuck = cf.r(dims(i)+1);    
    curcr = cf{dims(i)+1};
    curcr = curcr*nrm;
    % Orthogonal tucker factor is ready
    tucks{i} = chunk(cf, 1, dims(i));
    % Now prepare the core
    if (i<dphys)
        curcr = reshape(curcr, rtuck*rc1, rc2);
        [curcr, rv]=qr(curcr, 0);
        % Orthogonal tucker core is ready
        % Convolve the RV to the rest of linear TT
        cr_in_tt = tt{1};
        n2 = size(cr_in_tt, 2);
        rtt2 = size(cr_in_tt, 3);
        cr_in_tt = reshape(cr_in_tt, rc2, n2*rtt2);
        cr_in_tt = rv*cr_in_tt;
        rc2 = size(curcr,2);
        tt{1} = reshape(cr_in_tt, rc2, n2, rtt2);
    end;
    curcr = reshape(curcr, rtuck, rc1, rc2);
    qtucore{i} = permute(curcr, [2, 1, 3]);
end;

qtu.core = qtucore;
qtu.tuck = tucks;
qtu.sz = sz;

end