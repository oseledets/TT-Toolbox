function [out_val,out_ind,nevals]=tt_stat(tt, varargin)
% Compute statistics (min,max) of a tt_tensor via maxvol cross
%   [out_val,out_ind,nevals]=tt_stat(tt, varargin)
%
% Tries to compute values and indices of the sought quantities, listed in
% varargin. Possible choices are (case insensitive):
%   'LM', 'SM': largest, resp. smallest in magnitude entries,
%   'LR', 'SR': entries having largest, resp. smallest real parts, or
%   'LI', 'SI': entries having largest, resp. smallest imaginary parts.
% out_val is a 1 x N vector of output values, listed in the same order as
% in varargin, out_ind is a d x N matrix of indices, where d is the
% dimension of tt. The number of inspected elements is returned in nevals.
% 
% Additional parameters can be given in varargin in the form 
% 'param_name', param_value:
%       'kickrank': enrich tt by random subspaces of sizes kickrank. Larger
%       kickrank improves the accuracy in "difficult" cases, for the price
%       of larger neval. Default value is 1.
%       'nswp': number of alternating iterations (default 2)
%
% Examples:
%   Mimic the call to tt_abs_max, tt_max and tt_min simultaneously:
%       [out,ind,neval]=tt_stat(x, 'lm', 'lr', 'sr');
%   The same with kickrank set to 5:
%       [out,ind,neval]=tt_stat(x, 'lm', 'lr', 'sr', 'kickrank', 5);
%   Minimize the banana function on a grid 32768 x 32768:
%       n = 15;
%       x = tkron(tt_x(2,n)/2^n, tt_ones(2,n));
%       y = tkron(tt_ones(2,n), tt_x(2,n)/2^n);
%       f = (0.5-x).^2 + 100*(y-x.^2).^2;
%       v = tt_stat(f,'sm','nswp',15)
%   This should give values below 1e-8 in ~50% of runs.
%
% The difference with tt_*max* is the usage of dnr^2 entries 
% instead of dr^2 (outputs are sought over TT blocks, not density matrices).
% Besides, it can carry out several iterations with random enrichment.

kickrank = 1; 
nswp = 2; % these seem to be good params at least for 'LM'

if (numel(varargin)==0)
    error('at least one sought quantity is needed');
end;

% Distinguish soughts from parameters
soughts = cell(1,numel(varargin));
j = 1;
while (j<=numel(varargin))
    if    ((strcmpi(varargin{j}, 'lm'))||(strcmpi(varargin{j}, 'sm'))||...
           (strcmpi(varargin{j}, 'lr'))||(strcmpi(varargin{j}, 'sr'))||...
           (strcmpi(varargin{j}, 'li'))||(strcmpi(varargin{j}, 'si')))
       soughts{j} = varargin{j};
    end;
    if (strcmpi(varargin{j}, 'kickrank'))
        kickrank = varargin{j+1};
        j = j+1;
    end;
    if (strcmpi(varargin{j}, 'nswp'))
        nswp = varargin{j+1};
        j = j+1;
    end;
    j = j+1;
end;
soughts(cellfun('isempty', soughts))=[];

if (kickrank==0)
    nswp = 1;
end;

d = tt.d;
n = tt.n;
r = tt.r;
tt = core2cell(tt);
r0 = r; % r will become different below

ind_global = cell(d+1,1);
ind_local = cell(d+1,1); % for reducing units
out_val = nan*ones(1,numel(soughts));
out_ind = zeros(d,numel(soughts));
% reductions
XX0 = cell(d+1,1);
XX0{1} = 1; XX0{d+1} = 1;

% We need to subtract unitary vectors to exlude already considered indices
% and select fresh ones
% Prepare storages for them
units = cell(d,numel(soughts));
for j=1:numel(soughts)
    for i=1:d
        units{i,j} = zeros(n(i),1);
    end;
end;
XU = cell(d+1,numel(soughts));
for j=1:numel(soughts)
    XU{1,j}=1; XU{d+1,j}=1;
end;
nevals = 0;

% Start iteration
swp = 0; % Zeroth iteration will be just QR and MV
dir = -1;
i = d;
while (swp<=nswp)
    if (swp>0)
        % Compute statistics
        q = reshape(tt{i}, r0(i), n(i)*r0(i+1));
        q = XX0{i}*q;
        q = reshape(q, r(i)*n(i), r0(i+1));
        q = q*XX0{i+1};
        q = reshape(q, r(i)*n(i)*r(i+1), 1);
        nevals = nevals + r(i)*n(i)*r(i+1); % number of evaluations
        for j=1:numel(soughts)
            out_changed = false;
            if (strcmpi(soughts{j}, 'lm'))
                [val,ind]=max(abs(q));
                out_changed = (val>abs(out_val(j)))||(isnan(out_val(j)));
            end;
            if (strcmpi(soughts{j}, 'lr'))
                [val,ind]=max(real(q));
                out_changed = (val>real(out_val(j)))||(isnan(out_val(j)));
            end;
            if (strcmpi(soughts{j}, 'li'))
                [val,ind]=max(imag(q));
                out_changed = (val>imag(out_val(j)))||(isnan(out_val(j)));
            end;
            if (strcmpi(soughts{j}, 'sm'))
                [val,ind]=min(abs(q));
                out_changed = (val<abs(out_val(j)))||(isnan(out_val(j)));
            end;
            if (strcmpi(soughts{j}, 'sr'))
                [val,ind]=min(real(q));
                out_changed = (val<real(out_val(j)))||(isnan(out_val(j)));
            end;
            if (strcmpi(soughts{j}, 'si'))
                [val,ind]=min(imag(q));
                out_changed = (val<imag(out_val(j)))||(isnan(out_val(j)));
            end;
            if (out_changed)
                out_val(j) = q(ind);
                ind = tt_ind2sub([r(i), n(i), r(i+1)], ind);
                if (i>1)&&(i<d)
                    out_ind(:,j) = indexmerge(ind_global{i}(ind(1),:), ind(2), ind_global{i+1}(ind(3),:));
                elseif (i==1)
                    out_ind(:,j) = indexmerge(ind(2), ind_global{i+1}(ind(3),:));
                else % i==d
                    out_ind(:,j) = indexmerge(ind_global{i}(ind(1),:), ind(2));
                end;
                % Update units vectors
                for k=1:d
                    units{k,j}=zeros(n(k),1);
                    units{k,j}(out_ind(k,j))=1;
                end;
                % and their reductions
                for k=1:i-1                    
                    XU{k+1,j} = XU{k,j}*units{k,j}.';
                    XU{k+1,j} = reshape(XU{k+1,j}, r(k)*n(k), 1);
                    XU{k+1,j} = XU{k+1,j}(ind_local{k+1}, :);
                end;
                for k=d:-1:i+1
                    XU{k,j} = units{k,j}*XU{k+1,j};
                    XU{k,j} = reshape(XU{k,j}, 1, n(k)*r(k+1));
                    XU{k,j} = XU{k,j}(:, ind_local{k});
                end;
            end;
        end;
        
        % Subtract units to get fresh indices
        for j=1:numel(soughts)
            qu = XU{i,j}*reshape(units{i,j}*XU{i+1,j}, 1, n(i)*r(i+1));
            qu = reshape(qu, r(i)*n(i)*r(i+1),1);
            q = q-qu*out_val(j);
        end;
    end;
    
    if (dir>0)&&(i<d)
        q = reshape(q, r(i)*n(i), r(i+1));
        rr = randn(r(i)*n(i), kickrank);
        [q,~]=qr([q, rr],0);
        rnew = size(q,2);
        
        % Maxvol
        ind = maxvol2(q);
        qq = q(ind,:);
        q = q/qq;
        % Indices
        ind_local{i+1} = ind; % for units
        ind_global{i+1} = indexmerge(ind_global{i}, (1:n(i))');
        ind_global{i+1} = ind_global{i+1}(ind,:);
        r(i+1) = rnew;
        % Reduction
        XX0{i+1} = reshape(tt{i}, r0(i), n(i)*r0(i+1));
        XX0{i+1} = XX0{i}*XX0{i+1};
        XX0{i+1} = reshape(XX0{i+1}, r(i)*n(i), r0(i+1));
        XX0{i+1} = XX0{i+1}(ind,:);
        % For units
        for j=1:numel(soughts)
            XU{i+1,j} = XU{i,j}*reshape(units{i,j}, 1, n(i));
            XU{i+1,j} = reshape(XU{i+1,j}, r(i)*n(i), 1);
            XU{i+1,j} = XU{i+1,j}(ind, :);
        end;
    elseif (dir<0)&&(i>1)
        % QR and right indices
        if (swp==0)
            q = reshape(tt{i}, r(i), n(i)*r(i+1));
        else
            q = reshape(q, r(i), n(i)*r(i+1));
        end;
        rr = randn(n(i)*r(i+1), kickrank*double(swp>0));
        [q,rv]=qr([q.',rr], 0);
        rv = rv(:,1:r(i)).';
        q = q.';
        rnew = size(q,1);
        % Maxvol
        ind = maxvol2(q.');
        qq = q(:, ind);
        q = qq\q;
        rv = rv*qq;
        % Indices
        ind_local{i} = ind; % for units
        ind_global{i} = indexmerge((1:n(i))', ind_global{i+1});
        ind_global{i} = ind_global{i}(ind,:);
        if (swp==0)
            tt{i} = reshape(q, rnew, n(i), r(i+1));
            % Cast rv to the next block
            q = reshape(tt{i-1}, r(i-1)*n(i-1), r(i));
            tt{i-1} = reshape(q*rv, r(i-1), n(i-1), rnew);
            r0(i) = rnew;
        end;
        r(i) = rnew;
        % Reduction
        XX0{i} = reshape(tt{i}, r0(i)*n(i), r0(i+1));
        XX0{i} = XX0{i}*XX0{i+1};
        XX0{i} = reshape(XX0{i}, r0(i), n(i)*r(i+1));
        XX0{i} = XX0{i}(:,ind);   
        for j=1:numel(soughts)
            XU{i,j} = units{i,j}*XU{i+1,j};
            XU{i,j} = reshape(XU{i,j}, 1, n(i)*r(i+1));
            XU{i,j} = XU{i,j}(:, ind);
        end;
    end;
    
    i = i+dir;    
    if ((dir>0)&&(i==d))||((dir<0)&&(i==1))
        dir = -dir;
        swp = swp+1;
    end;
end;

end

function [J]=indexmerge(varargin)
% Merges two or three indices in the little-endian manner
sz1 = max(size(varargin{1},1),1);
sz2 = max(size(varargin{2},1),1);
sz3 = 1;
if (nargin>2) % Currently allows only 3
    sz3 = max(size(varargin{3}, 1), 1);
end;
% J1 goes to the fastest index, just copy it
J1 = repmat(varargin{1}, sz2*sz3, 1);
% J2 goes to the middle
J2 = reshape(varargin{2}, 1, []);
J2 = repmat(J2, sz1, 1); % now sz1 ones will be the fastest
J2 = reshape(J2, sz1*sz2, []);
J2 = repmat(J2, sz3, 1);
J = [J1,J2];
if (nargin>2)
    % J3 goes to the slowest
    J3 = reshape(varargin{3}, 1, []);
    J3 = repmat(J3, sz1*sz2, 1); % now sz1 ones will be the fastest
    J3 = reshape(J3, sz1*sz2*sz3, []);
    J = [J,J3];
end;
end
          