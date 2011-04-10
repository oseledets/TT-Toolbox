function s = latex(M,varargin)
%LATEX   Print a matrix in LaTeX tabular format.
%   LATEX(M) prints out the numeric matrix M in a LaTeX tabular
%   format. The '&' character appears between entries in a row, '\\'
%   is appended to the ends of rows, and each entry is set in math
%   mode. Complex numbers are understood, and exponentials will be
%   converted to a suitable format.
%
%   LATEX(M,'nomath') does not include the $$ needed to put each 
%   entry in math mode (e.g., for use with the amsmath matrix modes).
%   
%   LATEX(M,FMT) uses a format specifier FMT of the SPRINTF type for
%   each entry.
%   
%   LATEX(M,FMT1,FMT2,...) works through the given format specifiers
%   on each row of M. If fewer are given than the column size of M,
%   the last is used repeatedly for the rest of the row.
%   
%   S = LATEX(M,...) does not display output but returns a character
%   array S.
%   
%   Examples:
%     latex( magic(4) )
%     latex( magic(4), '%i', 'nomath' )
%     latex( magic(4), '%i', '%.2f' )
%   
%   See also SPRINTF, SYM/LATEX.

%   Copyright 2002 by Toby Driscoll. Last updated 12/06/02.

if ~isa(M,'double')
  error('Works only for arrays of numbers.')
elseif ndims(M) > 2
  error('Works only for 2D arrays.')
end

if nargin < 2
  fmt = {'%#.5g'};
  mathstr = '$';
else
  fmt = varargin;
  idx = strmatch('nomath',fmt);
  if isempty(idx)
    mathstr = '$';
  else  
    mathstr = '';
    fmt = fmt([1:idx-1 idx+1:end]);
    if isempty(fmt), fmt = {'%#.5g'}; end
  end 
end

% Extend the format specifiers.
[m,n] = size(M);
if n > length(fmt)
  [fmt{end:n}] = deal(fmt{end});
end
  
% Create one format for a row.
rowfmt = '';
for p = 1:n
  % Remove blanks.
  thisfmt = deblank(fmt{p});

  % Add on imaginary part if needed.
  if ~isreal(M(:,p)) 
    % Use the same format as for the real part, but force a + sign for
    % positive numbers. 
    ifmt = thisfmt;
    j = findstr(ifmt,'%');
    if ~any(strcmp(ifmt(j+1),['-';'+';' ';'#']))
      ifmt = [ifmt(1:j) '+' ifmt(j+1:end)];
    else
      ifmt(j+1) = '+';
    end
    ifmt = [ifmt 'i'];
    thisfmt = [thisfmt ifmt];
  end

  % Add to row.
  rowfmt = [rowfmt mathstr thisfmt mathstr ' & '];
end

% After last column, remove column separator and put in newline.
rowfmt(end-1:end) = [];
rowfmt = [rowfmt '\\\\\n'];

% Use it.
A = M.';
if isreal(M)
  S = sprintf(rowfmt,A);
else
  S = sprintf(rowfmt,[real(A(:)) imag(A(:))].');
end

% Remove extraneous imaginary part for real entries.
if ~isreal(M)
  zi = sprintf(ifmt,0);
  S = strrep(S,zi,blanks(length(zi)));
end

% Remove NaNs.
S = strrep(S,'$NaN$','--');
S = strrep(S,'NaN','--');

% Convert 'e' exponents to LaTeX form. This is probably really slow, but
% what can you do without regular expressions?
S = strrep(S,'e','E');
ex = min(findstr(S,'E'));
while ~isempty(ex)
  % Find first non-digit character. Where is ISDIGIT?
  j = ex+2;
  while ~isempty(str2num(S(j))) & ~strcmp(S(j),'i')
    j = j+1;
  end

  % This strips off leading '+' and zeros.
  num = sprintf('%i',str2num(S(ex+1:j-1)));
  
  ee = ['\times 10^{' num '}'];
  S = [S(1:ex-1) ee S(j:end)];
  
  ex = ex + min(findstr(S(ex+1:end),'E'));
end

% For good form, remove that last '\\'.
S(end-2:end-1) = '  ';

% Display or output?
if nargout==0
  disp(S)
else
  s = S;
end
