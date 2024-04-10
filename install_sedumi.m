function install_sedumi( varargin )

%SeDuMi installation script
%
%   install_sedumi
%
%      Basic usage.
%
%   install_sedumi -nopath
%
%      Do not add SeDuMi directories to the search path.
%
%   install_sedumi -rebuild
%
%      Rebuild binary mex-files SeDuMi.
%
%   install_sedumi -rebuild 'mex -O -largeArrayDims %s -lmwblas'
%
%      Rebuild binary mex-files SeDuMi using that mex command template
%      including ' %s ' as placeholder for input files.
%
%
% Heavy rewrite by Michael C. Grant for SeDuMi 1.3.4
% Now detects if binaries are already
%
% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA

need_rebuild = any (strcmp (varargin, '-rebuild'));
no_path = any (strcmp (varargin, '-nopath'));
mex_template = strncmp (varargin, 'mex ', 4);
if (any (mex_template))
  mex_template = varargin{mex_template};
  % Minimal validation.
  % "contains" would be better, but introduced in Matlab R2016b.
  if (isempty (strfind (mex_template, ' %s '))) %#ok<STREMP>
    error (['The mex compilation command must contain '' %s '' as ', ...
      'placeholder for the source files.']);
  end
else
  mex_template = '';
end

targets64 = {...
  'bwblkslv.c sdmauxFill.c sdmauxRdot.c',...
  'choltmpsiz.c',...
  'cholsplit.c',...
  ['dpr1fact.c auxfwdpr1.c sdmauxCone.c sdmauxCmp.c sdmauxFill.c', ...
  ' sdmauxScalarmul.c sdmauxRdot.c blkaux.c'],...
  'symfctmex.c symfct.c',...
  'ordmmdmex.c ordmmd.c',...
  'quadadd.c',...
  'sqrtinv.c sdmauxCone.c',...
  'givensrot.c auxgivens.c sdmauxCone.c',...
  'urotorder.c auxgivens.c sdmauxCone.c sdmauxTriu.c sdmauxRdot.c',...
  ['psdframeit.c reflect.c sdmauxCone.c sdmauxRdot.c sdmauxTriu.c', ...
  ' sdmauxScalarmul.c'],...
  ['psdinvjmul.c reflect.c sdmauxCone.c sdmauxRdot.c sdmauxTriu.c', ...
  ' sdmauxScalarmul.c blkaux.c'],...
  'bwdpr1.c sdmauxCone.c sdmauxRdot.c',...
  'fwdpr1.c auxfwdpr1.c sdmauxCone.c sdmauxScalarmul.c',...
  'fwblkslv.c sdmauxScalarmul.c',...
  'qblkmul.c sdmauxScalarmul.c',...
  'blkchol.c blkchol2.c sdmauxFill.c sdmauxScalarmul.c',...
  'vecsym.c sdmauxCone.c',...
  'qrK.c sdmauxCone.c sdmauxRdot.c sdmauxScalarmul.c',...
  'finsymbden.c sdmauxCmp.c',...
  'symbfwblk.c',...
  'ddot.c sdmauxCone.c sdmauxRdot.c sdmauxScalarmul.c',...
  'partitA.c sdmauxCmp.c',...
  'getada1.c sdmauxFill.c',...
  'getada2.c sdmauxCone.c sdmauxRdot.c sdmauxFill.c',...
  ['getada3.c spscale.c sdmauxCone.c sdmauxRdot.c sdmauxScalarmul.c', ...
  ' sdmauxCmp.c'],...
  'adendotd.c sdmauxCone.c',...
  'adenscale.c',...
  'extractA.c',...
  'sortnnz.c sdmauxCmp.c',...
  'iswnbr.c',...
  'incorder.c',...
  'findblks.c sdmauxCone.c sdmauxCmp.c',...
  ['invcholfac.c triuaux.c sdmauxCone.c sdmauxRdot.c sdmauxTriu.c', ...
  ' sdmauxScalarmul.c blkaux.c'],...
  };

% Get first file of "targets64" and replace ".c" with mex-extension.
% Keep order of "targets64".
mext = mexext;
mex_binaries = cellfun (@strtok, targets64, 'UniformOutput', false);
mex_binaries = cellfun (@(x) strrep (x, '.c', ['.', mext]), ...
  mex_binaries, 'UniformOutput', false);

fs = filesep();
sedumi_path = mfilename ('fullpath');
sedumi_path = sedumi_path(1:max ([1, strfind(sedumi_path, fs)]) - 1);

% Note the use of 0.01 here. That's because version 7 had more than 10
% minor releases, so 7.10-7.14 need to be ordered after 7.01-7.09.
VERSION  = [1, 0.01] * sscanf (version (), '%d.%d');
COMPUTER = computer ();
ISOCTAVE = (exist ('OCTAVE_VERSION', 'builtin') == 5);
mdir = '';
if (ISOCTAVE)
  prog = 'Octave';
  page_output_immediately (true, 'local');
  switch COMPUTER
  case 'x86_64-pc-linux-gnu'
      mdir = 'o_lin';
  case 'x86_64-apple-darwin21.6.0'
      mdir = 'o_maci';
  case 'aarch64-apple-darwin23.4.0'
      mdir = 'o_maca';
  case 'i686-w64-mingw32'
      mdir = 'o_win';
  case 'x86_64-w64-mingw32'
      mdir = 'o_win';
  otherwise
    error(sprintf('Unexpected computer type: %s', COMPUTER))
  end
else
  prog = 'Matlab';
end

line_sep = repmat ('-', 1, 75);
fprintf ('\n%s\nSeDuMi installation script\n', line_sep);
fprintf ('   Directory: %s\n', sedumi_path);
fprintf ('   %s %s on %s\n', prog, version (), computer ());
disp (line_sep);

% We don't want to rebuild the binaries if they're already present, unless
% the user has specifically asked for a rebuild. Frankly, we can't
% guarantee that rebuilding will work.
if (~need_rebuild)
  fprintf ('Looking for existing binaries...');
  nfound = [0, 0];
  for k = 1 : length(targets64)
      targ = mex_binaries{k};
      if exist( [ sedumi_path, fs, targ ], 'file' )
          nfound(1) = nfound(1) + 1;
      elseif ~isempty(mdir) && exist( [ sedumi_path, fs, mdir, fs, targ ], 'file' )
          nfound(2) = nfound(2) + 1;
      end
  end
  if sum(nfound) == 0
    fprintf ('none found; building...\n');
    need_rebuild = true;
  elseif sum(nfound) < length(targets64)
    fprintf ('incomplete set found.\n');
    disp (line_sep);
    error (['%s\n', ...
      'Some of the binaries for this platform were found, but some\n', ...
      'were missing as well. This may mean your download was corrupt;\n', ...
      'consider downloading and unpacking SeDuMi again. Otherwise, to\n', ...
      'try rebuilding the MEX files yourself, run this command:\n\n', ...
      '    install_sedumi -rebuild\n%s\n\n'], line_sep, line_sep);
  else
    fprintf ('found!\n');
    fprintf ('   If for some reason you need to rebuild the binaries,');
    fprintf (' use this command:\n');
    fprintf ('      install_sedumi -rebuild\n');
  end
end

if (need_rebuild)
  mdir = '';
  disp ('Attempting to recompile the SeDuMi binaries:');

  % Customization by providing a mex template.
  flags = {};
  libs = {};
  if ~isempty(mex_template)
    found = false;
    [cmd, remain] = strtok(mex_template);
    while ~isempty(remain)
      [tok, remain] = strtok(remain);
      if found
        libs{end+1} = tok;
      elseif tok == '%s'
        found = true;
      else
        flags{end+1} = tok;
      end
    end
  elseif (ISOCTAVE)
    % Matlab mex optimization '-O' corresponds to gcc '-O2',
    % whereas calling '-O' would result in gcc '-O1'.
    cmd = 'mkoctfile';
    flags{end+1} = '--mex';
    flags{end+1} = '-O2';
    flags{end+1} = '-DOCTAVE';
    flags{end+1} = '-Wall';
    if (ismac ())
      % Assume Homebrew (https://brew.sh/) installation.
      % https://stackoverflow.com/questions/50634727/dyld-library-not-loaded-usr-local-opt-openblas-lib-libopenblasp-r0-2-20-dylib
      % this glob covers /opt/homebrew and /usr/local installations
      homebrew = glob('/*/*/Cellar/openblas/*/include');
      flags{end+1} = ['-I', homebrew{1}];
      libs{end+1}  = ['-L', strrep(homebrew{1}, '/include', '/lib'), ' -lopenblas'];
    elseif (ispc ())
      libs{end+1}  = '-lopenblas';
    else
      % Including the default OpenBLAS path works for most Octave
      % installations and does not harm if not present.
      flags{end+1} = '-I/usr/include/openblas';
    end
  else % Matlab
    cmd = 'mex';
    % The last Matlab release with a 32 bit version was R2015b.
    IS64BIT  = strcmp (COMPUTER(end-1:end), '64');
    flags{end+1} = '-O';  % optimize
    if (IS64BIT && (VERSION >= 7.03)) % R2006b
      flags{end+1} = '-largeArrayDims';
    elseif (VERSION < 7.03)           % R2006b
      flags{end+1} = '-DmwIndex=int';
      flags{end+1} = '-DmwSize=int';
      flags{end+1} = '-DmwSignedIndex=int';
      flags{end+1} = '-DFWRAPPER';
    end
    if ((VERSION >= 7) && ispc ())  % R14 (2004)
      if IS64BIT
        dirval = 'win64';
      else
        dirval = 'win32';
      end
      libdir = fullfile (matlabroot (), 'extern', 'lib', dirval);
      if exist (fullfile (libdir, 'microsoft'), 'dir')
        libdir = fullfile (libdir, 'microsoft');
        found = true;
      elseif exist (fullfile (libdir, 'msvc60'), 'dir')
        libdir = fullfile (libdir, 'msvc60');
        found = true;
      elseif exist (fullfile (libdir, 'lcc'), 'dir')
        libdir = fullfile (libdir, 'lcc');
        found = true;
      end
      if found
        libs{end+1} = ['-L''', strrep(libdir, '\', '\\'), ''''];
      end
    end
    if (VERSION >= 7.05)  % R2007b
      libs{end+1} = '-lmwblas';
    else
      libs{end+1} = '-lmwlapack';
    end
  end

  mex_template = sprintf(' %s', cmd, flags{:}, '%s', libs{:});
  fprintf ('Template:%s\n', mex_template);

  % Move to SeDuMi root directory and compile.
  failed = false;
  olddir = cd (sedumi_path);
  for i = 1:length(targets64)
    files = {};
    remain = targets64{i};
    while ~isempty(remain),
      [fname, remain] = strtok(remain);
      files{end+1} = fname;
    end
    fprintf('  %s: %s\n', mex_binaries{i}, targets64{i});
    try
      if ISOCTAVE
        mkoctfile(flags{:}, files{:}, libs{:});
      else
        mex(flags{:}, files{:}, libs{:})
      end
    catch
      failed = true;
    end
  end
  cd (olddir);
  if failed
    error (['%s\nSeDuMi was not successfully installed.\n', ...
      'Please attempt to correct the errors and try again.\n\n', ...
      'For example, try setting the compilation command:\n\n', ...
      '  install_sedumi -rebuild ''%s'''], line_sep, mex_template);
  end
elseif (~isempty(mdir) && ~exist(fullfile(sedumi_path, mdir)))
  mdir = '';
end

% Add SeDuMi to search path, if not forbidden by "-nopath".
if (~no_path)
  disp (line_sep);
  disp ('Adding SeDuMi to the search path:');

  search_path_cstr = regexp (path (), pathsep (), 'split');
  if (ispc ())  % MS Windows is not case-sensitive.
    is_on_search_path = @(x) any (strcmpi (x, search_path_cstr));
  else
    is_on_search_path = @(x) any (strcmp (x, search_path_cstr));
  end

  fprintf ('   Base directory...       ');
  if (~is_on_search_path (sedumi_path))
    addpath (sedumi_path);
    fprintf ('added.\n');
  else
    fprintf ('already there.\n');
  end
  if (~isempty(mdir))
    fprintf ('   MEX directory... ');
    spath = fullfile(sedumi_path, mdir);
    if (~exist(spath, 'dir'))
      fprintf('not present.\n');
    elseif (~is_on_search_path (spath))
      addpath (spath);
      fprintf ('added.\n');
    else
      fprintf ('already there.\n');
    end
  end
  fprintf ('   Conversion directory... ');
  if (~is_on_search_path (fullfile (sedumi_path, 'conversion')))
    addpath (fullfile (sedumi_path, 'conversion'));
    fprintf ('added.\n');
  else
    fprintf ('already there.\n');
  end
  fprintf ('   Examples directory...   ');
  if (~is_on_search_path (fullfile (sedumi_path, 'examples')))
    addpath (fullfile (sedumi_path, 'examples'));
    fprintf ('added.\n');
  else
    fprintf ('already there.\n');
  end
end

disp (line_sep);
disp ('SeDuMi has been successfully installed.');
disp ('For more information, type "help sedumi" or see the user guide.');
disp (line_sep);
disp ('');
end
