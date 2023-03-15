function test_sedumi (do_rebuild, do_exit)
% TEST_SEDUMI  Verify minimal functionality of SeDuMi using an example SDP.
%
%   Input:
%     do_rebuild - Rebuild SeDuMi
%     do_exit    - Exit Matlab/GNU Octave on failure (suiteable for continuous
%                  integration CI).

try
  old_dir = cd ('..');
  if ((nargin == 1) && do_rebuild)
    install_sedumi -rebuild
  else
    install_sedumi
  end
  cd (old_dir)

  example_path = mfilename ('fullpath');
  example_path = example_path(1:max ([1, strfind(example_path, filesep ())]) - 1);

  % Name and optimal function value
  test_problems = { ...
    'arch0.mat',                  -5.665170e-01; ...
    'control07.mat',              -2.062510e+01; ...
    'nb.mat',                     -5.070309e-02; ...
    'OH_2Pi_STO-6GN9r12g1T2.mat',  7.946708e+01; ...
    'trto3.mat',                  -1.279999e+04; ...
    'quantum.mat', 				  -0.75395345};

  tol = 1e-6;

  for i = 1:size (test_problems, 1)
    line_sep = repmat ('-', 1, 75);
    fprintf ('\n%s\n-- Test ''%s'' (%2d/%2d)\n%s\n\n', line_sep, ...
      test_problems{i,1}, i, size (test_problems, 1), line_sep);

    load (fullfile (example_path, test_problems{i,1}), 'At', 'b', 'c', 'K');
    [x, y, info] = sedumi (At, b, c, K);
    cx = c' * x;
    by = b' * y;
    cx_relerr = abs (cx - test_problems{i,2}) / abs (test_problems{i,2});
    by_relerr = abs (by - test_problems{i,2}) / abs (test_problems{i,2});
    if ((cx_relerr >= tol) || (by_relerr >= tol) ...
        || (info.pinf ~= 0) || (info.dinf ~= 0) || (info.numerr == 2))
      fprintf ('\n\tTEST FAILED\n');
      fprintf ('\n\t\tExpected = %e (   tol = %.1e)', test_problems{i,2}, tol);
      fprintf ('\n\t\t     c*x = %e (relerr = %.1e)', cx, cx_relerr);
      fprintf ('\n\t\t     b*y = %e (relerr = %.1e)\n', by, by_relerr);
      if ((nargin == 2) && do_exit)
        exit (1);
      else
        return;
      end
    end
    if (info.numerr == 1)
      fprintf ('\n\tWARNING: SeDuMi reported reduced accuracy.\n');
    end
    fprintf ('\n\tTEST PASSED\n');
  end
catch
  fprintf ('\n\tTEST FAILED: Compilation or runtime error.\n');
  if ((nargin == 2) && do_exit)
    exit (1);
  else
    return;
  end
end

end

