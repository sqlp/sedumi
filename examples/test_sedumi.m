function test_sedumi (do_rebuild, do_exit)
% TEST_SEDUMI  Verify minimal functionality of SeDuMi using an example SDP.

if (nargin < 2)
  do_rebuild = false;
  do_exit = false;
end

try
  old_dir = cd ('..');
  if (do_rebuild)
    install_sedumi -rebuild
  else
    install_sedumi
  end
  cd (old_dir)
  
  example_path = mfilename ('fullpath');
  example_path = example_path(1:max ([1, strfind(example_path, filesep ())]) - 1);
  
  test_problems = {'arch0.mat', 'control07.mat', 'nb.mat', ...
    'OH_2Pi_STO-6GN9r12g1T2.mat', 'trto3.mat'};
  
  for i = 1:length (test_problems)
    line_sep = repmat ('-', 1, 75);
    fprintf ('\n%s\n-- Test ''%s'' (%2d/%2d)\n%s\n\n', line_sep, ...
      test_problems{i}, i, length (test_problems), line_sep);
    
    load (fullfile (example_path, test_problems{i}), 'At', 'b', 'c', 'K');
    [~,~,info] = sedumi (At, b, c, K);
    if (info.pinf ~= 0 || info.dinf ~= 0 || info.numerr == 2)
      fprintf ('\n\nSeDuMi TEST FAILED: Could not solve ''%s'' example.\n', ...
        test_problems{i});
      if (do_exit)
        exit (1);
      else
        return;
      end
    end
    if (info.numerr == 1)
      warning ('SeDuMi reported reduced accuracy.')
    end
  end
catch
  disp ('SeDuMi TEST FAILED: Compilation or runtime error.');
  if (do_exit)
    exit (1);
  else
    return;
  end
end

end
