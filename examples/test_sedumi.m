function test_sedumi (do_rebuild, do_exit)
% TEST_SEDUMI  Verify minimal functionality of SeDuMi using an example SDP.

if (nargin < 2)
  do_rebuild = false;
  do_exit = false;
end

try
  old_dir = cd ('..')
  if (do_rebuild)
    install_sedumi -rebuild
  else
    install_sedumi
  end
  cd (old_dir)
  load (fullfile ('examples', 'nb.mat'));
  [~,~,info] = sedumi (At, b, c, K);
  if (sum ([info.pinf, info.dinf, info.numerr])) % 0 = no error
    disp ('SeDuMi TEST FAILED: Could not solve ''nb.mat'' example.');
    if (do_exit)
      exit (1);
    else
      return;
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
