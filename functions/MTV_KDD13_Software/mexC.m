fprintf('make sure you have set up the C complier, see "help mex"\n');
fprintf('mex C files, please wait....\n');
current_dir = cd;
cd([current_dir '/Code']);
mex mtv2.c
mex mtvp2.c COMPFLAGS="$COMPFLAGS -openmp"  LINKFALGS="$LINKFALGS -openmp"
mex mtvp3.c COMPFLAGS="$COMPFLAGS -openmp"  LINKFALGS="$LINKFALGS -openmp"
cd(current_dir);