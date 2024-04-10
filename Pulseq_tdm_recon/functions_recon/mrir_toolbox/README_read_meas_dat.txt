
matlab supplies a MEX version of 'typecast.m' that is extremely fast
but it is located in a private directory.

to see the directory name, run the following on the MATLAB command line:
   fullfile(fileparts(which('typecast')), 'private')


'read_meas_dat.m' makes use of the 'typecastc.mex*' binary, but to do
so it needs to be in the MATLAB path.

when installing 'read_meas_dat.m' for the first time, create a
subdirectory named "private" in the same directory as
'read_meas_dat.m'. then create symbolic links to all of the
"typecastc" MEX files in this private directory. this will allow
'read_meas_dat.m' to call typecastc directly. (optionally these files
can be copied since the file sizes are negligibly small.) we copy all
of the files since different systems may call different versions of
the MEX library (e.g., MATLAB supplies 32- and 64-bit versions of the
MEX binaries).



NOTE: i have heard that typecastc may not be available in newer
versions of MATLAB (7.11?), or maybe 'typecast.m' is now faster.


NOTE 2: also included is the binary for typecastx, which can be found
here:

http://www.mathworks.com/matlabcentral/fileexchange/17476-typecast-and-typecastx-c-mex-functions

for MATLAB running on windows, download the MEX C code and compile for
the windows environment.
