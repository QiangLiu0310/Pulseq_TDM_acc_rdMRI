dpgfullpath = mfilename('/rfanfs/pnl-zorro/home/lipeng/routines/functions_recon/');
dpgname = mfilename;
dpgpath0 = '/rfanfs/pnl-zorro/home/lipeng/routines/functions_recon/';%regexprep(dpgfullpath,dpgname,'');
dpgpath=[dpgpath0 '/dpg_tools/'];
addpath([ dpgpath '/mbin/' ]);
addpath([ dpgpath '/mathlib/conjgrad/' ]);
addpath([ dpgpath '/mathlib/tensor/' ]);
addpath([ dpgpath '/mri/mtlib/' ]);
addpath([ dpgpath '/mri/pmri/' ]);
addpath([ dpgpath '/mri/sie_pcm/' ]);
addpath([ dpgpath '/mri/siemens/' ]);
addpath([ dpgpath '/mri/siemens/mapVBVD/' ]);
addpath([ dpgpath '/nifti/' ]);
addpath([ dpgpath '/mri/myepi/' ]);
addpath([ dpgpath '/mri/Yang_function/' ]);
addpath([dpgpath0 '/mrir_toolbox/']);
clear dpgpath0 dpgpath dpgname dpgfullpath
