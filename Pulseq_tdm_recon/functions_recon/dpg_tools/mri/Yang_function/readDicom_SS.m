% 2019/04/12 JiYang,read the dicom files under one folder which contains
% many dicom files. One sequence one foder.
% The matrix size(phase, frequency and slice number) for each sequence should be same.
% single sequence (SS)
function [img, info] =readDicom_SS(dir_folder1)

disp(['---',' Folder, Read data begin---',dir_folder1]);

img_path_list=dir(dir_folder1);
img_num=length(img_path_list);

for j=3:img_num
    temp_data=dicomread(fullfile(dir_folder1,img_path_list(j).name));
    img(:,:,j-2)=double(temp_data);
end


disp(['---',' Folder, Read data end---',dir_folder1]);
info=dicominfo(fullfile(dir_folder1,img_path_list(3).name),'UseDictionaryVR',true);
