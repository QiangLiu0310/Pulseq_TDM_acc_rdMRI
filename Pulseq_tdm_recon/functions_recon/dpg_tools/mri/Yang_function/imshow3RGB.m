function imshow3RGB(imgRGB,shape)
% imshow3(img, [ range, [shape )
%
% function to display series of images as a montage.
% 
% img - a 3D array representing a series of images
% range - window level (similarly to imshow
% shape - a 2x1 vector representing the shape of the motage
%
% Example: 
% 		im = repmat(phantom(128),[1,1,6]);
%		figure;
%		imshow3(im,[],[2,3]);
%
% (c) Michael Lustig 2012
for i=1:3
    img = imgRGB(:,:,:,i);
    [sx,sy,nc] = size(img);
    
    
    img = reshape(img,sx,sy*nc);
    img = permute(img,[2,3,1]);
    img = reshape(img,sy*shape(2),shape(1),sx);
    img = permute(img,[3,2,1]);
    img = reshape(img,sx*shape(1),sy*shape(2));
    
    img2(:,:,i)=img;
end

%imagesc(img,range); colormap(gray(256));axis('equal');
imshow(img2);
