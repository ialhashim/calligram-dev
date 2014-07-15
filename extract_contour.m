image_filename = 'bunny.jpg';
I = imread(image_filename);
BW = im2bw(I);
B = bwboundaries(BW);
C = flip(B,2)
csvwrite([image_filename, '.txt'], C);

% plot extracted boundary
imshow (BW);

