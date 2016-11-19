numsteps = 5;
varincrements = 1.5;
P = phantom(128); 
%P = rgb2gray(imread('small_beetle.jpg'));
P = double(P)./double(max(P(:)));
R = radon(P,0:179);

subplot(3, numsteps+1, numsteps+2), imshow(P), title('Original')
subplot(3, numsteps+1, 2*numsteps+3), imshow(P), title('Original')

for i = 0:numsteps-1
    var = varincrements*i;
    noise = var*randn(size(R));
    noisyRad = R + noise;
    
    I1 = iradon(noisyRad,0:179);
    I2 = iradon(noisyRad,0:179,'linear','none');
    subplot(3,numsteps+1,2+i), imshow(noisyRad,[]), title(sprintf('Radon %.2f', var))
    subplot(3,numsteps+1,numsteps+3+i), imshow(I1,[]), title(sprintf('Filtered backproj. %.2f', var))
    subplot(3,numsteps+1,2*numsteps+4+i), imshow(I2,[]), title(sprintf('Unfiltered backproj. %.2f', var))
end
