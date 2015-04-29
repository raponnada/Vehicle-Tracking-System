% hw1_2
% Chris Pennock
% Computer Vision, Fall 2004, NYU, Professor  Bregler
% 10/21/04

% Take as arguments the root of an image name (eg 'SRI/1bw')
% and the last number of the images (eg 19)
% This function will load SRI1.png and let the user choose a
% point of interest to track.
% It will sucessively load SRI1.png, SRI2.png, etc
% and track the square region of height and width 40 centered on the
% chosen pixel.
% It will save a series of images that have the
% region of interest highlighted.

% returns nothing
function r = hw1_2(imageroot, numimages)

% make a directory to hold the results
resultdirname = 'results';
mkdir(resultdirname);

% load first image, display to user, get their point of interest
imstring = strcat(imageroot, '1.jpg');
F = imread(imstring);
graymap;
image(F);

% make the X and Y lookup matrices for interpolation
imsize = size(F);
[X,Y] = meshgrid(1:1:imsize(2), 1:1:imsize(1));
X = double(X);
Y = double(Y);

% create region of interest - this is rounded to the nearest integer values
% of x and y
roicenter = ginput(1)
% r is half of the side length of the ROI box
r = 10; 
roitopy = round(roicenter(2)-r);
roiboty = round(roicenter(2)+r);
roileftx = round(roicenter(1)-r);
roirightx = round(roicenter(1)+r);
roi = F(roitopy:roiboty, roileftx:roirightx);

% draw the first image, with the roi selected
% and write out that image
Fbox = drawBox(F, roitopy, roiboty, roileftx, roirightx);
image(Fbox);
imtowrite = strcat(resultdirname, '/1.jpg');
imwrite(Fbox, imtowrite, 'jpg');

% make the convolved images Fx and Fy
% and all other constants necessary for calclutaion with sigma = 3
sigma = 3;
g = gauss(sigma);
dg = dgauss(sigma);
kernelx = g'*dg;
kernely = kernelx';

fx = conv2(double(F), kernelx, 'same');
fy = conv2(double(F), kernely, 'same');
fxroi = double(fx(roitopy:roiboty, roileftx:roirightx));
fyroi = double(fy(roitopy:roiboty, roileftx:roirightx));

C11 = sum(sum(fxroi.^2));
C12 = sum(sum(fxroi.*fyroi ));
C21 = C12;
C22 = sum(sum(fyroi.^2));
C = [C11, C12; C21, C22];

% make another set of all constants, for sigma = 1.5
% for the later iterations
sigma2 = 1.5;
g2 = gauss(sigma2);
dg2 = dgauss(sigma2);
kernelx2 = g2'*dg2;
kernely2 = kernelx2';

fx2 = conv2(double(F), kernelx2, 'same');
fy2 = conv2(double(F), kernely2, 'same');
fxroi2 = double(fx2(roitopy:roiboty, roileftx:roirightx));
fyroi2 = double(fy2(roitopy:roiboty, roileftx:roirightx));

C11_2 = sum(sum(fxroi2.^2));
C12_2 = sum(sum(fxroi2.*fyroi2 ));
C21_2 = C12_2;
C22_2 = sum(sum(fyroi2.^2));
C_2 = [C11_2, C12_2; C21_2, C22_2];

% if this is set out here, then every frame will start with the uv of the
% previous frame
uv = [0 0];

% For all images in the sequence, find where the roi has moved to.
for i = 2:numimages
    
    % open the image
    imname = strcat(imageroot, num2str(i))
    imstring = strcat(imname, '.jpg')
    G = imread(imstring);
    
    % initial value of the warped version of G is just G.
    warpedG = G;
    
    % Set these to have some values going into the 'while' test below.
    groi = warpedG(roitopy:roiboty, roileftx:roirightx);
    ftroi = double(groi)-double(roi);
    
    % double version of G, necessary for later interpolation.
    Gdub = double(G);
    
    % iteration counter
    j = 0;

    % ftcount is the amount of pixel in the roi
    [fty, ftx] = size(ftroi);
    ftcount = fty*ftx;
    
    % While error > some amount, and iterations < 10 (so it will always stop,
    % even if the error never reaces an ideal point).
    % Since ftroi ranges from 0 to 255, normalize it to 0->1 with ./255, 
    % then get the sqaured error with .^2
    while ( (j<10) && ( (sum(sum((ftroi./255).^2))/ftcount) > .0005 ) )
        
        % The warp step.
        % uv is not reset after each frame, it just keeps getting incremented.
        % So warp the original G by this total amount of uv at each iteration.
        XI = double(X+uv(1));
        YI = double(Y+uv(2));
        warpedG = interp2(X,Y,Gdub,XI,YI,'cubic');
        groi = warpedG(roitopy:roiboty, roileftx:roirightx);
        ftroi = double(groi)-double(roi);
        
        % for the first 2 iterations, use a large sigma
        if j<2
            D11 = sum(sum(ftroi.*fxroi));
            D12 = sum(sum(ftroi.*fyroi));
            D = [D11; D12];
            invC = inv(C);
            % incruv is the incremental addition to uv at each step
            incruv = invC*D;
        % then use a smaller sigma
        else 
            D11 = sum(sum(ftroi.*fxroi2));
            D12 = sum(sum(ftroi.*fyroi2));
            D = [D11; D12];
            invC = inv(C_2);
            incruv = invC*D;
        end
        % update uv by adding the incremental uv.
        uv = uv - incruv';        
        j = j+1;
    end
    
    % write an image, with tracking box, to the designated results directory
    Gbox = drawBox(G, roitopy+round(uv(2)), roiboty+round(uv(2)), roileftx+round(uv(1)), roirightx+round(uv(1)));
    imtowrite = strcat(resultdirname, '/', num2str(i), '.jpg')
    imwrite(Gbox, imtowrite, 'jpg');
end

% the final image, with a tracking box
figure
graymap;
Gbox = drawBox(G, roitopy+round(uv(2)), roiboty+round(uv(2)), roileftx+round(uv(1)), roirightx+round(uv(1)));
image(Gbox)

% return nothing
r = '0 kmph';



% This function takes a matrix (the image, in gray values 0->255)
% and the x and y positions of the box lines to draw.
% It returns the image with a box drawn
% It also interpolates the box for subpixel accuracy. This makes the
% animation look nicer and be more accurate.
function boxed = drawBox(M, ytop, ybot, xleft, xright)
boxed = M;
boxed(ytop, xleft:xright) = 255;
boxed(ybot, xleft:xright) = 255;
boxed(ytop:ybot, xleft) = 255;
boxed(ytop:ybot, xright) = 255;

    