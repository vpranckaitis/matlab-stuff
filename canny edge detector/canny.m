function canny(filename)
%CANNY Algorithm from http://docs.opencv.org/2.4/doc/tutorials/imgproc/imgtrans/canny_detector/canny_detector.html
  if (nargin == 0)
    filename = 'zebra-small.jpg';  
  end

  % Input image
  A = imread(filename);
  figure(1); image(A); boldTitle('Original image');
  
  % Grayscale image
  B = rgb2gray(A);
  figure(2); image(gray2rgb(B)); boldTitle('Grayscale image');
  
  % Applying Gaussian filter (Gaussian blur) 
  Bg = applyFilter(B, gaussianFilter(5, 2));
  figure(3); image(gray2rgb(Bg)); boldTitle('After Gaussian filter');
  
  % Finding gradient strength and direction
  gradFilter = [-1 0 1; -2 0 2; -1 0 1];
  Gx = applyFilter(Bg, gradFilter);
  Gy = applyFilter(Bg, gradFilter');
  G = sqrt(Gx.*Gx +Gy.*Gy);
  Th = atan(Gy ./ Gx);
  figure(4); image(gray2rgb(G)); boldTitle('Gradient strengths');
  figure(5); imageQuiver(Gx, Gy); boldTitle('Gradient directions');
  
  % Non-maximum suppression
  Gm = nonMaximumSupression(G, Th); 
  figure(6); image(gray2rgb(Gm)); boldTitle('After non-maximum supression');
  
  % Filter by upper and lower threshold
  upper = 150; lower = upper / 3;
  Gt = filterByThreshold(Gm, upper, lower); 
  figure(7); image(gray2rgb(Gt)); boldTitle('After threshold filtering');
end

function [K] = gaussianFilter(size, std)
  % 'size' describes matrix K size, standard deviation 'std' describes 
  % blurring amount
  range=double((-idivide(int32(size),2)):(idivide(int32(size),2)));
  [X, Y] = meshgrid(range, range);
  % https://en.wikipedia.org/wiki/Gaussian_filter#Definition
  K = exp(-(X.^2+Y.^2)/(2*std^2)) / (2 * 3.14 * std^2); 
  K = K ./ sum(sum(K)); % Normalize sum to 1
end

function [B] = applyFilter(A, K)
  A = double(A); 
  l = size(K, 1);
  l2 = idivide(int32(l), 2);
  B = zeros(size(A));
  for i=1:size(B,1)
    for j=1:size(B,2)
      ri = (i-l2):(i+l2);
      ri = min(max(ri, 1), size(A, 1));
      rj = (j-l2):(j+l2);
      rj = min(max(rj, 1), size(A, 2));
      B(i,j)=sum(sum(A(ri, rj).*K));
    end
  end
end

function [Gm] = nonMaximumSupression(G, Th)
  Thd = floor(Th / (pi/4)) * 45;
  Gm = zeros(size(G));
  for i=2:size(Gm,1)-1
      for j=2:size(Gm,2)-1
        switch (Thd(i, j))
          case -90
            if (G(i,j) > G(i+1,j) && G(i,j) >= G(i-1,j))
              Gm(i,j) = G(i,j);  
            end
          case -45
            if (G(i,j) > G(i+1,j-1) && G(i,j) >= G(i-1,j+1))
              Gm(i,j) = G(i,j);  
            end
          case 0
            if (G(i,j) > G(i,j+1) && G(i,j) >= G(i,j-1))
              Gm(i,j) = G(i,j);  
            end
          case 45
            if (G(i,j) > G(i+1,j+1) && G(i,j) >= G(i-1,j-1))
              Gm(i,j) = G(i,j);  
            end
        end
      end
  end
end

function [Gt] = filterByThreshold(G, upper, lower) 
  Gt = zeros(size(G));
  for i=2:size(G,1)-1
    for j=2:size(G,2)-1
      
      if G(i,j) > upper % If is above upper threshold 
        Gt(i,j) = G(i,j);
      elseif G(i,j) > lower % If is between upper and lower threshold
        for ii=-1:1 
          for jj=-1:1
            if G(i+ii,j+jj) > upper % If has adjacent pixel above upper threshold
              Gt(i,j) = G(i,j);  
            end
          end
        end
      end
      
    end
  end
end

function [B] = rgb2gray(A)
  c = [0.299; 0.587; 0.114];
  B = A(:,:,1)*c(1) + A(:,:,2)*c(2) + A(:,:,3)*c(3);
end

function [B] = gray2rgb(A)
  B = zeros([size(A) 3]); 
  B(:, :, 1) = A;
  B(:, :, 2) = A;
  B(:, :, 3) = A;
  B = uint8(B);
end

function imageQuiver(Gx, Gy)
  scale = 3;
  Gys = Gy(1:scale:end,1:scale:end);
  Gxs = Gx(1:scale:end,1:scale:end);
  [X, Y] = meshgrid(size(Gxs,1):-1:1, 1:size(Gxs,2));
  quiver(Y', X', -Gxs, Gys);
end

function boldTitle(str)
  title(str,'fontweight','bold','fontsize',14);
end