%loading the 4 images of different intensity functions
I1 = rgb2gray(imread('picture1.jpeg'));
I2 = rgb2gray(imread('picture2.jpeg'));
I3 = rgb2gray(imread('picture3.jpeg'));
I4 = rgb2gray(imread('picture4.jpeg'));
%% crop region of interest
[J1 RECT] = imcrop(I1);
J2 = imcrop(I2, RECT);
J3 = imcrop(I3, RECT);
J4 = imcrop(I4, RECT);
%% calculate phase
Num = double(J4) - double(J2);
Den = double(J1) - double(J3);
PHI = atan2(Num, Den);
%% Phase unwrapping
IM = PHI;
im_mag   = abs(IM);                  %Magnitude image
im_phase = angle(IM);                %Phase image
%% Replace with your mask (if required)
mag_max = max(im_mag(:));
indx1 = find(im_mag < 0.1*mag_max);  %Intensity = mag^2, so this = .04 threshold on the intensity
im_mask = ones(size(IM));
im_mask(indx1) = 0;                  %Mask
if(~exist('im_mask','var'))
  im_mask = ones(size(IM));          %Mask (if applicable)
end
figure; imagesc(im_mag.*im_mask),   colormap(gray), axis square, axis off, title('Initial masked magnitude'); colorbar;
figure; imagesc(im_phase.*im_mask), colormap(gray), axis square, axis off, title('Initial masked phase'); colorbar;
im_unwrapped = nan(size(IM));        %Initialze the output unwrapped version of the phase
adjoin = zeros(size(IM));            %Zero starting matrix for adjoin matrix
unwrapped_binary = zeros(size(IM));  %Binary image to mark unwrapped pixels
%%
im_phase_quality = PhaseDerivativeVariance_r1(im_phase);
%% Automatically (default) or manually identify starting seed point on a phase quality map 
minp = im_phase_quality(2:end-1, 2:end-1); minp = min(minp(:));
maxp = im_phase_quality(2:end-1, 2:end-1); maxp = max(maxp(:));
if(0)    % Chose starting point interactively
  figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), colorbar, axis square, axis off; title('Phase quality map'); 
  uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
  [xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm
  colref = round(xpoint);
  rowref = round(ypoint);
  close;  % close the figure;
else   % Chose starting point = max. intensity, but avoid an edge pixel
  [rowrefn,colrefn] = find(im_mag(2:end-1, 2:end-1) >= 0.99*mag_max);
  rowref = rowrefn(1)+1; % choose the 1st point for a reference (known good value)
  colref = colrefn(1)+1; % choose the 1st point for a reference (known good value)
end
%% Unwrap
IM_unwrapped(rowref,colref) = im_phase(rowref,colref);                          %Save the unwrapped values
unwrapped_binary(rowref,colref,1) = 1;
if im_mask(rowref-1, colref, 1)==1;  adjoin(rowref-1, colref, 1) = 1; end       %Mark the pixels adjoining the selected point
if im_mask(rowref+1, colref, 1)==1;  adjoin(rowref+1, colref, 1) = 1; end
if im_mask(rowref, colref-1, 1)==1;  adjoin(rowref, colref-1, 1) = 1; end
if im_mask(rowref, colref+1, 1)==1;  adjoin(rowref, colref+1, 1) = 1; end
im_unwrapped = GuidedFloodFill(im_phase, im_mag, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap
% Plot images
figure; imagesc(im_mag),       colormap(gray), colorbar, axis square, axis off; title('QG Magnitude image'); 
figure; imagesc(im_phase),     colormap(gray), colorbar, axis square, axis off; title('QG Wrapped phase'); 
figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis square, axis off, title('QG Phase quality map'); colorbar;
figure; imagesc(im_unwrapped), colormap(gray), colorbar, axis square, axis off; title('QG Unwrapped phase');

%%
function derivative_variance = PhaseDerivativeVariance_r1(IM_phase, varargin)
    [r_dim,c_dim] = size(IM_phase);
    if nargin>=2                                    %Has a mask been included? If so crop the image to the mask borders to save computational time
        IM_mask = varargin{1};
        [maskrows,maskcols] = find(IM_mask);          %Identify coordinates of the mask
        minrow = min(maskrows)-1;                     %Identify the limits of the mask 
        maxrow = max(maskrows)+1;
        mincol = min(maskcols)-1;
        maxcol = max(maskcols)+1;
        width = maxcol-mincol;                        %Now ensure that the cropped area is square
        height = maxrow-minrow;
        if height>width
            maxcol = maxcol + floor((height-width)/2) + mod(height-width,2);
            mincol = mincol - floor((height-width)/2);
        elseif width>height
            maxrow = maxrow + floor((width-height)/2) + mod(width-height,2);
            minrow = minrow - floor((width-height)/2);
        end
        if minrow<1;      minrow = 1; end
        if maxrow>r_dim;  maxrow = r_dim; end
        if mincol<1;      mincol = 1; end
        if maxcol>c_dim;  maxcol = c_dim; end
        IM_phase = IM_phase(minrow:maxrow, mincol:maxcol);        %Crop the original image to save computation time
    end
        
    [dimx, dimy] = size(IM_phase);
    dx = zeros(dimx,dimy);
    p = unwrap([IM_phase(:,1) IM_phase(:,2)],[],2);
    dx(:,1) = (p(:,2) - IM_phase(:,1))./2;                    %Take the partial derivative of the unwrapped phase in the x-direction for the first column
    p = unwrap([IM_phase(:,dimy-1) IM_phase(:,dimy)],[],2);
    dx(:,dimy) = (p(:,2) - IM_phase(:,dimy-1))./2;            %Take the partial derivative of the unwrapped phase in the x-direction for the last column
    for i=2:dimy-1
        p = unwrap([IM_phase(:,i-1) IM_phase(:,i+1)],[],2);
        dx(:,i) = (p(:,2) - IM_phase(:,i-1))./3;              %Take partial derivative of the unwrapped phase in the x-direction for the remaining columns
    end
    dy = zeros(dimx,dimy);
    q = unwrap([IM_phase(1,:)' IM_phase(2,:)'],[],2);
    dy(1,:) = (q(:,2)' - IM_phase(1,:))./2;                   %Take the partial derivative of the unwrapped phase in the y-direction for the first row
    %p = unwrap([IM_phase(dimx-1,:)' IM_phase(dimx,:)'],[],2);  % unused--appears to be a bug
    q = unwrap([IM_phase(dimx-1,:)' IM_phase(dimx,:)'],[],2);  % Corrected
    dy(dimx,:) = (q(:,2)' - IM_phase(dimx-1,:))./2;           %Take the partial derivative of the unwrapped phase in the y-direction for the last row
    for i=2:dimx-1
        q = unwrap([IM_phase(i-1,:)' IM_phase(i+1,:)'],[],2);
        dy(i,:) = (q(:,2)' - IM_phase(i-1,:))./3;             %Take the partial derivative of the unwrapped phase in the y-direction for the remaining rows
    end
    dx_centre = dx(2:dimx-1, 2:dimy-1);
    dx_left = dx(2:dimx-1,1:dimy-2); 
    dx_right = dx(2:dimx-1,3:dimy);
    dx_above = dx(1:dimx-2,2:dimy-1);
    dx_below = dx(3:dimx,2:dimy-1);
    mean_dx = (dx_centre+dx_left+dx_right+dx_above+dx_below)./5;
    dy_centre = dy(2:dimx-1, 2:dimy-1);
    dy_left = dy(2:dimx-1,1:dimy-2); 
    dy_right = dy(2:dimx-1,3:dimy);
    dy_above = dy(1:dimx-2,2:dimy-1);
    dy_below = dy(3:dimx,2:dimy-1);
    mean_dy = (dy_centre+dy_left+dy_right+dy_above+dy_below)./5;
    stdvarx = sqrt( (dx_left - mean_dx).^2 + (dx_right - mean_dx).^2 + ...
                  (dx_above - mean_dx).^2 + (dx_below - mean_dx).^2 + (dx_centre - mean_dx).^2 ); 
    stdvary = sqrt( (dy_left - mean_dy).^2 + (dy_right - mean_dy).^2 + ...
                  (dy_above - mean_dy).^2 + (dy_below - mean_dy).^2 + (dy_centre - mean_dy).^2 ); 
    derivative_variance = 100*ones(dimx, dimy);                         %Ensure that the border pixels have high derivative variance values
    derivative_variance(2:dimx-1, 2:dimy-1) = stdvarx + stdvary;
    if nargin>=2                                                      %Does the image have to be padded back to the original size?
        [orig_rows, orig_cols] = size(IM_mask);
        temp = 100*ones(orig_rows, orig_cols);
        temp(minrow:maxrow, mincol:maxcol) = derivative_variance;       %Pad the remaining pixels with poor phase quality values
        derivative_variance = temp;
    end
end
%%
function IM_unwrapped = GuidedFloodFill(IM_phase, IM_mag, IM_unwrapped, unwrapped_binary, derivative_variance, adjoin, IM_mask)
    [r_dim, c_dim] = size(IM_phase);
    % Include edge pixels
    while sum(sum(adjoin(:))) ~= 0  %Loop until there are no more adjoining pixels
      adjoining_derivative_variance = derivative_variance.*adjoin + 101.*~adjoin; %Derivative variance values of the adjoining pixels (pad the zero adjoining values with 100)
      min_deriv_var = min(adjoining_derivative_variance(:)); % the minimum derivative variance
      if(min_deriv_var >= 101)  % 101 is an indicator of already assigned variance
        break;  % finished
      end
      [r_adjoin, c_adjoin] = find(adjoining_derivative_variance==min_deriv_var); %Obtain coordinates of the adjoining unwrapped phase pixel w/ the minimum derivative variance
      for(ii = 1:length(r_adjoin))
        r_active = r_adjoin(ii);
        c_active = c_adjoin(ii);
        if(adjoin(r_active,c_active)==1)  % find a valid point
          break;
        end
      end
      if(adjoin(r_active,c_active)==0)  % No valid point was found, 
        break;
      end
      phasev   = nan(1,4);     % Initialize.  Will overwrite for valid pixels
      IM_magv  = nan(1,4);     % Initialize.  Will overwrite for valid pixels
      %First search below for an adjoining unwrapped phase pixel
      if(r_active+1<=r_dim)  % Is this a valid index?
        if unwrapped_binary(r_active+1, c_active)==1
          phase_ref = IM_unwrapped(r_active+1, c_active);       % Obtain the reference unwrapped phase
          % Itoh's Method (suggested by 'Eric' on MATLAB Central to use for a length 2 vector):
          D = IM_phase(r_active, c_active)-phase_ref;
          deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
          phasev(1) = phase_ref + deltap;  % This is the unwrapped phase
          IM_magv(1)= IM_mag(r_active+1, c_active);
        else % unwrapped_binary(r_active+1, c_active)==0
          if(IM_mask(r_active+1, c_active)==1)
            adjoin(r_active+1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
          end
        end
      end
      %Then search above
      if(r_active-1>=1)  % Is this a valid index?
        if unwrapped_binary(r_active-1, c_active)==1
          phase_ref = IM_unwrapped(r_active-1, c_active);                                   %Obtain the reference unwrapped phase
          D = IM_phase(r_active, c_active)-phase_ref;
          deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
          phasev(2) = phase_ref + deltap;  % This is the unwrapped phase
          IM_magv(2)= IM_mag(r_active-1, c_active);
        else % unwrapped_binary(r_active-1, c_active)==0
          if(IM_mask(r_active-1, c_active)==1)
            adjoin(r_active-1, c_active) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
          end
        end
      end
      %Then search on the right
      if(c_active+1<=c_dim)  % Is this a valid index?
        if unwrapped_binary(r_active, c_active+1)==1
          phase_ref = IM_unwrapped(r_active, c_active+1);                                   %Obtain the reference unwrapped phase
          D = IM_phase(r_active, c_active)-phase_ref;
          deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
          phasev(3) = phase_ref + deltap;  % This is the unwrapped phase
          IM_magv(3)= IM_mag(r_active, c_active+1);
        else % unwrapped_binary(r_active, c_active+1)==0
          if(IM_mask(r_active, c_active+1)==1)
            adjoin(r_active, c_active+1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
          end
        end
      end
      %Finally search on the left
      if(c_active-1>=1)  % Is this a valid index?
        if unwrapped_binary(r_active, c_active-1)==1
          phase_ref = IM_unwrapped(r_active, c_active-1);                                   %Obtain the reference unwrapped phase
          D = IM_phase(r_active, c_active)-phase_ref;
          deltap = atan2(sin(D),cos(D));   % Make it modulo +/-pi
          phasev(4) = phase_ref + deltap;  % This is the unwrapped phase
          IM_magv(4)= IM_mag(r_active, c_active-1);
        else % unwrapped_binary(r_active, c_active-1)==0
          if(IM_mask(r_active, c_active-1)==1)
            adjoin(r_active, c_active-1) = 1;  % Put the elgible, still-wrapped neighbors of this pixels in the adjoin set
          end
        end
      end
      idx_del = ~isnan(phasev);
      if(any(idx_del)) % Any valid adjoining pixels?
           % Use the strongest neighbor
        IM_max  = max(IM_magv(idx_del));
        idx_max = find((IM_magv >= 0.99*IM_max) & (idx_del==1));
        IM_unwrapped(r_active, c_active) = phasev(idx_max(1));  % Use the first, if there is a tie
        unwrapped_binary(r_active, c_active) = 1;      %Mark the pixel as unwrapped
        adjoin(r_active, c_active) = 0;              %Remove it from the list of adjoining pixels
      else  % no valid adjoining pixels found
        adjoin(r_active,c_active) = 0;  %Remove the current active pixel from the adjoin list
        continue;
      end
      %end
    end % while sum(sum(adjoin(2:r_dim-1,2:c_dim-1))) ~= 0  %Loop until there are no more adjoining pixels
    disp(['All of the valid interior pixels have been calculated']);
end


