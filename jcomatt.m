function [I out]=jcomatt(I,Q);

I=imresize(I,[512 512]);

input_image = im2double(I);

dct_8x8_image_of{1} = image_8x8_block_dct( input_image );

mean_matrix_8x8 = zeros( 8,8 );

  for m = 0:63
      for n = 0:63
         mean_matrix_8x8 = mean_matrix_8x8 +abs( dct_8x8_image_of{1}(m*8+[1:8],n*8+[1:8]) ).^2;
      end
  end
  
original_image = im2double(I);

coef_selection_matrix = zeros(8,8);

compressed_set = [Q];

for number_of_coefficient = 1:64
    % find the most energetic coefficient from the mean_matrix
    [y,x] = find(mean_matrix_8x8==max(max(mean_matrix_8x8)));    
    % select if for the compressed image
    coef_selection_matrix(y,x) = 1;    
    % replicate the selection matrix for all the parts of the dct transform
    % (remember that the DCT transform creates a set of 8x8 matrices, where
    %  in each matrix I need to choose the coefficients defined by the 
    %  <<coef_selection_matrix>> matrix )
    selection_matrix = repmat( coef_selection_matrix,64,64);    
    % set it as zero in the mean_matrix, so that in the next loop, we will
    % choose the "next-most-energetic" coefficient
    mean_matrix_8x8(y,x) = 0;    
    % choose the most energetic coefficients from the original image
    % (total of <<number_of_coefficient>> coefficients for this run in the loop)
%     size(I)
%     size(selection_matrix)
    compressed_image = image_8x8_block_dct(I) .* selection_matrix;    
    % restore the compressed image from the given set of coeficients
    restored_image = image_8x8_block_inv_dct( compressed_image );    
    % calculate the snr of this image (based on the original image)  
    if ~isempty(find(number_of_coefficient==compressed_set))        
         out=restored_image;
         %imshow(out)
   end
end




function out = pdip_dct2( in )

% get input matrix size
N = size(in,1);
% build the matrix
n = 0:N-1;
for k = 0:N-1
   if (k>0)
      C(k+1,n+1) = cos(pi*(2*n+1)*k/2/N)/sqrt(N)*sqrt(2);
   else
      C(k+1,n+1) = cos(pi*(2*n+1)*k/2/N)/sqrt(N);
   end   
end

out = C*in*(C');

% ---------------------------------------------------------------------------------
% pdip_inv_dct2 - implementation of an inverse 2 Dimensional DCT
%
% assumption: input matrix is a square matrix !
% ---------------------------------------------------------------------------------
function out = pdip_inv_dct2( in )

% get input matrix size
N = size(in,1);

% build the matrix
n = 0:N-1;
for k = 0:N-1
   if (k>0)
      C(k+1,n+1) = cos(pi*(2*n+1)*k/2/N)/sqrt(N)*sqrt(2);
   else
      C(k+1,n+1) = cos(pi*(2*n+1)*k/2/N)/sqrt(N);
   end   
end

out = (C')*in*C;

% ---------------------------------------------------------------------------------
% image_8x8_block_dct - perform a block DCT for an image
% ---------------------------------------------------------------------------------
function transform_image = image_8x8_block_dct( input_image )

transform_image = zeros( size( input_image,1 ),size( input_image,2 ) );
for m = 0:63
    for n = 0:63
        transform_image( m*8+[1:8],n*8+[1:8] ) = ...
            pdip_dct2( input_image( m*8+[1:8],n*8+[1:8] ) );
    end
end


% ---------------------------------------------------------------------------------
% image_8x8_block_inv_dct - perform a block inverse DCT for an image
% ---------------------------------------------------------------------------------
function restored_image = image_8x8_block_inv_dct( transform_image )

restored_image = zeros( size( transform_image,1 ),size( transform_image,2 ) );
for m = 0:63
    for n = 0:63
        restored_image( m*8+[1:8],n*8+[1:8] ) = ...
            pdip_inv_dct2( transform_image( m*8+[1:8],n*8+[1:8] ) );
    end
end


% ---------------------------------------------------------------------------------
% calc_snr - calculates the snr of a figure being compressed
%
% assumption: SNR calculation is done in the following manner:
%             the deviation from the original image is considered 
%             to be the noise therefore:
%
%                   noise = original_image - compressed_image
%
%             the SNR is defined as:  
%
%                   SNR = energy_of_image/energy_of_noise
%
%             which yields: 
%
%                   SNR = energy_of_image/((original_image-compressed_image)^2)
% ---------------------------------------------------------------------------------
function SNR = calc_snr( original_image,noisy_image )

original_image_energy = sum( original_image(:).^2 );
noise_energy = sum( (original_image(:)-noisy_image(:)).^2 );
SNR = original_image_energy/noise_energy;
