clear all;
close all;
clc;
original_image= imread('E:\Final Year\image and video coding\Parrots-680x680.jpg');%change this to check the code
figure
imshow(original_image) ; title('Original Image')

cropped_image= imcrop(original_image,[200 16 15 15]); % E number 404
figure ;
imshow (cropped_image);title('Cropped Image');
cropped_gray=rgb2gray(cropped_image); 
figure;
imshow(cropped_gray); title('Cropped GRAY Scale Image')
full_gray=rgb2gray(original_image);
full_gray2=full_gray;
figure;imshow(full_gray); title('Original GRAY Scale Image')


%quntization
thresh = multithresh(cropped_gray,7);
valuesMax = [thresh max(cropped_gray(:))];
[quant8_I, index] = imquantize(cropped_gray, thresh, valuesMax);
figure; imshow(quant8_I,[]); title('Quantized Cropped Image')
% probabilities
[g,~,Intensity_val] = grp2idx(quant8_I(:));
Frequency  = accumarray(g,1);




[Intensity_val Frequency];
probability=Frequency./256; 
prob2 =Frequency./256;
T = table(Intensity_val,Frequency,probability); 
T(1:length(Frequency),:);


%Algorithm for Huffman

symbol=cell(8,2); 
codebook=cell(8,2); 

symbol2= ["a";"b";"c";"d";"e";"f";"g";"h"]; % use these symbols to track

%initial placement
for m=1:8

symbol{m}=cellstr(symbol2(m)); 
codebook{m,1}= cellstr(symbol2(m));
end

 for a = 8:-1:2
     
     %probability sorting
     for i = 1:length(probability)
         for j = 1:length(probability)
            if (probability(i)>probability(j))
                  num=probability(i);
                  probability(i)=probability(j);
                  probability(j)=num;
                  nnn=symbol{i};
                  symbol{i}=symbol{j}; % symbol sorting the same way probability did
                  symbol{j}=nnn;
                  ttt=symbol2(i);
                  symbol2(i)=symbol2(j);
                  symbol2(j)=ttt;
            end
         end
     end
  sort_prob = probability; 
  probability(a-1)= probability(a-1)+probability(a); %huffmann iterative addition of probabilities
  probability(a)=[];% make last raw empty 
  
  sort_sym = symbol2; 
  symbol2(a-1)= strcat(symbol2(a-1),symbol2 (a));% concatination of symbols
  symbol2(a)=[];
  
  sorted= symbol;
  symbol{a-1}= horzcat(symbol{a-1},symbol{a}); %combine the vectors to a matrix
  symbol{a,1}=[];
  

  
 %bits assignig
 for p = 1:1:length(sorted{a,1})
     for r=1:8
     tf1 = ismember(sorted{a,1}{1,p},codebook{r,1}); % check equality  
        if (tf1==1)
          codebook{r,2}=[codebook{r,2} 1]; % assign 1f matched
        end 
     end
 end
 
  for q = 1:1:length(sorted{a-1,1})
     for t=1:8
     tf2 = ismember(sorted{a-1,1}{1,q},codebook{t,1}); % check equality
        if (tf2==1)
          codebook{t,2}=[codebook{t,2} 0]; % assign 1 if matched
        end 
     end
 end
  
 end
for y = 1:8
   codebook{y,2}=fliplr(codebook{y,2});% make asscending
   codebook{y,1}=[Intensity_val(y)]; % codebook assign
    
end

%Encoding code
cropped_gray3=cropped_gray;
huffman=cell(16,16);
for q=1:size(cropped_gray3,1)
    for s=1:size(cropped_gray3,2)
        for k=1:1:8
            if isequal(codebook{k,1},quant8_I(q,s))
                huffman(q,s)=codebook(k,2);
                break
            end
        end
    end
end

%Decoding algorithm from using codebook
huffman2=huffman(:)%making an array

d = 1;

for v = 1:1:16
    for w = 1:1:16
        for e = 1:1:length(huffman{v,w})
            encode_array(d,1) = huffman{v,w}(1,e);
            d = d+1;
        end
    end
end

decode=zeros(256,1);
d=1;m=1;
while d<length(encode_array)
    if encode_array(d,1)==0
        if encode_array(d+1,1)==0
            decode(m,1)=codebook{5,1};
            m=m+1;
            d=d+2;
        else
            decode(m,1)=codebook{4,1};
            m=m+1;
            d=d+2;
        end
    else if encode_array(d+1,1)==0%2nd
            if encode_array(d+2,1)==1%3rd
                if encode_array(d+3,1)==1 %4th
                    if encode_array(d+4,1)==1%5th
                        
                        decode(m,1)=codebook{8,1};
                        m=m+1;
                        d=d+5;
                    
                        else
                        decode(m,1)=codebook{7,1};
                        m=m+1;
                        d=d+5;
                        end
                else
                    decode(m,1)=codebook{1,1};
                    m=m+1;
                    d=d+4;
                end
            else
                decode(m,1)=codebook{3,1};
                m=m+1;d=d+3;
                
                end
        else
                if encode_array(d+2)==1;
                decode(m,1)=codebook{2,1};
                m=m+1;
                d=d+3;
                
            
        
            else decode(m,1)=codebook{6,1};
                m=m+1;
                d=d+3;
            end
        end
    end
end
v=1;
for m=1:1:16
    for n=1:1:16
        img(m,n)=decode(v,1);
        v=v+1;
    end
end
        
        
    
    
    



%Encoding for cropped using inbuilt
column=quant8_I(:);%making a column vector 
encode_cropped=huffmanenco(column,codebook);

croppedd = huffmandeco(encode_cropped,codebook);
%decoding for cropped using inbuilt
decode_cropped = vec2mat(croppedd,16)';
figure; imshow(decode_cropped,[]), title('Decoded GRAY Scale Cropped ')



%Encoding for Original grayscale
kkk=cropped_gray(:);
maximum=max(kkk);%max of codebooked image part
kkk2=cropped_gray(:);
minimum=min(kkk2);%min of it
%margin creation
for i=1:1:length(full_gray)
    for j=1:1:length(full_gray)
        if full_gray(i,j)>maximum
            full_gray(i,j)=maximum;
        else if full_gray(i,j)<minimum
                full_gray(i,j)=minimum;
            else
                full_gray(i,j)=full_gray(i,j);
            end
        end
    end
end



ttt = multithresh(full_gray,7);
valuesMax2 = [thresh max(full_gray(:))];
[quant8_I2, index] = imquantize(full_gray, ttt, valuesMax2);
figure; imshow(quant8_I2 ,[]),title('GRAY Scale Quantized')

quant8_I2=quant8_I2(:);
encode2=huffmanenco(quant8_I2,codebook);
decode2 = huffmandeco(encode2,codebook);

return_image = vec2mat(decode2,680)';
figure; imshow(return_image,[]), title('Decoded GRAY Scale ')


%Entropy calculation
Entropy1= entropy(original_image)
Entropy=entropy(full_gray)
Entropy2=entropy(full_gray2)
Entropy3=entropy(cropped_image)
Entropy4=entropy(cropped_gray)
Entropy5=entropy(return_image)
Entropy6=entropy(decode_cropped)

%Compression ratio valuating file size is obtained
%To check the code change the below paths
info_cropped_gray = imfinfo('E:\Final Year\image and video coding\lab1\croppedimage.jpg');
info_cropped_gray_decoded = imfinfo('E:\Final Year\image and video coding\lab1\decoded_gray_copied.jpg');
info_graydecoded = imfinfo('E:\Final Year\image and video coding\lab1\gray_decoded full.jpg');
info_grayscale = imfinfo('E:\Final Year\image and video coding\lab1\gray_fll.jpg');

size_cropped_gray = info_cropped_gray.FileSize;
size_cropped_gray_decoded= info_cropped_gray_decoded.FileSize;
size_grayscaleD = info_graydecoded.FileSize;
size_grayscale = info_grayscale.FileSize;

%Compression Ratio Value calculation
CR1 =(size_grayscale/size_grayscaleD) 
CR2=(size_cropped_gray/size_cropped_gray_decoded)

%PSNR value calculations
PSNR2 = psnr(decode_cropped,cropped_gray) 
PSNR3 = psnr(return_image,full_gray) 

