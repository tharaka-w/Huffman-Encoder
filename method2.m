clear all 
close all
clc;
original_image= imread('E:\Final Year\image and video coding\Parrots-680x680.jpg');
figure
imshow(original_image);title('Original Image');
cropped_image=imcrop(original_image,[200 16 15 15]);
figure;
imshow(cropped_image); title('Cropped Image');

cropped_gray=rgb2gray(cropped_image);
figure;
imshow(cropped_gray);title('Cropped GRAY Scale Image');

full_gray=rgb2gray(original_image);
figure;imshow(full_gray);title('Original GRAY Scale Image');

%Quantizing
[Width,Length]=size(cropped_image);
thresh=multithresh(cropped_gray,7);%to 8 levels
valuesMax=[thresh max(cropped_gray(:))];%??????????
[quant8_I,index]=imquantize(cropped_gray,thresh,valuesMax);%?????
figure;imshow(quant8_I,[]);title('Quantized Cropped Image')%?????

%find probabilities
[g,~,Intensity_val]=grp2idx(quant8_I(:));%(:)reshapes to a column vector
Frequency=accumarray(g,1);

[Intensity_val Frequency];
probability=double(Frequency)./256;
prob2=double(Frequency./256);

table3=table(Intensity_val, Frequency ,probability);% Edited this f0r 2019 version

T=[Intensity_val Frequency probability]%table form element count probability
T(1:length(Frequency),:);
% 
% p=length(Frequency)
% 
% w=min(probabilities);
table4=sortrows(table3,3);
%T[order(T$V1),]
B = sortrows(T,2)
% for i=1:1:p
%     for w>probabilities(i)
%         w2=probability(i)+w;
%     end
%     if w2=1
%         
B(:,2);
freq_value_ascend=(double(B(:,2)))./256;
f2=freq_value_ascend;

num=zeros(length(freq_value_ascend));
num(1)=  freq_value_ascend(1)+freq_value_ascend(2);
layer=zeros(length(freq_value_ascend));
layer(1)=num(1);
x=0;
for i=1:1:length(freq_value_ascend)
    if i<= length(freq_value_ascend)-3
        num(i+1)=freq_value_ascend(i+2)+freq_value_ascend(i+3);
         if num(i+1)>=(layer(i)+freq_value_ascend(i+2))
             layer(i+1)=layer(i)+freq_value_ascend(i+2);
         else
            count=0;
          
            layer_found(i)=i+2+x;
            x=x+1;
            freq_value_ascend(i+2)=freq_value_ascend(i+2)+freq_value_ascend(i+3);
            for j=(i+2):1:(length(freq_value_ascend))
                if (j+1)==length(freq_value_ascend)
                    break 
                else
                    
                    count=count+1;
                    freq_value_ascend(j+1)=freq_value_ascend(j+2)
                end
            end
                freq_value_ascend((length(freq_value_ascend)-1),:)=[];
                layer(i+1)=layer(i)+freq_value_ascend(i+2);
         end
    end
end



Depth=0;

b=length(probability);
while b>0
    code=repmat('1',[1 Depth]); 
    code2=string(code+string(0));
    if b~=1
        for r=1:1:length(layer_found)
            measure=0;
            if (layer_found(r)>0)&(b==(layer_found(r)+1))
                H(b,:)=string(code2+string(0));
                H(b-1,:)=string(code2+string(1));
                b=b-1;
                measure=measure+1;
                break
            
                
                %H(b,:)=string(code+string(0));
            end
            
        end
        if measure==0;
            H(b,:)=string(code+string(0));
        end
    else
        H(b,:)=string(code);
    
    end
   
    Depth=Depth+1;
    b=b-1;
end





Table=table(table4,H)
values=B(:,1)

%Encoding
for q=1:size(cropped_gray,1)
    for s=1:size(cropped_gray,2)
        for k=1:1:length(f2)
            if quant8_I(q,s)==values(k,1)
                cropped_gray2(q,s)=H(k,1);
                break
            end
        end
    end
end
        
% %Decoding
% for q=1:size(cropped_gray2,1)
%     for s=1:size(cropped_gray2,2)
%         for k=1:1:length(f2)
%             if cropped_gray2(q,s)==H(k,1)
%                 cropped_gray2_decoded(q,s)=values(k,1);
%                 break
%             end
%         end
%     end
% end
% figure;imshow(cropped_gray2_decoded,[]);title('Decoded');    
%     
    












