clear all;
clc;
%% Frame detection
videoObj = Videoreader('data/Walk1.mpg');
hight = videoObj.Height;
width = videoObj.Width;
nFrame = 0 ;
isGood = true ;

figure;
pause on;
while isGood
    nFrame = nFrame + 1 ;
    try
        thisFrame = read(videoObj, nFrame) ;
%         grayFrame = rgb2gray(thisFrame);    
    catch
        isGood = false ;
        nFrame = nFrame - 1 ; % Don't include in the count if fails
    end
    temp_r=thisFrame(:,:,1);% r value pixels of current frame 
    temp_g=thisFrame(:,:,2);% g value pixels of current frame 
    temp_b=thisFrame(:,:,3);% b value pixels of current frame 
    allFrames_r(:,nFrame)=temp_r(:);%converting matrix to vector 
    allFrames_g(:,nFrame)=temp_g(:);
    allFrames_b(:,nFrame)=temp_b(:);
   % image of some frames
  %{
    if(nFrame == 100)||(nFrame == 300)||(nFrame == 400)
        figure,imshow(thisFrame);
        if(nFrame == 100)
            title('Frame 100');
        elseif(nFrame == 300)
            title('Frame 300');
        else
            title('Frame 400');
        end
    end
   
   %}
   
subplot(1,2,1);
imshow(thisFrame);
pause(0.08);
end

%% Activity Mask Algorithm
N = nFrame-6;
for k=1 : N
    %inter-frame intensity difference of all pixels
   d_r =double(allFrames_r(:,k+1)-allFrames_r(:,k));
   d_g =double(allFrames_g(:,k+1)-allFrames_g(:,k));
   d_b =double(allFrames_b(:,k+1)-allFrames_b(:,k));
   % combined inter-frame difference of all pixels
   d_k=sqrt(d_r.^2+d_g.^2+d_b.^2);
   D(:,k)=(d_k);% matrix of all N-1 frame difference
  end
    D_pixel_t = D';
    kD_i= kurtosis(D_pixel_t);%4th order mean of each pixel along frames
    kD_i(isnan(kD_i))=0;
    kD_mean = mean(kD_i);%mean value of kurtosis of each pixel
    kD_thsld = kD_mean*1.5;%threslod is the 15% of KD_mean
    kD_mat =  reshape(kD_i',[hight width]);%vector to matrix conversion
    Activity_Frame = zeros(hight,width); 
    %Binary Activity Mask
            for i=1:hight
                for j=1:width
                    if(kD_mat(i,j) >=kD_thsld)
                        Activity_Frame(i,j) = 1;
                    else
                        Activity_Frame(i,j) = 0;
                    end
                end
            end
         
 subplot(1,2,2);
 imshow((Activity_Frame));
 title('Activity mask');
           
%% Sequential likelyhood change testing
 K = hight*width;%No. of pixels
 % Finding Active No. of Pixels
 Act_N=0;
 for k=1:K
      if kD_i(k) >= kD_thsld
         Act_N=Act_N+1;
      end
 end
 % Finding Intensity difference value of active pixels
 D_1_k = zeros(N,Act_N);
 i=1;
  for k=1:K
      if kD_i(k) >= kD_thsld
         D_1_k(:,i) = D_pixel_t(:,k);
         i=i+1;
      end
  end
 wFrames = uint8(0.1*N);% initial 10% frames
  D_1_w0= zeros(wFrames,Act_N);
  for j=1:wFrames
      D_1_w0(j,:) = D_1_k(j,:);% intensity diff. value of pixels in first 10% frame
  end
  mu_0 = mean(mean(D_1_w0));%mean of mean
  var_0 = mean(var(D_1_w0));%mean of variance
  
  kFrames = 200;% upto frame
  hFrames = 150;%past frames
  H1_Frames = (kFrames-hFrames) : kFrames;%frames in motion
  len = length(H1_Frames);
  D_1_w1 = zeros(len,Act_N);
  for l= 1:len
      D_1_w1(l,:)= D_1_k(l,:);% intensity diff. value of pixels in motion frames
  end
   mu_1 = mean(mean(D_1_w1));%mean of mean
   var_1= mean(var(D_1_w1));%mean of variance
   
   % Log-likelyhood test
   num_pixels = length(D_1_k);
   uptoFrame = kFrames;
%    sum=zeros(1,uptoFrame);
  l_k =zeros(uptoFrame,num_pixels);
 for p=1: num_pixels
          sum =0;
     for n=1:uptoFrame
%           if(var_1(p)~=0 && var_0(p)~=0 )
        temp= ( -((D_1_k(n,p) - mu_1)^2/(2*var_1))+...
              ((D_1_k(n,p) - mu_0)^2/(2*var_0)));
%           temp(isnan(temp))=0;
          sum=sum+temp;
          l_k(n,p) = log(var_0/var_1)+sum;
%           end
     end
 end
  
g_k=gradient(l_k);%gradient of LLRT

%Sequential Change Detection
mean_lk =mean(mean(l_k));
pi_l = 1.25*mean_lk;
pi_u=0.75*mean_lk;
for p=1:num_pixels
    for n=1: uptoFrame
       if(l_k(n,p)>pi_u)
          SCD(n,p) =1;
       elseif (l_k(n,p)<pi_l)
          SCD(n,p) = 0;
       else
        SCD(n,p)=2;
       end
    end
end
%Results
r_0_pel = 449;% particular pixel number
m=1;

for n=1:(uptoFrame-1)
    if(SCD(n,r_0_pel)~=SCD(n+1,r_0_pel))
        changeFrame(m) =n+1
        m=m+1;
    end
end
% inter-frame luminace variation
i=1;
for n=1:uptoFrame
  lum_r0(i,r_0_pel) =D_1_k(n,r_0_pel);
  i=i+1;     
end
% plots
y=1:uptoFrame;
figure,plot(y,lum_r0(:,r_0_pel));
xlabel('Frame numbers');
ylabel('interframe illumination variation');
title('D_1 for pixel r0');
grid on
figure,plot(y,l_k(:,r_0_pel));
xlabel('Frame numbers');
ylabel('log-likelyhood ratio for active pixel');
title('Log-likelyhood ratio over all frames');
grid on
figure,plot(y,g_k(:,r_0_pel));
xlabel('Frame numbers');
ylabel('gradiant of log-likelyhood ratio for active pixel');
title('Gradiant of log-likelyhood ratio');
grid on
gtext('change');
 
