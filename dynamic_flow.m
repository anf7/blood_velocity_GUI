function dynamic_flow(framespacing,sequencelength)

scale_factor = 0.322; %brain study

[filename,pathname] = uigetfile('*.mat','Select .mat files to process...','MultiSelect','off');
mkdir(pathname,'Dynamic Flow')


load(strcat(pathname,filename));

z = size(imagestack,1);

startframe = 1;
endframe = sequencelength;

while endframe < z 
    
    stackpart = imagestack(startframe:endframe,:,:);
    
    mean_img = squeeze(mean(stackpart,1)); 
    
    [zz r c] = size(stackpart);
    zimg = zeros(zz,r,c);
    for n = 1:zz
        zimg(n,:,:) = (squeeze(stackpart(n,:,:)) - mean_img)./mean_img;
    end
    
    fluctimg = squeeze(sum((zimg.^2),1)); 
    clear zimg

    [movement_matrix togetherimg] = flowmaps(0,0,[],[],scale_factor,0,0,stackpart);
    newfilename = strcat('dynamicflow_frame',num2str(startframe,'%04u'),'to',num2str(endframe,'%04u'),'.mat');
    
    save(strcat(pathname,'Dynamic Flow\',newfilename),'movement_matrix','togetherimg','fluctimg')
    
    startframe = startframe + framespacing;
    endframe = endframe + framespacing;
end