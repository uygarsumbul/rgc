function chatstruct=processchatNN(chat)



load NNParams.mat;

smallx=floor(size(chat,2)/10);
smally=floor(size(chat,1)/10);

smallchatcomplete=zeros(smally,smallx,size(chat,3));

for i=1:size(chat,3)
    smallchatcomplete(:,:,i)=imresize(chat(:,:,i),[smally smallx],'bilinear');%*(1+0.8*i/size(a.chat,3));
end

smallchat=[];

for i=2:2:size(chat,3)
smallchat(:,:,i/2)=smallchatcomplete(:,:,i);
end
smallchat=permute(smallchat,[1 3 2]);

smallchat=double(smallchat)/max(smallchat(:));

output=NNpass(W,B,smallchat,num_of_maps_in_layer);
CHAT=output{3}>0.5*max(output{3}(:));
CHAT=bwareaopen(CHAT,1000);

% if there are more than two components
l=bwlabeln(CHAT);
if max(l(:))>2
    for i=1:max(l(:))
        avg(i)=mean(output{3}(l==i));
    end
    [v ind]=min(avg);
    CHAT(l==ind)=0;
end

% fix small holes
smallchat=permute(smallchat,[1 3 2]);
CHAT=permute(CHAT,[1 3 2]);

outchat=zeros(size(smallchat));
for i=7:size(smallchat,3)-6
    outchat(7:end-6,7:end-6,i)=CHAT(:,:,i-6);
end


%CHAT=permute(CHAT,[1 3 2]);

histfg=smallchat(outchat(:)>0);
a=sort(histfg);
thres=a(round(length(histfg)*0.98));

% get rid of bright somas

enhanced=smallchat;
enhanced=permute(enhanced,[1 3 2]);

for i=1:size(enhanced,3)
    enhanced(:,:,i)=medfilt2(enhanced(:,:,i),[2 10],'symmetric');
%    enhanced(:,:,i)=filter2(ones(3,20)/15,enhanced(:,:,i),'same');%-filter2(ones(3,5)/5,enhanced(:,:,i),'same');
end
enhanced=permute(enhanced,[1 3 2]);

enhanced=permute(enhanced,[3 2 1]);
for i=1:size(enhanced,3)
    enhanced(:,:,i)=medfilt2(enhanced(:,:,i),[10 2],'symmetric');
%    enhanced(:,:,i)=filter2(ones(3,20)/15,enhanced(:,:,i),'same');%-filter2(ones(3,5)/5,enhanced(:,:,i),'same');
end
enhanced=permute(enhanced,[3 2 1]);


for i=1:size(enhanced,3)
   enhanced(:,:,i)=medfilt2(enhanced(:,:,i),[2 2])-medfilt2(enhanced(:,:,i),[6 6],'symmetric');
%    filter2(ones(3,3)/9,smallchat(:,:,i),'same')-filter2(ones(7,7)/225,smallchat(:,:,i),'same'); 
end

enhancedfull=[];
for i=1:size(smallchatcomplete,3)
   enhancedfull(:,:,i)=enhanced(:,:,min(ceil(i/2),size(enhanced,3)));
%    filter2(ones(3,3)/9,smallchat(:,:,i),'same')-filter2(ones(7,7)/225,smallchat(:,:,i),'same'); 
end




% pick out the centers of the cells
l=enhancedfull>max(enhanced(outchat(:)>0));
%lo=enhanced>max(enhanced(
lbig=imdilate(l,strel('arbitrary',ones([5 5 5])));

outsmally=round(smally/5);
outsmallx=round(smallx/5);
smallchatcomplete(smallchatcomplete>thres*max(smallchatcomplete(:)))=thres*max(smallchatcomplete(:));
smallchatNN=[];
for i=1:1:size(chat,3)
    smallchatNN(:,:,i)=imresize(smallchatcomplete(:,:,i),[outsmally outsmallx],'bilinear');%*(1+0.8*i/size(a.chat,3));
end

somas=zeros(size(smallchatcomplete));
props=regionprops(bwlabeln(l));
% for i=1:length(props)
%     props(i).Centroid=round(props(i).Centroid);
%     xextend=max(props(i).Centroid(1)-4,1):min(props(i).Centroid(1)+4,size(enhanced,2));
%     yextend=max(props(i).Centroid(2)-4,1):min(props(i).Centroid(2)+4,size(enhanced,1));
%     zextend=max(props(i).Centroid(3)-10,1):min(props(i).Centroid(3)+10,size(enhancedfull,3));
%     somas(yextend,xextend,zextend)=1;
% end


somas(lbig)=1;

smallchatsoma=[];
for i=1:1:size(chat,3)
    smallchatsoma(:,:,i)=imresize(somas(:,:,i),[outsmally outsmallx],'bilinear');%*(1+0.8*i/size(a.chat,3));
end

                % bring back chat pixels
% suppress very bright pixels

smallchat=smallchatNN;


outchat=zeros(size(smallchatNN));
for i=13:size(chat,3)-12
    chatslice=zeros(smally,smallx);
    chatslice(7:end-6,7:end-6)=CHAT(:,:,min(ceil((i-12)/2),size(CHAT,3)));
    outchat(:,:,i)=imresize(chatslice,[outsmally outsmallx],'bilinear');
end
%outchat=imdilate(outchat>0,strel(ones(5,5,1)));


% outchat=permute(outchat,[1 3 2]);
% for i=1:size(outchat,3)
%     outchat(:,:,i)=imdilate(outchat(:,:,i),strel('rectangle',[3 8]));
% end
% outchat=permute(outchat,[1 3 2]);

% extrapolate to the edges
outchat(1,:,:)=reshape(imdilate(squeeze(outchat(2,:,:)),strel('square',4)),[1 size(outchat,2) size(outchat,3)]);
outchat(outsmally,:,:)=reshape(imdilate(squeeze(outchat(outsmally-1,:,:)),strel('square',4)),[1 size(outchat,2) size(outchat,3)]);
outchat(:,1,:)=reshape(imdilate(squeeze(outchat(:,2,:)),strel('square',4)),[size(outchat,1) 1 size(outchat,3)]);
outchat(:,outsmallx,:)=reshape(imdilate(squeeze(outchat(:,outsmallx-1,:)),strel('square',4)),[size(outchat,1) 1 size(outchat,3)]);

outchatthin=outchat;


outchat=permute(outchat,[1 3 2]);
for i=1:size(outchat,3)
    outchat(:,:,i)=imdilate(outchat(:,:,i),strel('rectangle',[4 4]));
end
outchat=permute(outchat,[1 3 2])>0;

smallchatNN(outchat==0)=0;
smallchat1=smallchatNN;

smallchatsoma(outchat>0)=0;

% smallchat1=smallchat;


% 
% for i=1:size(smallchat,3)
%     smallchat(:,:,i)=  filter2(ones(10,10)/100,smallchat(:,:,i),'same');
% end

% smallchat1=[];
% for i=1:size(smallchat,3)
%     smallchat1(:,:,i)=imresize(smallchat(:,:,i),[20 20],'bilinear'); %'box'
% end

% smalldapi1=[];
% for i=1:size(smalldapi,3)
%     smalldapi1(:,:,i)=imresize(smalldapi(:,:,i),[20 20],'bilinear'); %'box'
% end


%  zprofdapi=squeeze(mean(mean(dapi,1),2));
%  zmin=find(imregionalmin(zprofdapi));
%  % find the one closest to the center
%  [val zcentermin]=min(abs(zmin-size(dapi,3)/2));
%  zcentermin=zmin(zcentermin);


chatzmin=[];
chatzmax=[];
% get the chat bands
for i=1:size(smallchatNN,1)
    for j=1:size(smallchatNN,2)
        
        % does outchat has two nice peaks indicating success of NN
        % algorithm?
        
        zprof=smallchat(i,j,:);
        zprof=smooth(zprof,5);
        [maxtab mintab]=peakdet(zprof,0.5);

        peaks=zeros(size(zprof));
        for k=1:size(maxtab,1)
           peaks(maxtab(k,1))=maxtab(k,2);
        end
        peaks(outchat(i,j,:)==0)=0;
        % are the peaks 10 slices apart?
        if sum(peaks>0)>=2 
            a=sort(peaks,'descend');
            ind=find(peaks>=a(2));
            if (ind(2)-ind(1))>0
                chatzmin(i,j)=ind(1);
                chatzmax(i,j)=ind(2);
                continue;
            end
        end
        
        
        % suppress really broad peaks
        zprof1=zprof-smooth(zprof,15);
        [maxtab mintab]=peakdet(zprof1,0.5);

        peaks=zeros(size(zprof1));
        for k=1:size(maxtab,1)
           peaks(maxtab(k,1))=maxtab(k,2);
        end
        peaks(smallchatsoma(i,j,:)>0)=0;
        
      %  maxtab=maxtab(maxtab(:,2)>max(zprof)*0.5
        
        
            
        
        
%         zprof=outchatthin(i,j,:);
%         [maxtab mintab]=peakdet(zprof,0.5);
%         
%         
%         
%         
%         zprof=smallchat1(i,j,:);
%         zprof=zprof-min(zprof);
% %         zprofdapi=smalldapi1(i,j,:);
% %         zprofdapi=zprofdapi-min(zprofdapi);
% %         zprof(zprofdapi(1:zcentermin)>max(zprofdapi(1:zcentermin))*0.5)=0;
% %         zprof(zcentermin+find(zprofdapi(zcentermin+1:end)>max(zprofdapi(zcentermin+1:end))*0.5))=0;
% %         %        zprof=medfilt2(zprof(:),[10 1]);
%      
%         zprof(zprof<max(zprof)*0.3)=0;
        [maxtab mintab]=peakdet(peaks,0.5);
        zmax=maxtab(:,1);
        % get rid of peaks near boundary
        zmax(zmax<5)=[];
        zmax(zmax>length(zprof)-5)=[];
        if length(zmax)==0
            chatzmin(i,j)=0;
            chatzmax(i,j)=1;
        else
            if length(zmax)==1
                chatzmin(i,j)=zmax-8;
                chatzmax(i,j)=zmax+8;
            else

                if length(zmax)>2
%                     i
%                     j
                end

                % two highest peaks
                [peaks ind]=sort(zprof1(zmax),'descend');
                chatzmin(i,j)=zmax(min(ind(1:2)));
                chatzmax(i,j)=zmax(max(ind(1:2)));
            end
        end
    end
end

%rehabitate the problem spots
chatzminfilt=medfilt2(chatzmin,'symmetric');
ind=abs(chatzmin-chatzminfilt)>5;
chatzmin(ind)=chatzminfilt(ind);

chatzmaxfilt=medfilt2(chatzmax,'symmetric');
ind=abs(chatzmax-chatzmaxfilt)>5;
chatzmax(ind)=chatzmaxfilt(ind);

chatzminfilt=medfilt2(chatzmin,'symmetric',[5 5]);
ind=abs(chatzmin-chatzminfilt)>10;
chatzmin(ind)=chatzminfilt(ind);

chatzmaxfilt=medfilt2(chatzmax,'symmetric',[5 5]);
ind=abs(chatzmax-chatzmaxfilt)>10;
chatzmax(ind)=chatzmaxfilt(ind);



% chatzmin= size(chat,3)-chatzmin+1;
% chatzmax= size(chat,3)-chatzmax+1;
stepx=round(size(chat,2)/(1*size(chatzmin,2)));
stepy=round(size(chat,1)/(1*size(chatzmin,1)));


chatx = round(stepx/2):stepx:size(chat,2);
chaty = round(stepy/2):stepy:size(chat,1);
[chatX chatY]=meshgrid(chatx, chaty);

smoothfactor = 20;
x = 1:3:size(chat,2);
y = 1:3:size(chat,1);
[X Y] = meshgrid(x,y);

zminmesh = gridfit(chatX, chatY,chatzmin, x, y, ...
        'regularizer', 'gradient', 'smooth', smoothfactor);
zmaxmesh = gridfit(chatX, chatY,chatzmax, x, y, ...
        'regularizer', 'gradient', 'smooth', smoothfactor);

enhanced=[];


zminmesh=max(1,min(zminmesh,size(chat,3)));
zmaxmesh=max(1,min(zmaxmesh,size(chat,3)));


% assign nuclear centers



% create a normalized plot

chatstruct.chatX=chatX;
chatstruct.chatY=chatY;
chatstruct.zminmesh=zminmesh;
chatstruct.zmaxmesh=zmaxmesh;
chatstruct.chatzmin=chatzmin;
chatstruct.chatzmax=chatzmax;
centroids = cat(1, props.Centroid);

% following code has problems
chatstruct.nuclgcl=centroids(zminmesh((1+floor((centroids(:,1)-1)/50))*25+floor((-1+centroids(:,2))/50)+1)>centroids(:,3),:);
chatstruct.nuclinl=centroids(zmaxmesh((1+floor((centroids(:,1)-1)/50))*25+floor((-1+centroids(:,2))/50)+1)<centroids(:,3),:);
chatstruct.nuclgcl(:,1)=chatstruct.nuclgcl(:,1)*5;
chatstruct.nuclinl(:,1)=chatstruct.nuclinl(:,1)*5;
chatstruct.nuclgcl(:,2)=chatstruct.nuclgcl(:,2)*5;
chatstruct.nuclinl(:,2)=chatstruct.nuclinl(:,2)*5;



projchatx=squeeze(mean(chat,2));
projchaty=squeeze(mean(chat,1));

visualize=0;
if (visualize)
                    figure(14), clf;
                   
                    surf(X, Y, chatstruct.zminmesh, 'EdgeColor', 'none');
                    hold on
                    surf(X, Y, chatstruct.zmaxmesh, 'EdgeColor', 'none');
                    shading interp, camlight left;
                    scatter3(chatstruct.chatX(:), chatstruct.chatY(:), chatstruct.chatzmin(:), 'bo', 'filled');
                    scatter3(chatstruct.chatX(:), chatstruct.chatY(:), chatstruct.chatzmax(:), 'ro', 'filled');
                    % axis image
                    set(gca, 'XLim', [1 1500]);
                    set(gca, 'YLim', [1 1500]);
                    set(gca, 'ZLim', [1 size(chat,3)]);
                    title('DAPI IPL Boundary Localization')
                    colorbar
end
                
browse=0;
x = 1:3:size(chat,2);
y = 1:3:size(chat,1);
[X Y] = meshgrid(x,y);

%calculate an index for every pixel
[xi,yi]=meshgrid(1:size(chat,2),1:size(chat,1));
meanmin=mean(zminmesh(:));
meanmax=mean(zmaxmesh(:));
VZminmesh=interp2(X,Y,zminmesh,xi,yi,'*linear',meanmin);
VZmaxmesh=interp2(X,Y,zmaxmesh,xi,yi,'*linear',meanmax);


newstack=zeros(size(chat),'single');
newstack((1:size(chat,1)*size(chat,2))'+floor(VZmaxmesh(:))*size(chat,1)*size(chat,2))=3;
newstack((1:size(chat,1)*size(chat,2))'+floor(VZminmesh(:))*size(chat,1)*size(chat,2))=4;   

chatstruct.VZminmesh=VZminmesh;
chatstruct.VZmaxmesh=VZmaxmesh;
chatstruct.chatstack=newstack;

if (browse)
BrowseComponents('ii',chat,newstack);
end
