function GPUrecon_bands(inputhdf5file,outputhdf5file)

if exist('/tmp/gpu0lock')
    if exist('/tmp/gpu1lock')
        return
    else
        f=fopen('/tmp/gpu1lock','w');
        if ~f
            return
        end;
        fclose(f);
        gpunum=1;
    end
else
     f=fopen('/tmp/gpu0lock','w');
         if ~f
        return
    end;
         fclose(f);
     gpunum=0;
end

try

display(['using GPU:' num2str(gpunum)]);
gpu=['gpu' num2str(gpunum)];


rawdata=single(permute(hdf5read(inputhdf5file,'/TRITC'),[2 3 1]));
rawsize=size(rawdata);
load bandCNN.mat
onesidex = 22; onesidey = 22; onesidez = 3;

% prepare reconstruction by blocks
xpos = 1; ypos = 1;
CHUNKSIDELENGTH = floor(sqrt(3.3e6/(size(rawdata,3)+2*onesidez))); disp(strcat('near optimal chunk size: ',num2str(CHUNKSIDELENGTH)))
xposend = CHUNKSIDELENGTH; yposend = CHUNKSIDELENGTH; step=CHUNKSIDELENGTH;

%normalization
rawdata = rawdata - mean(rawdata(:)); rawdata = rawdata/sqrt(cov(rawdata(:)));
%%%%%%%%%%%%%%%%%%%%
temp = min(rawdata(:))*ones(size(rawdata,1)+2*onesidex, size(rawdata,2)+2*onesidey, size(rawdata,3)+2*onesidez);
temp(onesidex+1:end-onesidex, onesidey+1:end-onesidey, onesidez+1:end-onesidez) = rawdata;

rawdata = temp; clear temp;
[xsize ysize zsize] = size(rawdata);
allout = zeros(rawsize,'uint8');

split_size=[ CHUNKSIDELENGTH CHUNKSIDELENGTH zsize];


% initialize on gpu
m = cnpkg2_mapdim_layers_fwd(m,split_size,1);
cns('init',m,gpu,'mean')

outxpos=1; outypos=1; outxposend=xposend-2*onesidex; outyposend=yposend-2*onesidey;
input=zeros(split_size,'single');
% reconstruct
while xpos < xsize -2*onesidex
  while ypos < ysize -2*onesidey
    tic;
    disp([xpos ypos])
    % pad input
    input=single(rawdata(xpos:xposend,ypos:yposend,1:zsize));
    %input=input*0; input(1:xposend-xpos+1,1:yposend-ypos+1,1:zsize)=single(rawdata(xpos:xposend,ypos:yposend,1:zsize));
    cns('set',{m.layer_map.input,'val',reshape(input,1,CHUNKSIDELENGTH,CHUNKSIDELENGTH,zsize,1)});
    output = cns('step',m.step_map.fwd(1)+1,m.step_map.fwd(end),{m.layer_map.output,'val'});
    %output = output(:,onesidex+1:end-onesidex,onesidey+1:end-onesidey,onesidez+1:end-onesidez);
    output=squeeze(output);
    allout(outxpos:outxposend,outypos:outyposend,:) = uint8(255*output(1:outxposend-outxpos+1,1:outyposend-outypos+1,1:zsize-2*onesidez));
    ypos = ypos + step-2*onesidey;
    yposend = min(yposend+step-2*onesidey, ysize);
    outypos = outyposend+1;
    outyposend = yposend-onesidey*2;
toc;
  end
  xpos = xpos + step-2*onesidex; xposend = min(xposend+step-2*onesidex, xsize);
  outxpos = outxposend+1; outxposend = xposend-onesidex*2; %outxposend = xposend-2*onesidex;
  ypos = 1; yposend = step; outypos=onesidey+1; outyposend=step-onesidey*2; %outyposend=step-2*onesidey;
end
cns done;
delete(['/tmp/' gpu 'lock']);

TRITC=permute(allout,[2 3 1]);
save(outputhdf5file,'TRITC','-v7.3');

catch
   cns done;
   delete(['/tmp/' gpu 'lock']);
   lasterr
end
