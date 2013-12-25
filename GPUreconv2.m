function GPUreconv2(inputhdf5file,outputhdf5file,gpunum)

% choose the GPU with the most memory
% need at least 3294705152 to work

[status,result]=system('./gpumeminfo');
display(result);
if status==0
  A= strread(result, '%s', 'delimiter', sprintf('\n')); 
  B= strread(A{2}, '%s', 'delimiter', ' ');
  
  fm=[];
  for i=1:length(B)
   fm=[fm str2num(B{i})];
  end
  [val,ind]=max(fm);
  if val<3294705152 
    return
  end
else
  return
end

try
display(['using GPU:' gpunum]);
gpu=['gpu' gpunum];
% calculate optimal partition of data
rawdata=single(permute(hdf5read(inputhdf5file,'/GFP'),[3 2 1]));
load currentnetworkv5_net249.mat
onesidex = 8; onesidey = 8; onesidez = 2;
rawdata = (rawdata*15/max(rawdata(:)))-7.5;
% prepare reconstruction by blocks
xpos = 1; ypos = 1;
CHUNKSIDELENGTH = floor(sqrt(3.3e6/(size(rawdata,3)+2*onesidez))); disp(strcat('near optimal chunk size: ',num2str(CHUNKSIDELENGTH)))
xposend = CHUNKSIDELENGTH; yposend = CHUNKSIDELENGTH; step=CHUNKSIDELENGTH;
%%%%%%%%%%%%%%%%%%%%
temp = min(rawdata(:))*ones(size(rawdata,1)+2*onesidex, size(rawdata,2)+2*onesidey, size(rawdata,3)+2*onesidez);
temp(onesidex+1:end-onesidex, onesidey+1:end-onesidey, onesidez+1:end-onesidez) = rawdata;
rawdata = temp; clear temp;
[xsize ysize zsize] = size(rawdata);
allout = zeros(size(rawdata),'uint8');
split_size=[ CHUNKSIDELENGTH CHUNKSIDELENGTH zsize];
% initialize on gpu
m = cnpkg2_mapdim_layers_fwd(m,split_size,1);
cns('init',m,gpu,'mean')

outxpos=onesidex+1; outypos=onesidey+1; outxposend=xposend-onesidex; outyposend=yposend-onesidey;
input=zeros(split_size,'single');
% reconstruct
while xpos < xsize -2*onesidex
  while ypos < ysize -2*onesidey
    input(1:xposend-xpos+1,1:yposend-ypos+1,1:zsize)=single(rawdata(xpos:xposend,ypos:yposend,1:zsize));
    cns('set',{m.layer_map.input,'val',reshape(input,1,CHUNKSIDELENGTH,CHUNKSIDELENGTH,zsize,1)});
    output = cns('step',m.step_map.fwd(1)+1,m.step_map.fwd(end),{m.layer_map.output,'val'});
    output=squeeze(output);
    allout(outxpos:outxposend,outypos:outyposend,onesidez+1:end-onesidez) = uint8(255*output(1:outxposend-outxpos+1,1:outyposend-outypos+1,1:zsize-2*onesidez));
    ypos = ypos + step-2*onesidey;
    yposend = min(yposend+step-2*onesidey, ysize);
    outypos = outyposend+1;
    outyposend = yposend-onesidey;
  end
  xpos = xpos + step-2*onesidex; xposend = min(xposend+step-2*onesidex, xsize);
  outxpos = outxposend+1; outxposend = xposend-onesidex; %outxposend = xposend-2*onesidex;
  ypos = 1; yposend = step; outypos=onesidey+1; outyposend=step-onesidey; %outyposend=step-2*onesidey;
end
cns done;
allout = allout(onesidex+1:end-onesidex,onesidey+1:end-onesidey,onesidez+1:end-onesidez);


projectionXY=max(allout,[],3)';
projectionYZ=squeeze(max(allout,[],1))';
projectionXZ=squeeze(max(allout,[],2))';
GFP=permute(allout,[3 2 1]);

save(outputhdf5file,'GFP','projectionXY','projectionXZ','projectionYZ','-v7.3');

catch
   cns done;
   display('error encountered')
   lasterr
end
