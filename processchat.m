function processchat(inputhdf5file,outputhdf5file)
warning off all
a.chat=permute(hdf5read(inputhdf5file,'/TRITC'),[3 2 1]);
[chat]=processchatNN(a.chat);
clear a
chatX=chat.chatX';
chatY=chat.chatY';
zminmesh=chat.zminmesh';
zmaxmesh=chat.zmaxmesh';
chatzmin=chat.chatzmin';
chatzmax=chat.chatzmax';
nuclgcl=chat.nuclgcl';
nuclinl=chat.nuclinl';
VZminmesh=chat.VZminmesh';
VZmaxmesh=chat.VZmaxmesh';
%a.gfp=permute(hdf5read(inputhdf5file,'/GFP'),[3 2 1]);
%for i=1:size(a.gfp,3)
%    a.gfp(:,:,i)=medfilt2(a.gfp(:,:,i),[5 5]);
%end
chatstack=uint16(permute(chat.chatstack,[3 2 1]));
clear chat
clear a
save(outputhdf5file,'-v7.3')
