function activity=NNpass(W,B,sRawImg,num_of_maps_in_layer);
	activity{1}(:,:,:,1)=sRawImg;
    	nonlinearity=inline('1./(1+exp(-x))'); % sigmoid function
	for l=2:size(W,1)
		for feature_map=1:num_of_maps_in_layer{l}
			activity{l}(:,:,:,feature_map)=nonlinearity(convn(activity{l-1}, W{l, feature_map}, 'valid') + B{l,feature_map});
		end
	end
end

