function [ precoder ] = Precoder_selection_MC(codebook,H,Da_Str,noise_power)


codebook_size = 8;
%% Maximum capacity criterion

MC = zeros(1,codebook_size);
for ii = 1:codebook_size
    Fi = codebook(:,:,ii);
    MC(ii) = log2(det( eye(Da_Str)+ ( (H*Fi)'*(H*Fi) )/noise_power  )  );
end
[~,F_index] = max(MC);

precoder = codebook(:,:,F_index);


end
