function globfilt = globfilt_struct_func()

% author: Gabriel Benigno

%%% preallocate globfilt structure array, which, for each channel, will contain the following data from the global filtering step:
% - the various y_scsa vectors for different Nh-values
% - h and Nh vectors for each channel
% - hSF and NhSF values for each channel
% - mse vector (Delta(Nh) ie equation 7a from paper) for each channel

globfilt.y_scsa = cell(30,1);
globfilt.h_and_Nh=cell(30,1);
globfilt.h_and_Nh_SF=zeros(30,2);
globfilt.mse=cell(30,1);

end

