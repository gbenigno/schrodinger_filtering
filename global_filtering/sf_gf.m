function [EEG_sf, globfilt] = sf_gf(ch_vec,EEG_aas,EEG_input_for_gf,EEG_raw,globfilt,trigs)

% author: Gabriel Benigno

% performs the global filtering step of Schrodinger filtering

%%% inputs:
% ch_vec: vector of channels to loop through
% EEG_aas: dataset following AAS
% EEG_input_for_gf: the dataset to be inputted for global filtering
% EEG_raw: raw dataset
% globfilt: structure array containing various data from global filtering step; input as blank
% trigs: vector of slice acquisition trigger timings

%%% outputs:
% EEG_sf: output following entire schrodinger filtering algorithm
% globfilt: structure array containing various data from global filtering step


EEG_sf = EEG_aas;
EEG_sf.data( : , (trigs(1) : trigs(end)) - 75 ) = NaN;

EEG_sf_data = EEG_sf.data;

slicesPerSubts=6;

globfilt_h_and_Nh = globfilt.h_and_Nh;
globfilt_h_SF = globfilt.h_and_Nh_SF(:,1);
globfilt_Nh_SF = globfilt.h_and_Nh_SF(:,2);
globfilt_mse = globfilt.mse;
globfilt_y_scsa = globfilt.y_scsa;

for ch = ch_vec
    
    disp( ['gf ch ', num2str(ch)] )
    
    EEG_sf_data_ch = EEG_sf_data(ch,:);

    for subts = 1:140

        a = 1+(subts-1)*slicesPerSubts; % epoch to start at

        if subts == 1

            % on first subts, generate structure array globfilt, which
            % contains y_scsa for different Nh and matrix of (h,Nh) pairs
            [globfilt_h_and_Nh{ch}, globfilt_y_scsa{ch}] = sf_scsa(EEG_raw, slicesPerSubts, EEG_input_for_gf, a, ch, trigs);

            % next, select Nh_SF
            [Nh_SF, mse] = select_NhSF(EEG_raw, slicesPerSubts, a, trigs, EEG_input_for_gf.data(ch,:), globfilt_y_scsa{ch});
            
            if Nh_SF < 60
                % find closest Nh to 90
                Nh_vec = globfilt_h_and_Nh{ch}(:,2);
                Nh_diff = abs(Nh_vec - 90);
                Nh_SF = find(Nh_diff==min(Nh_diff));
                if length(Nh_SF) > 1
                      Nh_SF = Nh_SF(end);
                end
            end

            globfilt_mse{ch} = mse;
            globfilt_h_SF(ch) = Nh_SF;
            
            h_SF = globfilt_h_and_Nh{ch}( Nh_SF, 1);
            globfilt_Nh_SF(ch) = h_SF;

            t = (trigs(1):trigs(7))-75;
            
            EEG_sf_data_ch(t) = globfilt_y_scsa{ch}(:,Nh_SF);
            

        else

            thr=0.025;
            [t, ~, ~, ~, ~, ~] = int_thresh(EEG_raw, thr, slicesPerSubts, a, trigs);
            dx=1;
            gm=0.5;
            [EEG_sf_data_ch(t), ~, ~, ~, ~]  = SCSA(EEG_input_for_gf.data(ch,t), h_SF, dx, gm);

        end
    end
    
    EEG_sf_data(ch,:) = EEG_sf_data_ch;
    

end

globfilt.h_and_Nh = globfilt_h_and_Nh;
globfilt.h_and_Nh_SF(:,1) = globfilt_h_SF;
globfilt.h_and_Nh_SF(:,2) = globfilt_Nh_SF;
globfilt.mse = globfilt_mse;
globfilt.y_scsa = globfilt_y_scsa;

EEG_sf.data = EEG_sf_data;

end


