function read_CA_data(mice_name,session,area_info,output_name)
% 
% Reads data from .eeg file the specific area recordings and stores it in hd5 format for
% diffeent channels.
% 
%       input:
%       mice_name           string, name of the mice according to
%                           documentation
%       session             the recorded session according to
%                           documenatiation
%       
%       area_info           indices manually read from excel file 'epos'
%                           sheet
%
%       outputname          the name of h5 file to save all the channels
%                           relevant to area
%
%     file_address = strcat(mice_name,'/',session,'/',session);
    
    [xml, rxml] = LoadXml(strcat(file_address,'.xml'));
    
    %counter for the total number channels of a specific area
    counter=1;
    
    %setting ch_idx including area specific channel indices to empty array
    ch_idx=[];
    %reading the index for channels of specific area
    %first going through the indices stored in AnatGrps
    
    for i = area_info
        %looking to data for indices associated in AnatGrps
        for j = xml.AnatGrps(i+1).Channels
            ch_idx=cat(1,ch_idx, j+1);
        end        
    end
    disp(size(ch_idx))
    [data originIndex] = LoadBinary(strcat(file_address,'.eeg'),ch_idx);
    
    hdf5write(strcat(mice_name,'_',output_name,'_',session,'.h5'), '/data', data);