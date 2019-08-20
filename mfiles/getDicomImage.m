function [imgData, hdr, acqTime, flipAngles] = getDicomImage(dcmDir)
% Simple DICOM file reader, designed for Sarcoma DCE data
% Usage: [imgData, hdr, acqTime] = getDicomImage(dcmDir)
% dcmDir - Path to folder containing DICOM files
% imgData - Array containing loaded DICOM images
% hdr - hdr for the most recently read dicom files
% acqTime - Time at which frame was acquired (useful for DCE acquisitions)

    dcmFiles = dir([dcmDir '/*.dcm']);

    tmpImg = dicomread(fullfile(dcmDir,dcmFiles(1).name));
    [sX, sY] = size(tmpImg);
    sZ = length(dcmFiles);

    imgData = zeros(sX,sY,sZ);

    acqTime = zeros(sZ,1);
    flipAngles = acqTime;
    for i=1:length(dcmFiles)
       hdr = dicominfo(fullfile(dcmDir,dcmFiles(i).name));
       imgData(:,:,hdr.InstanceNumber) = dicomread(hdr);
       if isfield(hdr,'TriggerTime')
           acqTime(hdr.InstanceNumber) = hdr.TriggerTime;
       else
           acqTime(hdr.InstanceNumber) = str2num(hdr.AcquisitionTime);
       end
       if isfield(hdr,'FlipAngle')
           flipAngles(hdr.InstanceNumber)=hdr.FlipAngle;
       end
    end

end
