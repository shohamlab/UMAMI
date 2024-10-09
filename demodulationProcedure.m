function [ModulImageabs,FilterFFTImage] = demodulationProcedure(imagi,ACFerq,GPU)
% imagi is the 2D matrix that will be filtered in row 
% ACFerq is the [a b] frequency band that will be filtered
% GPU, 0 or 1, is whether the computer has a GPU available

    fftImage = fftshift(fft(fftshift(imagi')))'; 
    if GPU == 1
        FilterFFTImage = gpuArray(zeros(size(fftImage)));
    else
        FilterFFTImage = (zeros(size(fftImage)));
    end
    c = floor((ACFerq(2)+ACFerq(1))./2);   %center
    m = floor((ACFerq(2)-ACFerq(1))./2);
    
    FilterFFTImage(:,c-m:c+m) = fftImage(:,ACFerq(1):ACFerq(2));
    firstModulImage = fftshift(ifft(fftshift(FilterFFTImage')))';
    if GPU == 1
        ModulImageabs = gather(abs(firstModulImage));
        FilterFFTImage = gather(FilterFFTImage);
    else
        ModulImageabs = (abs(firstModulImage));
    end
end 