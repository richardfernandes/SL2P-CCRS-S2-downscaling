% Downscale SAFE S2 images in a given folder using ATPRK
% Parameters
inputDirectoryName = "F:\novascotia\SAFEimagesJulyAug\";
workingDirectoryName = "F:\novascotia\code\Code-for-S2-fusion-master\";

% Get list of input SAFE zip archives
inputDirectoryList = dir(strcat(inputDirectoryName,'*.SAFE.zip'));

% change over to working directory
cd(workingDirectoryName) 

% Loop thorough 
for k=1:length(inputDirectoryList)
    
    %unzip safe file into working directory
    inputFileName = inputDirectoryList(k).name()
    filenames = unzip(strcat(inputDirectoryName,inputFileName));

    %read in required files (we assume there is only one granule?)
    B02=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B02_10m')))));
    B03=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B03_10m')))));
    B04=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B04_10m')))));
    B08=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B08_10m')))));
    B05=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B05_20m')))));
    B06=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B06_20m')))));
    B07=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B07_20m')))));
    B8A=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B8A_20m')))));
    B11=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B11_20m')))));
    B12=imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'B12_20m')))));
    MSK_CLDPRB_20m = imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'MSK_CLDPRB_20m')))));
    SCL_20m = imread(strcat(workingDirectoryName,filenames(find(contains(filenames,'SCL_20m')))));
    %save the images as cubes
    S2_10m = cat(3,B02,B03,B04,B08);
    S2_20m = cat(3,B05,B06,B07,B8A,B11,B12);

    %clean up workspace
    clear('B*')
    system(strcat('rm -r ',workingDirectoryName,inputFileName(1:(end-4))))

    %downscale 
    S2_20m_downscaled =BestBandATPRK(S2_10m,S2_20m,MSK_CLDPRB_20m,SCL_20m);

    %save downscaled result band by band
    outputFilePath = strcat(workingDirectoryName,filenames(find(contains(filenames,'B02_10m'))));
    imwrite(squeeze(S2_20m_downscaled(:,:,1)),strrep(outputFilePath,'B02_10m','B05_10m'));
    imwrite(squeeze(S2_20m_downscaled(:,:,2)),strrep(outputFilePath,'B02_10m','B06_10m'));
    imwrite(squeeze(S2_20m_downscaled(:,:,3)),strrep(outputFilePath,'B02_10m','B07_10m'));
    imwrite(squeeze(S2_20m_downscaled(:,:,4)),strrep(outputFilePath,'B02_10m','B8A_10m'));
    imwrite(squeeze(S2_20m_downscaled(:,:,5)),strrep(outputFilePath,'B02_10m','B11_10m'));
    imwrite(squeeze(S2_20m_downscaled(:,:,6)),strrep(outputFilePath,'B02_10m','B12_10m'));


    end