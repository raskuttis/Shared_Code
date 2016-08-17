function [FV, LB] = feature_out(featname, nclust, nex, fisher, pcaval);

datadir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/voxResources';
run('/Users/sudhirraskutti/Desktop/Source/vlfeat/toolbox/vl_setup');
outdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Features';

%load all songs into a single struct
[DAT, LB, FNS] = loadAll(datadir);

% extract a single feature
feat = cell(1,1000);
for i = 1:length(DAT)
    switch featname
        case 'mfc'
            feat{i} = DAT{i}.mfc;
        case 'eng'
            feat{i} = DAT{i}.eng;
        case 'chroma'
            feat{i} = DAT{i}.chroma;
        case 't'
            feat{i} = DAT{i}.t;
        case 'keystrength'
            feat{i} = DAT{i}.keystrength;
        case 'brightness'
            feat{i} = DAT{i}.brightness;
        case 'zerocross'
            feat{i} = DAT{i}.zerocross;
        case 'roughness'
            feat{i} = DAT{i}.roughness;
        case 'inharmonic'
            feat{i} = DAT{i}.inharmonic;
        case 'hcdf'
            feat{i} = DAT{i}.hcdf;
        case 'key'
            feat{i} = DAT{i}.key;
        case 'tempo'
            feat{i} = DAT{i}.tempo;
        otherwise
    end
end

%create the structure used as input into the fisher vector extraction
if fisher == 1
    GENDATA.data = feat;
    GENDATA.class = LB;
    GENDATA.classnames = {'Blues', 'Classical', 'Country', 'Disco', 'Hiphop',...
        'Jazz', 'Metal', 'Pop', 'Reggae', 'Rock'};

    %run fisher vector
    FV = find_fv(GENDATA, nclust, nex, pcaval);
    FVdim = size(FV);
    pcadim = FVdim(1);

    % save the fisher vector to the output directory
    if pcaval >= 99.9
        outfv = strcat('_', featname, '_GMM', int2str(nclust), '_Nex', int2str(nex), '_All', '.mat')
    else
        outfv = strcat('_', featname, '_GMM', int2str(nclust), '_Nex', int2str(nex), '_PCA', num2str(pcaval), '.mat')
    end
else
    FV = feat;
    outfv = strcat('_', featname, '_All', '.mat');
end

save(strcat(outdir,'/FV',outfv),'FV');
save(strcat(outdir,'/LB',outfv),'LB');


