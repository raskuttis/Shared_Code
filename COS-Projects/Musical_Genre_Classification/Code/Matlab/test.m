function [FV, LB] = feature_out(featname, nclust, nex, pcaval)

datadir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/voxResources';
run('/Users/sudhirraskutti/Desktop/Source/vlfeat/toolbox/vl_setup');
outdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/HW1/Resources/Features';

nclust = 3;
nex = 3;

%load all songs into a single struct
[DAT, LB, FNS] = loadAll(datadir);

DAT

% extract a single feature
feat = cell(1,1000);
for i = 1:length(DAT)
    feat{i} = DAT{i}.mfc;
end

%create the structure used as input into the fisher vector extraction
GENDATA.data = feat;
GENDATA.class = LB;
GENDATA.classnames = {'Blues', 'Classical', 'Country', 'Disco', 'Hiphop',...
    'Jazz', 'Metal', 'Pop', 'Reggae', 'Rock'};

%run fisher vector
FV = find_fv(GENDATA, nclust, nex);
FVdim = size(FV);
pcadim = FVdim(1);

% save the fisher vector to the output directory
outfv = strcat('_MFCC_GMM', int2str(nclust), '_Nex', int2str(nclust), '_PCA', int2str(pcadim), '.mat')
% save(strcat(outdir,'/FV',outfv),'FV');

