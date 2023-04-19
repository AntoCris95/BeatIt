
%% Use this function to retrieve info about data directories

% Created in July 2020
% Written by Antonio Criscuolo

function datadir = getSessdir(sess)

datadir.data = 'E:';
datadir.main = 'E:\DATA\Humans - Beatit';

if nargin <1
    sess= input('Sess number?');
end

if sess == 1
    datadir.sess= fullfile(datadir.main);
    datadir.output= fullfile(datadir.main, 'DATAproc');
% elseif sess == 2
%     datadir.sess= fullfile(datadir.main, 'xxx');
%     datadir.output= fullfile(datadir.main, 'xxx', 'DATAproc');
end

% Set up SubjNames

files = dir(fullfile(datadir.sess, '*.hdr'));
for ff = 1:length(files)
    SubjName = strsplit(files(ff).name, '.'); SubjName = SubjName{1};
    datadir.SubjNames{ff} = SubjName;
end

datadir.NSubj= length(datadir.SubjNames);
datadir.sessN = sess;
datadir.NSeq = 96;


