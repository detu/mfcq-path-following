% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%addpath('C:\Data\work\casadi-matlabR2014b-v3.1.0-rc1');
%addpath('C:\Data\work\casadi-matlabR2014b-v3.2.0');
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64');
%addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\examples\src\matlab');

% disable warning message
warning off;
