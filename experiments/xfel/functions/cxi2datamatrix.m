function [ data ] = cxi2datamatrix( filename )
%converts cxi file from Condor to a p x n matrix
% source: http://xfel.icm.uu.se/condor/static/docs/cxi.html#file-structure

fid = H5F.open(filename);
dset_id = H5D.open(fid,'/entry_1/data_1/data');

intensity_patterns = H5D.read(dset_id);

H5D.close(dset_id);
H5F.close(fid);

sz = size(intensity_patterns);

data = reshape( intensity_patterns, [], sz(end));


end

