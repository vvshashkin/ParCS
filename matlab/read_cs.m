function D=read_cs(file, nx)

fid = fopen(file)
D = fread(fid, 'float');
D = squeeze(reshape(D, nx, nx, 6, []));

end
