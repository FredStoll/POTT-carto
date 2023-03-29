function [idx,idxs] = findenough(data,thr,nbRequired,op)
% Find the position in a matrix (data) where a given number of values (nbRequired) exceed a given threshold (thr)
%- Author: Fred M. Stoll, Icahn School of Medicine at Mount Sinai, NY

eval(['oo=[-1 find(data',op,'thr)];']);
%eval(['oo=[-1 find(data>thr)];']);
xx = [find(diff(oo)~=1) length(oo)];
yy=[xx(1) diff(xx)];
zz=find(yy>=nbRequired);
%zz=find(yy>=nbRequired,1,'first');
idx=oo(xx(zz-1)+1) ;

idxs = [];parts = [];
for i = 1 : length(idx)
    idxs = [idxs, idx(i) : idx(i)+yy(zz(i))-1];
    parts = [parts, i*ones(1,length(idx(i) : idx(i)+yy(zz(i))-1))];
end
