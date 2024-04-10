function ROI = comm_ROI(varargin)

[Sx Sy Sz] = size(varargin{1});

ROI=zeros(Sx,Sy,Sz);

tmp=0;
kk=1;
for kk=1:nargin
    tmp=tmp+varargin{kk};
end

ratio=nargin;

ind=find(tmp==nargin);

ROI(ind)=1;
