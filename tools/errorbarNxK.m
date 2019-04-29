function [h, hb] = errorbarNxK(mx,ex,facnames,legnames,colmap,basevalue,barwidth)

% example
% mx = [ [2;2.5;2.8]  [.05;.25;.9] ];
% ex = .5*rand(size(mx));
% facnames = {'x1','x2'};
% legnames = {'leg1','leg2','leg3'};
% figure;errorbarKxN(mx,ex,facnames);
% another example
% mx = [ [2;2.5;]  [.05;.25;] ];
% ex = .5*rand(size(mx));
% facnames = {'x1','x2'};
% legnames = {'leg1','leg2'};

if nargin<4, legnames = []; end
if nargin<6, basevalue = 0; end
if nargin<7, barwidth = []; end

cmap = [.5 .5 .5];
[h,hb]=errorbarKxN(mx,ex,facnames,legnames,cmap,basevalue,[]);
for k=1:length(hb)
%     set(h1b(k),'barlayout','stacked');
    set(hb(k),'FaceColor',colmap(k,:));
    if ~isempty(barwidth)
        set(hb(k),'BarWidth',barwidth);
        
    end
end

end