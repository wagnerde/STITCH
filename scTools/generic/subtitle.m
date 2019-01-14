function p = subtitle(string, Pos, dxy)
% p = subtitle(string, Pos, dxy)
% 
% Adds a text string to plot, locating string by default at inner-top-left
% Uses fontsize of current axes (gca)
%
% output: p returns pointer to the text string
% input:    string - the text
%           Pos (optional) - 'TopLeft', 'TopRight',
%           'BottomLeft','BottomRight'
%
%           dxy (optional) - 2x1 vector for manual adjustments to subtitle
%           position. Values from 0-1 give fractional of axis

if(~exist('Pos','var'))
    Pos = 'TopLeft';
end
if(~exist('dxy'))
    dxy = [0 0];
end

xlim = get(gca,'xlim');
ylim = get(gca,'ylim');

dxy(1) = diff(xlim)*dxy(1);
dxy(2) = diff(ylim)*dxy(2);

if(~isempty(strfind(Pos,'Top')))
    y = ylim(2)-dxy(2);
    va = 'Top';
else
    y = ylim(1)+dxy(2);
    va = 'Bottom';
end

if(~isempty(strfind(Pos,'Left')))
    x = xlim(1)+dxy(1);
    ha = 'Left';
else
    x = xlim(2)-dxy(1);
    ha = 'Right';
end



p = text(x, y, string, 'HorizontalAlignment',ha,'VerticalAlignment',va,...
    'FontSize', get(gca,'FontSize'));

