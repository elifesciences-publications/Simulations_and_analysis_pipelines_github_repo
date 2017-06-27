function plotLDBoxes(hAx, startDark, endDark,TONORM)
%Plot light-dark boxes on passed axes handle hAx.

    % plot light/dark as dark boxes
    axes(hAx);
    if TONORM ==1 %if data is normalized to 0, stdev=100
        axlim=get(hAx,'ylim');
        low = axlim(1) + 0.05*(axlim(2)-axlim(1)); %-240;
        high = low+0.05*(axlim(2)-axlim(1)); %-200;
    else %if raw %P KaiC
        low=-7;
        high=3;
    end
    
    %plot light/dark as dark boxes
    yellowCol = [242 230 194]/255;
    greyCol=[212 222 224]/255;
      
    v=[startDark low; startDark high; endDark high; endDark low];
    f=[1 2 3 4];
    patch('Faces',f,'Vertices',v,...
        'FaceColor',greyCol,'EdgeColor','k');
    
    v2=[0 low; 0 high; startDark high; startDark low];
    patch2=patch('Faces',f,'Vertices',v2,...
        'FaceColor',yellowCol,'EdgeColor','k');
    
    v3=[endDark low; endDark high; ...
        max(get(hAx,'xlim')) high; max(get(hAx,'xlim')) low];
    patch3=patch('Faces',f,'Vertices',v3,...
        'FaceColor',yellowCol,'EdgeColor','k');
    
    xline=[0:1:max(get(hAx,'xlim'))];
    plot(hAx,xline,low.*ones(size(xline)),'k','linewidth',1);
    plot(hAx,xline,high.*ones(size(xline)),'k','linewidth',1);
    yline=[low:1:high];
    plot(hAx,xline(end)*ones(size(yline)),yline,'k','linewidth',1);
end

