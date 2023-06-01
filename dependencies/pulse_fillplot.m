function pulse_fillplot(pt,pm,pstd,clr)

    curve_1 = (pm + pstd);
    curve_2 = (pm - pstd);
    x = pt;    
    
    plot(pt,pm,'color',clr)
    hold on
    x2 = [x, fliplr(x)];
    inBetween = [curve_1, fliplr(curve_2)];
    fill(x2,inBetween,'-','FaceColor',clr,'FaceAlpha',0.1,'EdgeColor',clr,'EdgeAlpha',0.1);

end