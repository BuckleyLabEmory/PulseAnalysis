function graph_beautify(axis_handle,text_size,line_width,box_width)
    box(axis_handle,'on')
    set(axis_handle,'fontsize',text_size) % 10
    set(axis_handle,'linewidth',box_width) % 3
    set(findobj(axis_handle,'type','line'),'linew',line_width)
    set(gcf,'color','w');
end
