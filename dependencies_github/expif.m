function expif(TOEXPORT, varargin)
%If TOEXPORT==1, export_fig(varargin)
% dependency: export_fig
if TOEXPORT == 1
    export_fig(varargin);
    feval(@export_fig, varargin);
end

end

