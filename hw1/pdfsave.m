function h = pdfsave(fileName, varargin)
%PDFSAVE A function to save the current plot as a PDF.
%   Detailed explanation goes here
    p = inputParser;
    p.addRequired('fileName', @ischar);
    p.addParamValue('orientation', 'landscape', @ischar);
    p.addParamValue('papersize', [11 8.5], @(x) isvector(x) && length(x) == 2);
    p.addParamValue('type', 'pdf', @ischar);
    p.parse(fileName, varargin{:});
    r = p.Results;
    h = gcf;
    set(h, 'PaperOrientation', r.orientation, ...
           'PaperSize', r.papersize, ...
           'PaperPosition', [0 0 r.papersize]);
    saveas(h, r.fileName, r.type);
end

