function [t0, x0, pQ, pC] = measureInitialValueParam(tmeasure, xmeasure, mpciter, pQ, pC)
%MEASUREINITIALVALUEPARAM Summary of this function goes here
% 
% [OUTPUTARGS] = MEASUREINITIALVALUEPARAM(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2018/05/04 01:18:39 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2018

%         pC = 2;
%         pQ = 0.5;
%         %x0 = [0.636232215874424;	0.431505717792292]; %hardcode!
%         t0 = tmeasure;
%         x0 = xmeasure;

switch mpciter
    
    case 1
        pC = 2;
        pQ = 0.5;
        %x0 = [0.636232215874424;	0.431505717792292]; %hardcode!
        t0 = tmeasure;
        x0 = xmeasure;
    case 25
        pC = 4;
        pQ = 0.25;
        t0 = tmeasure;
        x0 = xmeasure;
    case 80
        pC = 1;
        pQ = 1;
        t0 = tmeasure;
        x0 = xmeasure;
        
    otherwise
        t0 = tmeasure;
        x0 = xmeasure;
end

end
