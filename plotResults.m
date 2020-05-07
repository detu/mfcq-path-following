function plotResults(state)
%PLOTRESULTS Summary of this function goes here
% 
% [OUTPUTARGS] = PLOTRESULTS(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: suwartad $	$Date: 2017/06/21 13:56:10 $	$Revision: 0.1 $
% Copyright: Process Control Group - NTNU Trondheim 2017


plot(xmeasureAll_1pct(state,:));   hold on; 
plot(xmeasureAll_pfWoHc(state,:)); hold on; 
plot(xmeasureAll_pf(state,:))

end
