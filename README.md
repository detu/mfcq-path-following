Path-following for dual-degenerate cases

examples:
predictorCorrectorDualDegenerate(@(p)case1(p), 0,1,[0;0;0],[0;0.5;0;0.5;0;0;0],1.2)
predictorCorrectorDualDegenerate(@(p)case1(p), 0,1,[0;0;0],[0;0;1;0;0;0;0],1.2)

predictorCorrector(@(p)case1(p), 0,1,[0;0;0],[0;0.5;0;0.5;0;0;0],1.2)
predictorCorrector(@(p)case1(p), 0,1,[0;0;0],[0;0;1;0;0;0;0],1.2)
