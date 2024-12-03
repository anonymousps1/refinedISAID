name = "Skinny64_lastHalf"
rMax = 2
S = [12,6,9,0,1,10,2,11,3,8,5,13,4,14,7,15,]
P = [0,5,10,15,4,9,14,3,8,13,2,7,12,1,6,11,]
M = matrix(GF(2),[
    [1,0,1,1,],
    [1,0,0,0,],
    [0,1,1,0,],
    [1,0,1,0,],
    ])
linAsSbox = True
keyAfterMC = False
startAfterSSB = False
noLastMC = True
dedicatedPRESENTlastLayer = False
