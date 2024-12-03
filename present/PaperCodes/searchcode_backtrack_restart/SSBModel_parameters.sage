name = "Skinny64_SSB"
S = [12,6,9,0,1,10,2,11,3,8,5,13,4,14,7,15,]
M = matrix(GF(2),[
    [1,0,1,1,],
    [1,0,0,0,],
    [0,1,1,0,],
    [1,0,1,0,],
    ])
sboxSize = 4
blockSize = 16
linAsSbox = True
keyAfterMC = False
dedicatedPRESENTlastLayer = False
