genericRecoverySim(simPars[1, ], cuPar=CUpars, srDat=recDatTrim,
        variableCU=FALSE, ricPars=mcmcOut, cuCustomCorrMat = corMatrix,
         nTrials=nProj, makeSubDirs=FALSE, random=FALSE, outDir="C:/github/SalmonLRP_CJFAS/WCVIChinookStudy")

genericRecoverySim(x, cuPar=CUpars, srDat=recDatTrim,
                   variableCU=FALSE, ricPars=mcmcOut,
                   cuCustomCorrMat = corMatrix,
                   nTrials=nProj, makeSubDirs=FALSE,
                   random=FALSE, outDir=outDir)



simPar <- simPars[1,]
cuPar <- CUpars
srDat <- NULL#recDatTrim#NULL#recDatTrim%>%mutate(rec2=NA, rec3=NA, rec4=NA, rec5=NA, rec6=NA)
variableCU <- FALSE
ricPars <- mcmcOut
cuCustomCorrMat <- corMatrix
nTrials <- nProj
makeSubDirs <- FALSE
random <- FALSE
outDir <- outDir

catchDat <- NULL
larkPars <- NULL
erCorrMat <- NULL
uniqueProd <- TRUE
uniqueSurv <- FALSE


recDatTrim <- recDatTrim%>%mutate(rec2=NA, rec3=NA, rec4=NA, rec5=NA, rec6=NA)
