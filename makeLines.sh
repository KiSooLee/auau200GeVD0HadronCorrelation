ln -s `pwd`/StPicoD0AnaMaker StRoot
ln -s `pwd`/auau200GeVRun14/StRoot/StPicoD0EventMaker StRoot
ln -s `pwd`/auau200GeVRun14/StRoot/StPicoPrescales StRoot
ln -s `pwd`/star-picoDst/StPicoDstMaker StRoot
ln -s `pwd`/Run14AuAu200GeV_StRefMultCorr/VPDMB5/StRefMultCorr StRoot
ln -s `pwd`/auau200GeVRun14/StRoot/StPicoKFVertexFitter StRoot
cp -r -p auau200GeVRun14/run14AuAu200GeVPrescales/ .
cp `pwd`/StPicoD0AnaMaker/efficiency.txt .
cp `pwd`/StPicoD0AnaMaker/loadKFLibraries.C .
cp `pwd`/StPicoD0AnaMaker/runPicoD0AnaMaker.C .
cp `pwd`/StPicoD0AnaMaker/phiFix.root .
