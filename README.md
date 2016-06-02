# auau200GeVD0HadronCorrelation
```bash
mkdir myAnalysis
cd myAnalysis

# Replace address below with your own fork if you have one
git clone git@github.com:amcw7777/auau200GeVD0HadronCorrelation.git
cd auau200GeVD0HadronCorrelation

# Clone LBNL PicoHFLib
git clone git@github.com:rnc-lbl/auau200GeVRun14.git

# Now you need to get StPicoDstMaker
git clone git@github.com:rnc-lbl/star-picoDst.git


# Clone StRefMultCorr
git clone git@github.com:GuannanXie/Run14AuAu200GeV_StRefMultCorr.git

# Link all needed code under one StRoot directory:
mkdir StRoot
sh makeLinks.sh

# Compile
starver SL16d
cons
root4star -l -b -q -x runPicoD0AnaMaker.C\(\"$FILELIST\",\"OUTPUTFILENAME.root\"\)
