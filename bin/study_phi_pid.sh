#!/opt/local/bin/zsh



# ./pico.app config/phi_histo.xml --p.nSigmaDeltaY=3 --p.nSigmaDeltaZ=3 --p.DeltaTOF=1.0 --name="-TRA331p0"
# ./pico.app config/phi_histo.xml --p.nSigmaDeltaY=2 --p.nSigmaDeltaZ=2 --p.DeltaTOF=0.5 --name="-TRA220p5"
# ./pico.app config/phi_histo.xml --p.nSigmaDeltaY=2 --p.nSigmaDeltaZ=2 --p.DeltaTOF=0.25 --name="-TRA220p25"

# ./pico.app config/phi_mvahisto.xml --p.pid=1.2 --name="-DNN1p2"
# ./pico.app config/phi_mvahisto.xml --p.pid=1.3 --name="-DNN1p3"
# ./pico.app config/phi_mvahisto.xml --p.pid=1.35 --name="-DNN1p35"
# ./pico.app config/phi_mvahisto.xml --p.pid=1.36 --name="-DNN1p36"
# ./pico.app config/phi_mvahisto.xml --p.pid=1.37 --name="-DNN1p37"
# ./pico.app config/phi_mvahisto.xml --p.pid=1.38 --name="-DNN1p38"
# ./pico.app config/phi_mvahisto.xml --p.pid=1.39 --name="-DNN1p39"
# ./pico.app config/phi_mvahisto.xml --p.pid=1.4 --name="-DNN1p4"


# analysis

./pico.app config/phi_fitter.xml --name="-TRA331p0" --title="Traditional-J/Psi"
./pico.app config/phi_fitter.xml --name="-TRA220p5" --title="Traditional-Tighter"
./pico.app config/phi_fitter.xml --name="-TRA220p25" --title="Traditional-Tightest"

./pico.app config/phi_fitter.xml --name="-DNN1p2" --title="DNN>1.2"
./pico.app config/phi_fitter.xml --name="-DNN1p3" --title="DNN>1.3"
./pico.app config/phi_fitter.xml --name="-DNN1p35" --title="DNN>1.35"
./pico.app config/phi_fitter.xml --name="-DNN1p36" --title="DNN>1.36"
./pico.app config/phi_fitter.xml --name="-DNN1p37" --title="DNN>1.37"
./pico.app config/phi_fitter.xml --name="-DNN1p38" --title="DNN>1.38"
./pico.app config/phi_fitter.xml --name="-DNN1p39" --title="DNN>1.39"
./pico.app config/phi_fitter.xml --name="-DNN1p4" --title="DNN>1.4"