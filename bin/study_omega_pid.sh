#!/opt/local/bin/zsh


# run study_phi_pid.sh first to gen files

# analysis

./pico.app config/phi_fitter.xml --p.reso="omega" --name="-TRA331p0" --title="Traditional-J/Psi"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-TRA220p5" --title="Traditional-Tighter"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-TRA220p25" --title="Traditional-Tightest"

./pico.app config/phi_fitter.xml --p.reso="omega" --name="-DNN1p2" --title="DNN>1.2"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-DNN1p3" --title="DNN>1.3"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-DNN1p35" --title="DNN>1.35"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-DNN1p36" --title="DNN>1.36"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-DNN1p37" --title="DNN>1.37"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-DNN1p38" --title="DNN>1.38"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-DNN1p39" --title="DNN>1.39"
./pico.app config/phi_fitter.xml --p.reso="omega" --name="-DNN1p4" --title="DNN>1.4"