# ALLEGRO_PandoraPFA
Instructions and scripts to run the PandoraPFA reconstruction within ALLEGRO full simulation

## Setup
### Setup the environment
(I am working with release 2024-10-03)
```
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-10-03
```

### Checkout the detector related packages:
```
git clone https://github.com/Archil-AD/k4RecTracker.git
git clone https://github.com/Archil-AD/k4RecCalorimeter.git
git clone https://github.com/Archil-AD/k4geo.git
```

### Checkout the PandoraPFA related packages:
```
git clone https://github.com/Archil-AD/LCContent.git
git clone https://github.com/Archil-AD/DDMarlinPandora.git
```

### Compile packages:
```
cd k4geo/
# checkout the pandora branch
git checkout pandora
mkdir build install
cd build/
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ..
k4_local_repo
export K4GEO=$PWD/
cd ../
```

```
cd k4RecCalorimeter/
mkdir build install
cd build/
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ..
k4_local_repo
cd ../
```

```
cd k4RecTracker/
# checkout the pandora branch
git checkout pandora
mkdir build install
cd build/
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ..
k4_local_repo
cd ../
```

```
cd LCContent/
# checkout the ALLEGRO branch
git checkout ALLEGRO
mkdir build install
cd build/
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ..
k4_local_repo
cd ../
```

```
cd DDMarlinPandora/
# checkout the ALLEGRO branch
git checkout ALLEGRO
mkdir build install
cd build/
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ..
k4_local_repo
cd ../
```

## Simulation

Run the simulation of single photon/Kaon0L/electron/pion
```
ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "60*deg" --gun.thetaMax "60*deg" --gun.particle gamma --numberOfEvents 100 --outputFile ALLEGRO_sim_photon.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v04/ALLEGRO_o1_v04.xml
ddsim --enableGun --gun.distribution uniform --gun.energy "50*GeV" --gun.thetaMin "60*deg" --gun.thetaMax "60*deg" --gun.particle kaon0L --numberOfEvents 100 --outputFile ALLEGRO_sim_kaon0L.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v04/ALLEGRO_o1_v04.xml
ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.thetaMin "60*deg" --gun.thetaMax "60*deg" --gun.particle e- --numberOfEvents 100 --outputFile ALLEGRO_sim_electron.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v04/ALLEGRO_o1_v04.xml
ddsim --enableGun --gun.distribution uniform --gun.energy "50*GeV" --gun.thetaMin "60*deg" --gun.thetaMax "60*deg" --gun.particle pi- --numberOfEvents 100 --outputFile ALLEGRO_sim_pi.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v04/ALLEGRO_o1_v04.xml
```

## Reconstruction

Run the reconstruction:
```
k4run ALLEGROReconstruction.py --inputFiles ALLEGRO_sim_photon.root  --outputFile ALLEGRO_reco_photon.root
k4run ALLEGROReconstruction.py --inputFiles ALLEGRO_sim_kaon0L.root  --outputFile ALLEGRO_reco_kaon0L.root
k4run ALLEGROReconstruction.py --inputFiles ALLEGRO_sim_electron.root  --outputFile ALLEGRO_reco_electron.root
k4run ALLEGROReconstruction.py --inputFiles ALLEGRO_sim_pi.root  --outputFile ALLEGRO_reco_pi.root
```

The content of the output file can be checked using podio-dump:
```
podio-dump ALLEGRO_reco_kaon0L.root
```

An example of drawing reonstructed energy using Particle Flow Objects (PFO):
```
root -l	ALLEGRO_reco_kaon0L.root
events->Draw("Sum$(PandoraPFANewPFOs.energy)","","")
# compare it with the energy reconstruction using all cells of HCal and ECal barrels
events->Draw("(Sum$(HCalBarrelReadoutPositioned.energy)+Sum$(ECalBarrelModuleThetaMergedPositioned.energy))","")
```
