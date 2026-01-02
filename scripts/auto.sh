#!/bin/bash

#sh autorunImp.sh processTraj ImpX20
#sleep 10
#sh autorunImp.sh rate ImpX20
#sleep 10
sh autorunImp.sh conductivity ImpX20
sleep 10

#sh autorunAll.sh rate LiPF6inACN
#sleep 10
sh autorunAll.sh conductivity LiPF6inACN
sleep 10
#sh autorunAll.sh rate LiPF6inH2O
#sleep 10
sh autorunAll.sh conductivity LiPF6inH2O
sleep 10

sh autorunAll.sh conductivity KPF6inACN
sleep 10
