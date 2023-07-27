#!/bin/sh
cmsenv
cmsRun Scouting_cfg.py GlobalTagData=123X_dataRun3_HLT_v14 output=dataNtupleTrial.root inputFile=file:/vols/cms/pb4918/StoreNTuple/Scouting/2022FSample.root > anlzr.log 2>&1 
