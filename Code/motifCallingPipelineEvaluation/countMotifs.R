library(data.table)

dd55=fread("wc -l `find interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM45/ | grep unsorted` | perl -ne 's/ *// && print' | perl -ne 's/ /\t/ && print' | grep unsorted")
dd5=fread("wc -l `find interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM5/ | grep unsorted` | perl -ne 's/ *// && print' | perl -ne 's/ /\t/ && print' | grep unsorted")
dd6=fread("wc -l `find interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM6/ | grep unsorted` | perl -ne 's/ *// && print' | perl -ne 's/ /\t/ && print' | grep unsorted")

dd55[,names:=sub("interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM45/","",V2)]
setkey(dd55,names)
dd5[,names:=sub("interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM5/","",V2)]
setkey(dd5,names)
dd6[,names:=sub("interimData/motifInstances/HOCOMOCOVARIABLECONF/HOCOMOCO1EM6/","",V2)]
setkey(dd6,names)
tot=merge(dd55,dd5)
tot=merge(tot,dd6)

fewer=tot[,V1.y/V1]
more=tot[,V1.x/V1.y]
