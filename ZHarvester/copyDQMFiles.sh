#! /bin/bash
cat ~/.globus/usercert.pem /etc/ssl/certs/CERN*pem >  /tmp/$USER/certs.pem
certU=YOUR_PROXY_FROM_PROMPT
cat filelist.log | while read LINE
do
wget -r -k $LINE --certificate=$certU --ca-certificate=/tmp/$USER/certs.pem -P /eos/cms/store/group/comm_luminosity/ZCounting/DQMFiles2017/ 
echo $LINE
done
