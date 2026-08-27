#!/bin/bash
echo -n 'Enter key for slb-swt:  '
key=$(read token junk; echo $token)
echo -n 'Enter key for slb-swt1: '
key1=$(read token junk; echo $token)

for seq in 1 3 5
do
    echo "=== starting azure-slb-swt-${seq} ==="
    docker update --restart=no azure-slb-swt-${seq}
    docker kill azure-slb-swt-${seq}
    docker rm azure-slb-swt-old-${seq}
    docker rename azure-slb-swt-${seq} azure-slb-swt-old-${seq}
    docker run -d --name azure-slb-swt-${seq} \
           -e AZP_URL=https://dev.azure.com/slb-swt -e AZP_TOKEN="${key}" \
           -e AZP_AGENT_NAME=`uname -n`-${seq} \
           -e AZP_POOL=Paal-Kvamme-Test \
           -v /var/run/docker.sock:/var/run/docker.sock dockeragent:latest
done

for seq in 2 4 6
do
    echo "=== starting azure-slb1-swt-${seq} ==="
    docker update --restart=no azure-slb1-swt-${seq}
    docker kill azure-slb1-swt-${seq}
    docker rm azure-slb1-swt-old-${seq}
    docker rename azure-slb1-swt-${seq} azure-slb1-swt-old-${seq}
    docker run -d --name azure-slb1-swt-${seq} \
           -e AZP_URL=https://dev.azure.com/slb1-swt \
           -e AZP_TOKEN="${key1}" \
           -e AZP_AGENT_NAME=`uname -n`-${seq} \
           -e AZP_POOL=Paal-Kvamme-Test \
           -v /var/run/docker.sock:/var/run/docker.sock dockeragent:latest
done

# Note! restarting the container might not work.
# Delete and create a new one after a machine reboot.
# Should probably fix "start.sh".

# Note: Should update this to a newer version.
