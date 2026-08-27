# Copyright 2017-2021, Schlumberger
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Based on https://docs.microsoft.com/en-us/azure/devops/pipelines/agents/docker?view=azure-devops as of 2024-07-17.

FROM ubuntu:22.04

RUN apt update
RUN apt upgrade -y
RUN apt install -y curl git jq libicu70

WORKDIR /azp/

COPY ./start.sh ./
RUN chmod +x ./start.sh

# Create agent user and set up home directory
RUN groupadd --system docker && \
    useradd -m -d /home/agent -G docker agent

#############################################################################
## Local customization, originally for Ubuntu 18.04, may need tweaks. #######
#############################################################################

## Add Docker's official GPG key:
#RUN apt-get install ca-certificates curl &&
#    install -m 0755 -d /etc/apt/keyrings &&
#    curl -fsSL https://download.docker.com/linux/ubuntu/gpg \
#         -o /etc/apt/keyrings/docker.asc &&
#    chmod a+r /etc/apt/keyrings/docker.asc
#
## Add the repository to Apt sources:
#RUN echo \
#  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
#  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
#  tee /etc/apt/sources.list.d/docker.list > /dev/null && \
# apt-get update
#
## Probably a massive overkill. I just need the client
#RUN apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

COPY --from=docker:dind /usr/local/bin/docker /usr/local/bin/

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        iputils-ping \
        libcurl4 \
        libunwind8 \
        netcat \
        libssl1.0 \
  && rm -rf /var/lib/apt/lists/*

RUN curl -LsS https://aka.ms/InstallAzureCLIDeb | bash \
  && rm -rf /var/lib/apt/lists/*

#############################################################################
## END local#################################################################
#############################################################################

# "downloadUrl": "https://vstsagentpackage.azureedge.net/agent/<AGENT_VERSION>/vsts-agent-linux-x64-<AGENT_VERSION>.tar.gz"
# Note: This was in the old instructions but not in the new one.
# 1) To get the latest version, start with a known good (3.212.0) then let
#    Azure update it (which might not work, but it shows the version to use).
# 2) For some reason the azure install does a chown (which I override).
# 3) Not clear whether TARGETARCH should be amd64 or linux-x64 or x64.
ARG TARGETARCH=amd64
#ARG AGENT_VERSION=2.185.1
ARG AGENT_VERSION=3.241.0
WORKDIR /azp
RUN if [ "$TARGETARCH" = "amd64" ]; then \
      AZP_AGENTPACKAGE_URL=https://vstsagentpackage.azureedge.net/agent/${AGENT_VERSION}/vsts-agent-linux-x64-${AGENT_VERSION}.tar.gz; \
    else \
      AZP_AGENTPACKAGE_URL=https://vstsagentpackage.azureedge.net/agent/${AGENT_VERSION}/vsts-agent-linux-${TARGETARCH}-${AGENT_VERSION}.tar.gz; \
    fi; \
    curl -LsS "$AZP_AGENTPACKAGE_URL" | tar -xz

RUN chown -R agent:agent /azp /home/agent
USER agent

# Another option is to run the agent as root.
# ENV AGENT_ALLOW_RUNASROOT="true"

ENTRYPOINT [ "./start.sh" ]

# building:
#   docker build -t dockeragent:latest .
# running:
#   Create a PAT with "Agent pools read/write" scope and maybe "Code read",
#   replace __REDACTED__ with that token. Apparently you don't need username.
#   docker run -d --restart=always --name azure-slb1-swt -e AZP_URL=https://dev.azure.com/slb1-swt -e AZP_TOKEN=__REDACTED__ -e AZP_AGENT_NAME=`uname -n` -e AZP_POOL=Paal-Kvamme-Test -v /var/run/docker.sock:/var/run/docker.sock dockeragent:latest
