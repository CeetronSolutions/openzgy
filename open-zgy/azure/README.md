## Contents of azure/

Files used for building with Azure DevOps.

On Linux, most pipelines run the build and the tests tests inside Docker.
A container named ${TAG} will be explicitly created by the build steps.
Not by turning in some "docker" option in Azure. ${TAG} is different
in each target but is the same in each build of one specific target.
This is done to reduce the need for manually cleaning up old images.

${TAG}:latest and ${TAG}:old can lead to a few race conditions if
there are parallel builds of the same target.

The first race is harmless. The tags are used to clean up old images,
while making sure that the running build can still find the previous
build in the docker cache. A collision might cause a build to take
longer time and/or leave an unnamed image that eventually needs manual
cleanup.

The second race applies to pipelines that build and test in separate
steps. Tests pick up the docker container by its tag. So there a risk
that tests might run for the wrong build. N/A if there is just one step.

The second race can be avoided by somehow using Azure Devops features
to pass on the container id to the test steps. Maybe someday do this.
